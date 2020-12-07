"""Import SARS-CoV-2 sequence batches into the master JSON sequence record

This utility script reads in JSON files (scraped from GISAID or otherwise suitably prepared), 
parses sequence headers using attributed regular expressions and converts them to an
N-field sequence record (see the `expected_fields` for the currently supported fields). Records are
keyed on unique ACCESSION NUMBERS (provided by GISAID )

an example GISAID record looks like 

"ACCESSION" : {
        "Virus name": "hCoV-19/USA/WA-UW77/2020",
        "Accession ID": "EPI_ISL_416433",
        "Type": "betacoronavirus",
        "Passage details/history": "Original",
        "Collection date": "2020-03-10",
        "Location": "North America / USA / Washington",
        "Host": "Human",
        "Gender": "unknown",
        "Patient age": "unknown",
        "Sequencing technology": "Illumina NextSeq",
        "Assembly method": "custom",
        "Coverage": "19867x",
        "Originating lab": "UW Virology Lab",
        "Address": "1100 Fairview Ave N,\n98109 Seattle",
        "Submitting lab": "UW Virology Lab",
        "Authors": "Pavitra Roychoudhury, Hong Xie, Keith Jerome, Alexander Greninger",
        "Submitter": "Roychoudhury, Pavitra",
        "Submission Date": "2020-03-21",
        "FASTA" : sequence
}

Some fields may be missing or unknown;


and storing sequences with ACCESSION NUMBERS to a flat JSON file, in the format::

    json [ACCESSION] = {'descriptor' :
                                'attributes' : {'attr1' : ... 
                                                'attr2' : ..., }
                                'sequence' : UPPER case nucleotide sequence
                          }

if a JSON file is to be updated (via the --update flag), an archived copy of
its previous contents will be written to the same directory. The script outputs
diagnostics to ``stderr`` and summary statments in Markdown table rows to ``stdout``
or an existing log file::

    | `batches/sequences.json` | November 11 2015 (15:22) | **0** new sequences added. **1** duplicate sequences. **0** sequences errored | |

Author:
    Sergei L Kosakovsky Pond (spond@temple.edu)

Version:
    v0.0.1 (2020-03-22)


"""

import argparse
import sys
import json
import re
import datetime
import os
from   os import  path
import bz2
import tempfile

#-------------------------------------------------------------------------------
def process_fasta (fasta_in):
    fasta_in = fasta_in.splitlines ()
    return "".join ([k.strip() for k in fasta_in[1:]]).upper ()
    
def process_coverage (coverage):
    try:
        m = re.search('([0-9]+)x', coverage)
        return int (m.group(1))
    except:
        return coverage

def process_date (fmt):
    def date_handler (date):
        try:      
            return datetime.datetime.strptime (date, fmt).strftime("%Y%m%d")
        except:
            return date
    return date_handler

    
def process_location (location):
    pieces = [k.strip() for k in location.split ('/')]
    location = {
        'subregion' : pieces[0],
        'country' : None,
        'state' : None,
        'locality' : None
    }
    try:
        location['country'] = pieces[1]
        location['state'] = pieces[2]
        location['locality'] = pieces[3]
    
    except:
        pass
    return location
    
    
#-------------------------------------------------------------------------------

def update_json (json_object, file_name):

    file_directory = path.dirname(path.abspath(file_name))

    if path.isfile (file_name): # need to archive and compress
        with open (file_name, 'r') as existing:
            old_json = json.load (existing)
            old_json = json.dumps (old_json)

        print ("Compressing existing JSON object to a temporary file", file = sys.stderr)
        temp_file = tempfile.mkstemp(suffix=".bz2", prefix="pirc_db", dir=file_directory)[0]
        os.write (temp_file, bz2.compress (bytes(old_json, 'UTF-8')))
        os.close(temp_file)


    with open (file_name, 'w') as json_dump:
        json.dump (json_object, json_dump, sort_keys=True, indent=1)
    

#-------------------------------------------------------------------------------


arguments = argparse.ArgumentParser(description='Import FASTA sequences into the PIRC sequence database.')

arguments.add_argument('-f', '--file', help = 'Read sequences from this JSON file', required = True, type = argparse.FileType('r'), nargs = '*')
arguments.add_argument('-j', '--json', help = 'The master JSON cache file where sequence records live', required = True, type = str)

arguments.add_argument('-t', '--dformat',  help = 'The strptime for the date field (e.g. %%Y%%m%%D)', required = False, default = '%Y-%m-%d')
arguments.add_argument('-u', '--update', help = 'Parse the sequences and update the JSON master and logs (defaults to dry run)', default = False, action = 'store_true')
arguments.add_argument('-C', '--Comment', help = 'Optional comment text to add to the summary output', type = str, required = False, default = '')
arguments.add_argument('-L', '--log',  help = 'Append a transaction log line to this file', required = False, type = argparse.FileType('r+'))

import_settings = arguments.parse_args()
default_values  = vars (import_settings)


'''
These are the fields that will be imported from the JSON
    [
        id in database json,
        id in import json,
        T/F : is required
        None | lambda : default transform
    ]
'''



expected_fields = [['id', 'Accession ID',   True, lambda x : x.lower()],
                   ['name', 'Virus name', True, None],
                   ['sequence', 'FASTA',True, process_fasta],
                   ['collected', 'Collection date', False, process_date (default_values["dformat"])],
                   ['submitted', 'Submission Date', False, process_date (default_values["dformat"])],
                   ['type', 'Type', False, None],
                   ['passage', 'Passage details/history', False, None],
                   ['host', 'Host', False, None],
                   ['gender','Gender', False, None],
                   ['age','Patient age', False, None],
                   ['technology','Sequencing technology', False, None],
                   ['assembly','Assembly method', False, None],
                   ['coverage', 'Coverage', False, process_coverage],
                   ['lab','Originating lab', False, None],
                   ['authors','Authors', False, None],
                   ['submitter','Submitter', False, None],
                   ['address','Address', False, None],
                   ['location', 'Location', False, process_location]
                   ]


current_sequence_db = {}
try:
    with open (import_settings.json, 'r') as cache:
        try:
            current_sequence_db = json.load (cache)
        except ValueError:
               print ('Error parsing JSON; please fix the file format and try again', file = sys.stderr)
               sys.exit (1)
except:
    pass

failed_ids = set()

successfully_imported   = 0
already_in_db           = 0
json_changed            = False

for input_file in import_settings.file:
    for seq_id, match_dict in json.load(input_file).items():

        components = {}

        try:
            for field in expected_fields:
                if field[1] in match_dict:
                   components [field[0]] = field[3] ( match_dict[field[1]]) if field[3] else match_dict[field[1]]
                else:
                    if field[2]:
                        print ('Required field "%s" was not extracted from %s' % (field[1], seq_id), file = sys.stderr)
                        failed_ids.add (seq_id)
                        break
                    else:
                        if field[0] in default_values:
                            components [field[0]] = field[3] ( default_values[field[0]]) if field[3] else default_values[field[0]]
                        else:
                            components [field[0]] = None



        except Exception as e:
           print ('Sequence %s generated a processing error(%s)' % (seq_id, e), file = sys.stderr)
           failed_ids.add (seq_id)

        except:
            raise

        if seq_id not in failed_ids:

            seq_index = components['id']
            seq_label = seq_index
            if seq_index not in current_sequence_db:
                current_sequence_db [seq_index] = {}
                json_changed = True
            
            seq_data = components['sequence']

            if 'sequence' in current_sequence_db[seq_index]:
                if seq_data == current_sequence_db[seq_index]['sequence']:
                    already_in_db += 1
                else:
                    #failed_ids.add (seq_id)
                    print ("Sequence %s has a duplicate tag, but different nucleotides than what's already in the database "% (seq_id), file = sys.stderr)     
                    current_sequence_db[seq_index]['sequence'] = seq_data
            else:
                current_sequence_db [seq_index] = components
                json_changed = True
                successfully_imported += 1
                print ("Imported %s (%d nucs)" % (seq_id, len (seq_data)), file = sys.stderr)

if import_settings.log:
    import_settings.log.seek (0,2)

if json_changed and import_settings.update:
    update_json (current_sequence_db, import_settings.json)

print ("| `%s` | %s | % s | % s|" % (",".join([k.name for k in import_settings.file]), datetime.datetime.now().strftime ("%B %d %Y (%H:%M)"),
       "**%d** new sequences added. **%d** duplicate sequences. **%d** sequences errored" % (successfully_imported, already_in_db, len (failed_ids)), import_settings.Comment),
       file = import_settings.log if import_settings.log and import_settings.update else sys.stdout)

