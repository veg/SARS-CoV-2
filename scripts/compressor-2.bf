RequireVersion ("2.5.19");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary     ("libv3/IOFunctions.bf");
LoadFunctionLibrary     ("libv3/convenience/math.bf");



filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment with duplicate sequences (JSON), a CSV with sequence variants by site, and a JSON with sequence variants by sequence.
                            Filter minority variants (replace with gaps), and remove sequences with outlier-like variant count behavior.
                            ",
                            terms.io.version :          "0.1",
                            terms.io.reference :        "TBD",
                            terms.io.authors :          "Sergei L Kosakovsky Pond",
                            terms.io.contact :          "spond@temple.edu",
                            terms.io.requirements :     "An MSA and, optionally, a tree"
                          };


io.DisplayAnalysisBanner (filter.analysis_description);

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", FALSE);


KeywordArgument                     ("msa", "The MSA to filter rare variants from");

DataSet filter.dataset              = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter filter.input          = CreateFilter (filter.dataset,1);

console.log ("> Loaded an alignment with `filter.input.species` sequences and `filter.input.sites` sites from `LAST_FILE_PATH`");

KeywordArgument     ("duplicates", "Load sequence duplicate information from here", None);
filter.dups = io.PromptUserForString ("Load sequence duplicate information from here");
fscanf (filter.dups, "Raw", filter.dup_data);
filter.dup_data = Eval (filter.dup_data);

filter.total = 0;

for (k, v; in; filter.dup_data) {
    filter.total += Abs (v);
}

console.log ("> Total # of sequences, counting duplicates = " + filter.total);

KeywordArgument     ("csv", "CSV with variants", None);
filter.variants = io.ReadDelimitedFile (null, ",", TRUE);

/*

"header":{
   "0":"Site",
   "1":"Consensus",
   "2":"A",
   "3":"G",
   "4":"C",
   "5":"T",
   "6":"ambig",
   "7":"N",
   "8":"gap"
  }
  
"1005":
{
 "0":"298",
 "1":"T",
 "2":"0",
 "3":"7",
 "4":"0",
 "5":"75173",
 "6":"1",
 "7":"0",
 "8":"0"
}
*/


KeywordArgument     ("byseq", "Load the .JSON with by sequence information", None);
filter.byseq = io.PromptUserForString ("Load the .JSON with by sequence information");
fscanf (filter.byseq, "Raw", filter.byseq_data);
filter.byseq_data = Eval (filter.byseq_data);
KeywordArgument     ("p", "Binomial probability cutoff", 0.9);
filter.p = io.PromptUser ("> Binomial probability cutoff",0.9,0,1,FALSE);
filter.by_site_cutoff = filter.binomial_cutoff (filter.total, 1./10000, filter.p);


variant_count   = {filter.total , 1};
variant_count_by_seq = {};
i = 0;

for (k, v; in; filter.byseq_data) {
    c = Abs (filter.dup_data[k]);
    vc = 0;
    for (k2, cc; in; v) {
        if (k2 != 'copies') {
            if (+cc <= filter.by_site_cutoff) {
                vc+=1;
            }
        }
    }
    variant_count_by_seq[k] = vc;
    for (j = 0; j < c; j+=1) {
        variant_count[i] = vc;
        i+=1;
    }
}


filter.by_seq_cutoff = math.GatherDescriptiveStats (variant_count);
filter.by_seq_cutoff = (filter.by_seq_cutoff["mean"] + filter.by_seq_cutoff["Std.Dev"] * 5 + 0.5)$1;


console.log ("> Setting sequence minority variant cutoff to " + filter.by_seq_cutoff);
console.log ("> Setting minority variant cutoff to " + filter.by_site_cutoff);

// record all pairs of sites that appear to have minority variants in the same sequence

filter.all_pairs = {};
filter.pairs_by_seq = {};

for (seq = 0; seq < filter.input.species; seq += 1) {
    GetString (seq_name, filter.input, seq);
    filter.vc = filter.byseq_data[seq_name];   
    filter.mv = {};
    
    for (v, c; in; filter.vc) {
        if (v != 'copies') {
            if (+c < filter.by_site_cutoff) {
                 filter.mv [+v] = filter.vc['copies'];
            }
        }    
    }
    
    n = Abs (filter.mv);
    if (n > 1) {
        filter.paired = {n,2};
        i = 0;
        for (v, c; in; filter.mv) {
            filter.paired[i][0] = +v;
            filter.paired[i][1] = +c;
            i+=1;
        }
        
        filter.paired = filter.paired % 0;
        filter.pairs_by_seq [seq_name] = filter.paired;
        
        for (i = 0; i < n; i += 1) {
            for (j = i + 1; j < n; j += 1) {
                if (filter.paired[j][0] - filter.paired[i][0] >= 3) {
                    key = "" + filter.paired[i][0] + "|" + filter.paired[j][0];
                    filter.all_pairs [key] += filter.paired[j][1];
                }
            }
        }   
    
    }
}



filter.all_duplicates = {};

filter.fasta_string = ""; filter.fasta_string * 256000;
filter.input_seqs = {};
filter.edits = {};

for (seq = 0; seq < filter.input.species; seq += 1) {
    GetString (seq_name, filter.input, seq);
    GetDataInfo (seq_chars, filter.input, seq);
    filter.input_seqs [seq_name] = seq_chars;
    filter.vc = variant_count_by_seq[seq_name]; 
    
    filter.protected_sites = {};
    
    if (filter.pairs_by_seq / seq_name) {
        filter.paired = filter.pairs_by_seq [seq_name];
         for (i = 0; i < n; i += 1) {
            for (j = i + 1; j < n; j += 1) {
                key = "" + filter.paired[i][0] + "|" + filter.paired[j][0];
                
                if (filter.all_pairs[key] > 1) {
                    //console.log (key);
                    //console.log (filter.all_pairs[key]);
                    filter.protected_sites [filter.paired[i][0]] = 1;
                    filter.protected_sites [filter.paired[j][0]] = 1;
                }
            }
        }           
     }
    
    if (Abs (filter.protected_sites)) {
        console.log (">Found `Abs(filter.protected_sites)` protected sites (covariation): " +  Join (",", Rows(filter.protected_sites)));
    }
      
    if (filter.vc > filter.by_seq_cutoff && Abs (filter.protected_sites) == 0) {
        console.log (">Removing `seq_name` because is has " + filter.vc + " minority variants (variable outlier)"); 
        filter.edits [seq_name] = "removed";
        //filter.all_duplicates - seq_name;
    } else {
        filter.all_duplicates[seq_name] = filter.dup_data [seq_name];
        filter.to_filter = {};
        filter.edits [seq_name] = {};
        filter.has_protected_pairs = 0;
        
        for (v, c; in; filter.byseq_data[seq_name]) {
            if (v != 'copies') {
                if (+c < filter.by_site_cutoff) {
                    if (filter.protected_sites[+v] == 0) {
                        filter.to_filter [+v] = +c;
                    } 
                }
            }
        }
        
        if (Abs ( filter.to_filter)) {
            sl = Abs(seq_chars);
            filter.filtered_string = "";  filter.filtered_string * sl; 
            filter.punched = "";
            for (i = 0; i < sl; i+=1) {
                if (filter.to_filter / i) {
                    filter.filtered_string * "-";
                    filter.punched += seq_chars[i];
                    (filter.edits [seq_name])[i] = seq_chars[i];
                } else {
                    filter.filtered_string * seq_chars[i];
                }
            }
            filter.filtered_string * 0;
            console.log (">Punching `Abs(filter.to_filter)` (`filter.punched`) positions from `seq_name`: " + filter.punched); 
            filter.fasta_string * ">`seq_name`\n";
            filter.fasta_string * filter.filtered_string;
            filter.fasta_string * "\n";
        } else {
            filter.fasta_string * ">`seq_name`\n";
            filter.fasta_string * seq_chars;
            filter.fasta_string * "\n";
        }
    }
 }

filter.fasta_string * 0;

DataSet filtered.reduced = ReadFromString (filter.fasta_string);
DataSetFilter filtered.reducedF = CreateFilter (filtered.reduced,1);

filter.reduced_name_map = {};
for (i = 0; i < filtered.reducedF.species; i+=1){
    GetString (seq_name, filtered.reducedF, i);
    filter.reduced_name_map [i] = seq_name;
}


console.log (">Compressing error-corrected sequences");
GetDataInfo (filter.duplicate_info, filtered.reducedF, -2);

console.log (">Retained `filter.duplicate_info['UNIQUE_SEQUENCES']` unique haplotypes post-filtering");


filter.uniqueID2Name       = {};
filter.uniqueID2compressed = {};

filter.group_by_id = {};

for (i, j; in; filter.duplicate_info["SEQUENCE_MAP"]) {
    if (filter.uniqueID2Name / j == FALSE) {
        //console.log ("STORE `j` => `(filter.duplicate_info["UNIQUE_INDICES"])[j]`");
        filter.uniqueID2Name[j] = filter.reduced_name_map[(filter.duplicate_info["UNIQUE_INDICES"])[j]];
        filter.group_by_id [j] = {};
    }  
    filter.unique_ref = (filter.duplicate_info["UNIQUE_INDICES"])[j];
    filter.group_by_id [j] + (">`filter.reduced_name_map[i]`\n" + filter.input_seqs [filter.reduced_name_map[i]]);
    //console.log (i);
    if ( filter.reduced_name_map[i] != filter.uniqueID2Name[j]) {
        //console.log ("Merging `(filter.reduced_name_map[i])` into `filter.uniqueID2Name[j]`");
        for (seq; in; filter.all_duplicates[(filter.reduced_name_map[i])]) {
            filter.all_duplicates[filter.uniqueID2Name[j]] + seq;
        }
        filter.all_duplicates - (filter.reduced_name_map[i]);
    } 
}

DataSetFilter filtered.reducedF2 = CreateFilter (filtered.reducedF, 1, "", Join (",", filter.duplicate_info["UNIQUE_INDICES"]));

KeywordArgument     ("output", ".fasta for compressed data", None);
filter.out = io.PromptUserForFilePath(".fasta for compressed data");

KeywordArgument     ("output-edits", ".json for performed edit operations", None);
filter.edit_out = io.PromptUserForFilePath(".json for performed edit operations");
fprintf (filter.edit_out, CLEAR_FILE, filter.edits);

utility.SetEnvVariable ("DATA_FILE_PRINT_FORMAT", 9);
fprintf (filter.out, CLEAR_FILE, filtered.reducedF2);

/*
for (i, j; in; filter.group_by_id) {
    fprintf (filter.out + "_" + i, CLEAR_FILE, Join ("\n",j));
}
*/

KeywordArgument     ("json", ".json for updated duplicate counts", None);
filter.json = io.PromptUserForFilePath(".json for updated duplicate counts");

fprintf (filter.json, CLEAR_FILE,  filter.all_duplicates);


lfunction filter.binomial_cutoff (N,p,t) {
    s = 0;
    term = (1-p)^N;
    i = 0;
    while (s < t) {
        s += term;
        term = term * p/(1-p)*(N-i)/(i+1);
        i+=1;
    }
    return i;
}


