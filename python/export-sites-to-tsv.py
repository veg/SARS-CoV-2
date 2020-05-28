import sys
import json
import argparse
import operator
import csv

arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-f', '--file',   help = 'File to process', required = True, type = str)
import_settings = arguments.parse_args()

'''
"26522": {
  "G": "M",
  "S": 1,
  "bCFEL": {
   "p": 1,
   "a": 0,
   "b-nCOV": 0,
   "b": 0
  },
  "bFEL": {
   "a": 0,
   "b": 0,
   "p": 1
  },
  "bMEME": {
   "p": 1,
   "a": 0,
   "b+": 0,
   "w+": 0,
   "b-": 0,
   "w-": 1,
   "br": 0
  },
  "bSC2": "ATG",
  "bSC2-aa": "M",
  "bcdn": {
   "nCOV": {
    "ATG": 11
   },
   "others": {
    "ATG": 58
   }
  },
  "baa": {
   "nCOV": {
    "M": 11
   },
   "others": {
    "M": 58
   }
  },
  "evo": null,
  "cdn": {
   "ATG": 30821,
   "ANG": 1
  },
  "aa": {
   "M": 30821
  },
  "SLAC": {
   "N": 0,
   "S": 0,
   "EN": 3,
   "ES": 0,
   "p": 0,
   "NmS": null
  },
  "FEL": {
   "a": 0.1786478280845536,
   "b": 0,
   "p": 1
  },
  "MEME": {
   "p": 0.6666666666666666,
   "a": 0,
   "b+": 0,
   "w+": 0.25,
   "b-": 0,
   "w-": 0.75,
   "br": 0
  },
  "fubar": {
   "a": 4.08086035878136,
   "b": 0.4081288875520895,
   "p+": 0.2780466402160124,
   "p-": 0.6966442027975785,
   "w": 1.1647232610679439,
   "wl": 0.0,
   "wu": 26.0
  },
  "trend": 0.23682222678687434,
  "subs": {
   "cdn": {},
   "aa": {},
   "lcdn": {
    "ANG|ATG": 1
   },
   "laa": {
    "KTRM|M": 1
   }
  }
 }
'''

null = "null"
nucs = set (["A","C","G","T"])
pv = 0.1

writer = csv.writer (sys.stdout, delimiter = "\t", )

writer.writerow (["genomic_coordinate","gene",
                  "codon_in_gene","consensus","consensus_aa",
                  "genomes","syn_sites",
                  "nonsyn_sites","variants","aa_variants",
                  "selection","epidosic_selection","multiple_clades","trend","epitopes","related_selection","predicted_variants","unusual_variants","score"])

with open (import_settings.file, "r") as fh:
    data = json.load (fh)
    for site in sorted (data.keys(), key = lambda x: int (x)):
        info = data[site]
        if "aa" in info:
            record = ["%d" % (int(site)+1),info["G"],info["S"],
                        max(info["cdn"].items(), key=operator.itemgetter(1))[0],
                        max(info["aa"].items(), key=operator.itemgetter(1))[0], "%d" % sum (info["cdn"].values())]
            if "SLAC" in info:
                record.extend (["%g" % info["SLAC"]["ES"], "%g" % info["SLAC"]["EN"]])
            else:
                record.extend ([null,null])

            record.append ("%s" % "|".join (["%s:%s" % (k[0], k[1]) for k in sorted(info["cdn"].items(), key=operator.itemgetter(1), reverse = True)]))
            record.append ("%s" % "|".join (["%s:%s" % (k[0], k[1]) for k in sorted(info["aa"].items(), key=operator.itemgetter(1), reverse = True)]))
               
            score = 0   
               
            if "FEL" in info:
                if info["FEL"]["p"] <= pv:
                    record.append ("1" if info["FEL"]["a"] < info["FEL"]["b"] else "-1")
                    score += 1
                else:
                    record.append ("0")
            else:
                record.append (null)
               

            if "MEME" in info:
                if info["MEME"]["p"] <= pv:
                    record.append ("1")
                    score += 0 if score > 0 else 1
                    if info["MEME"]["br"] > 1:
                        score += 1
                        
                    record.append ("%d" % info["MEME"]["br"])
                else:
                    record.append ("0")
            else:
                record.append (null)


            if "trend" in info:
                if info["trend"] > 3:
                    record.append ("1")
                    score += 1
                else:
                    record.append ("0")
            else:
                record.append (null)

            if "hla" in info:
                record.append ("%d" % len (info["hla"]))
                score += 1 if len (info["hla"]) else 0
            else:
                record.append (null)

            if "bFEL" in info:
                if info["bFEL"]["p"] <= pv:
                    record.append ("1" if info["bFEL"]["a"] < info["FEL"]["b"] else "-1")
                    score += 1
                else:
                    record.append ("0")
            else:
                record.append (null)

            if "evo" in info and info["evo"]:
                record.append ("%s" % ":".join (info["evo"].keys()))
                unusual_variants = [k for k in info["cdn"] if (len ([n for n in k if n not in nucs]) == 0 and k not in info["evo"])]
                record.append ("%s" % ":".join (unusual_variants))
                if len (unusual_variants):
                    score += 1
            else:
                record.extend ([null,null])
                
            record.append ("%d" % score)

            writer.writerow (record)
        
        
 