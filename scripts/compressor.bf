RequireVersion ("2.5.19");

LoadFunctionLibrary     ("libv3/tasks/alignments.bf");
LoadFunctionLibrary     ("libv3/tasks/trees.bf");
LoadFunctionLibrary     ("libv3/UtilityFunctions.bf");

filter.analysis_description = {terms.io.info :
                            "
                            Read an alignment with duplicate sequences (JSON), create a CSV with sequence variants by site, and a JSON with sequence variants by sequence.
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

for (v; in; filter.dup_data) {
    filter.total += Abs (v);
}

console.log ("> Total # of sequences, counting duplicates = " + filter.total);
filter.variants_by_site  = {};
filter.consensus_variant = {};
filter.duplicates_by_seq = {};
filter.all_duplicates = {};

GetDataInfo (filter.site_patterns, filter.input);
filter.unique_patterns = utility.Array1D (filter.input.site_freqs);

filter.seq_names = {}; // id => name

for (seq = 0; seq < filter.input.species; seq += 1) {
    GetString (seq_name, filter.input, seq);
    filter.duplicates_by_seq [seq] = Abs (filter.dup_data[seq_name]);
    filter.all_duplicates [seq] = filter.dup_data[seq_name];
    filter.seq_names [seq] = seq_name;
}

KeywordArgument     ("output", ".CSV for nucleotide counts", None);
filter.out = io.PromptUserForFilePath(".CSV for nucleotide counts");
fprintf (filter.out, KEEP_OPEN, "Site,Consensus,A,G,C,T,ambig,N,gap\n");

KeywordArgument     ("json", ".json for by-sequence variants", None);
filter.json = io.PromptUserForFilePath(".json for by-sequence variants");


filter.patter2site = {};

for (i,j,v; in; filter.site_patterns) {
    index = i+j;
    if (filter.patter2site / v == FALSE ) {
        filter.patter2site [v] = {};
    }  
    filter.patter2site [v] + index;
}


GET_DATA_INFO_RETURNS_ONLY_THE_INDEX = TRUE;
COUNT_GAPS_IN_FREQUENCIES = FALSE;

variantsBySequence = {};
/* 
id=> {"pattern" : count,...}
*/

for (seq = 0; seq < filter.input.species; seq += 1) {
   variantsBySequence[filter.seq_names[seq]] = {'copies' : filter.duplicates_by_seq[seq]};
}

//#profile START;

filter.letters = "ACGTNN-";

for (pattern = 0; pattern < filter.unique_patterns; pattern += 1) {
    io.ReportProgressBar ("filter","Processing pattern " + pattern);
      
    counts_by_nucs = {7,1};
    position_by_sequence = {};
    by_resolution = {};
    for (seq = 0; seq < filter.input.species; seq += 1) {
        GetDataInfo (pattern_info, filter.input, seq, pattern);  
        if (pattern_info < 0) {
            GET_DATA_INFO_RETURNS_ONLY_THE_INDEX = FALSE;
            GetDataInfo (pattern_info, filter.input, seq, pattern);  
            if (+pattern_info == 0) { index = 6;} else {
                if (+pattern_info == 4) { index = 5; } else { index = 4; }
            }
            GET_DATA_INFO_RETURNS_ONLY_THE_INDEX = TRUE;
        } else {
            index = pattern_info;
        }
        by_resolution[index] = by_resolution[index] + filter.duplicates_by_seq[seq];
        position_by_sequence [seq] = index;
        counts_by_nucs [index] += filter.duplicates_by_seq[seq];
    }
    consensus = +(Max (by_resolution,1))["key"]; 
    /*
        { "key":"0","value":75174}
    */

    for (s,r; in; position_by_sequence) {
       if (r >= 0 && r < 6 && r != consensus) {
         for (s2; in; filter.patter2site [pattern]) {
            (variantsBySequence [filter.seq_names[s]])[s2] = by_resolution[r];
         }
       }
    }
    
    for (s; in; filter.patter2site [pattern]) {
        fprintf (filter.out, s, ",", filter.letters[consensus],",", Join (",",counts_by_nucs), "\n");
    }
}


USE_JSON_FOR_MATRIX = 1;
fprintf (filter.json, CLEAR_FILE, variantsBySequence);

//utility.FinishAndPrintProfile ();


