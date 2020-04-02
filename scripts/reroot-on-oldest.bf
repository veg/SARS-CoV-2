LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/convenience/regexp.bf");

KeywordArgument ("tree", "A phylogenetic tree to reroot");
KeywordArgument ("csv", "The CSV with sequence ID to date mapping (last column)");
KeywordArgument ("output", "Write the rerooted tree to");

root.tree = trees.LoadAnnotatedTopology (TRUE);
root.data = io.ReadDelimitedFile (null, ",", TRUE);

root.tree.string = root.tree[terms.trees.newick_with_lengths];

Topology T = root.tree.string;

id_to_date = {};

for (r; in; root.data["rows"]) {
    date = r[Abs(r)-1];
    if (Abs (date) == 8) {
        id_to_date[r[0]] = date;
    }
}

root.oldest          = "9999999";
root.oldest_sequence = "";

for (seq; in; T) {
    seq_s = regexp.Replace (seq, "_[0-9]+$", "");
    seq_s = regexp.Find (seq_s,"^[^_.]+_[^_.]+_[^_]+");
    if (id_to_date/seq_s) {
        if (id_to_date[seq_s] < root.oldest) {
            root.oldest_sequence = seq;
            root.oldest = id_to_date[seq_s];
        }
    }
}

root.output = io.PromptUserForFilePath ("Write the rerooted tree to");
fprintf (root.output, CLEAR_FILE, ((trees.RootTree (root.tree , root.oldest_sequence))["tree"])[terms.trees.newick_with_lengths])


