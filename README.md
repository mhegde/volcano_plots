To run this code type:

python hypergeom.py --input-file <Path to inputfile> --chip-file <Path to chip file> --sing-pert <Y/N> --fraction <1-100> --ctrls-num <> --min-pert <> --max-pert <> --label-num <> --disp-nosite <Y/N> --disp-onesite <Y/N>


Code arguments:
1. Input file, .txt file with construct barcodes in first column, score for each sample in the following columns
2. Chip file, File mapping construct barcodes to genes, construct barcodes in first column, gene target in second column
3. Option to include/exclude singletons in output file;Default: Singletons included
4. Option to specify n% for mean p-val, LFC calculation, Default: 100
5. Option of number of control guides to be included in a random set. This number is common for both types of controls, NO_SITE and ONE_NON-GENE_SITE. Default: 4
6. Option of minimum number of perturbations for a gene to be included in the volcano plot, Default: 3
7. Option of maximum number of perturbations for a gene to be included in the volcano plot, Default: 8
8. Option of number of genes to be labeled on the plot, Default: 3
9. Option to determine whether no-site controls will be displayed on the resulting volcano plot
10.Option to determine whether one-non-gene-site controls will be displayed on the resulting volcano plot 
