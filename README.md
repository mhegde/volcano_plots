<b>Updates</b>: v2.3.2: Requires pandas-0.17.0 or higher

Requires:
Python 2.7;
pandas;
scipy;
numpy

To run this code type:

python hypergeom-2.3.2.py --input-file <Path to inputfile> --chip-file <Path to chip file> --sing-pert <Y/N> --fraction <1-100> --ctrls-num <> --min-pert <> --max-pert <> --label-num <>


Code arguments:
1. Input file, .txt file with construct barcodes in first column, score for each sample in the following columns
2. Chip file, File mapping construct barcodes to genes, construct barcodes in first column, gene target in second column
3. Option to include/exclude singletons in output file;Default: Singletons included
4. Option to specify n% for mean p-val, LFC calculation, Default: 100
5. Option of number of control guides to be included in a random set, Default: 4
6. Option of minimum number of perturbations for a gene to be included in the volcano plot, Default: 3
7. Option of maximum number of perturbations for a gene to be included in the volcano plot, Default: 8
8. Option of number of genes to be labeled on the plot, Default: 3
