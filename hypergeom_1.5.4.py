'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
Performs hypergeometric analysis for genetic perturbation screens
Input: 1. Input file with column headers
       2. Chip File
       3. Option to include/exclude singletons in output file;Default: Singletons included
       4. Option to specify n% for mean p-val, LFC calculation
'''
import pandas as pd
import numpy as np
from scipy import stats
from math import log10
import csv, argparse, os, sys, re
from datetime import datetime
from decimal import Decimal

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file',
        type=str,
        help='File containing sgRNA sequence and score(s) columns with headers')
    parser.add_argument('--chip-file',
        type=str,
        help='Chip file')
    parser.add_argument('--sing-pert',
        type=str,
        default='Y',
        help='Y to include genes with single perturbations in output file, N to exclude; Default: Y')
    parser.add_argument('--fraction',
        type=int,
        default=100,
        help='Fraction of perturbations to be included in mean LFC, pval calculation, 1 - 100; Default: 100')
    return parser

'''
Sorts a data frame on the column specified and re-indexes
Argument: Dataframe, column to sort on, direction of sorting
Return Value: Sorted and re-index dataframe
'''
def sort_reindex(df, col, direction):
    if direction == 'P' or direction == 'p':
        df = df.sort(columns=col, ascending=False)
        df.index = range(0, len(df))
        #df['Rank'] = map(lambda x:x+1, range(len(df)))
    elif direction == 'N' or direction == 'n':
        df = df.sort(columns=col, ascending=True)
        df.index = range(0, len(df))
    else:
        print 'Please enter a relevant direction; P for positive and N for negative'
        sys.exit(1)
    return df

def calc_hypergeom_scores(merged, st_in, ge, frac):
    grouped = merged.groupby('Gene Symbol')
    sps = list(st_in['Guide Sequence'])
    tot_sps = len(sps)
    sp_rank = {}
    for i in range(len(sps)):
        sp_rank[sps[i]] = i+1
    g_p_val = {}
    for g in ge:
        ge_df = grouped.get_group(g)
        ge_df = ge_df.drop_duplicates('Guide Sequence')
        gene = g
        guides = list(ge_df['Guide Sequence'])
        guide_list = ';'.join(guides)
        length_sps = len(guides)
        rank = list()
        rank_hash = {}
        for s in guides:
            r=sp_rank[s]
            rank.append(r)
            rank_hash[s] = r
        rank.sort()
        lfcs = list(ge_df['Score'])
        op_lfcs = ';'.join(str(round(l,2)) for l in lfcs)
        #op_lfcs = ';'.join([str(l) for l in lfcs])
        all_ranks = ';'.join(str(r) for r in rank)
        p_values = [-log10(stats.hypergeom.pmf(rank.index(x)+1, tot_sps, length_sps, x)) for x in rank]
        p_values.sort(reverse=True)
        num_perts = int(frac/100.0*length_sps)
        avg_p_val = np.mean(p_values[0:num_perts])
        avg_lfc = np.mean(lfcs[0:num_perts])
        op_p_vals = ';'.join("%.2E" %Decimal(p) for p in p_values)
        g_p_val[g] = op_lfcs+'_'+str(avg_lfc)+'_'+op_p_vals+'_'+str(avg_p_val)+'_'+guide_list+'_'+all_ranks
    return g_p_val


if __name__ == '__main__':
    args = get_parser().parse_args()
    inputfile = args.input_file
    input_df = pd.read_table(inputfile)
    cols = list(input_df.columns)[1:]
    cols = [re.sub('[^a-zA-Z0-9 \n\.]', '_', x) for x in cols]
    cols.insert(0,'Guide Sequence')
    input_df.columns = cols
    cols_iter = cols[1:]
    ref = pd.read_table(args.chip_file)
    ref_colnames = list(ref.columns)
    ref_colnames[0:2] = ['Guide Sequence', 'Gene Symbol']
    ref.columns = ref_colnames
    include = args.sing_pert
    frac = args.fraction
    inputname = inputfile.split('/')[-1]
    o_folder = 'Volcano-plots_'+inputname[:-4]+'_'+str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
    if not os.path.exists(o_folder):
        os.makedirs(o_folder)
    for ci,c in enumerate(cols_iter):
        outputfile = o_folder+'/'+c+'_'+str(frac)+'_'+str(datetime.now().strftime("%y-%m-%d-%H-%M"))+'.txt'
        st_in = input_df[['Guide Sequence',c]]
        st_in = st_in.rename(columns={c:'Score'})
        merged = pd.merge(st_in, ref, on='Guide Sequence')
        ge=list(merged['Gene Symbol'])
        ge=list(set(ge))

        merged = sort_reindex(merged, 'Score', 'P')
        st_in = sort_reindex(st_in, 'Score', 'P')
        g_p_val_P = calc_hypergeom_scores(merged, st_in, ge, frac)
        merged = sort_reindex(merged, 'Score', 'N')
        st_in = sort_reindex(st_in, 'Score', 'N')
        g_p_val_N = calc_hypergeom_scores(merged, st_in, ge, frac)

        with open(outputfile,'w') as o:
            w = csv.writer(o, delimiter='\t', lineterminator='\n')
            w.writerow(('Gene Symbol', 'Average LFC', 'Average -log(p-values)', 'Number of perturbations', 'Perturbations', 'Individual LFCs', 'Ascending ranks','Individual ascending -log(p-values)', 'Descending ranks', 'Individual descending -log(p-values)'))
            print 'Analyzing '+c
            for g in ge:
                in_lfcs, avg_lfc_p, p_p_vals, avg_p_val, guide_list, ranks_p = g_p_val_P[g].split('_')
                in_lfcs, avg_lfc_n, n_p_vals, avg_n_val, guide_list, ranks_n = g_p_val_N[g].split('_')

                if float(avg_p_val) > float(avg_n_val):
                    avg_p_val = float(avg_p_val)
                    avg_lfc = float(avg_lfc_p)
                else:
                    avg_p_val = float(avg_n_val)
                    avg_lfc = float(avg_lfc_n)
                if include == 'N':
                    if len(guide_list.split(';')) != 1:
                        w.writerow((g, avg_lfc, avg_p_val, len(in_lfcs.split(';')), guide_list, in_lfcs, ranks_n, n_p_vals, ranks_p, p_p_vals))
                else:
                    w.writerow((g, avg_lfc, avg_p_val, len(in_lfcs.split(';')), guide_list, in_lfcs, ranks_n, n_p_vals, ranks_p, p_p_vals))  
    with open(o_folder+'/README.txt','w') as o:
        w = csv.writer(o,delimiter='\t')
        w.writerow((['Code Version: 1.5.4']))
        w.writerow((['Input file:'+inputfile]))
        w.writerow((['Chip file:'+args.chip_file]))
        w.writerow((['Output folder:'+o_folder]))
        w.writerow((['Fraction of perturbations for mean p-value, LFC:'+str(frac)]))
                

