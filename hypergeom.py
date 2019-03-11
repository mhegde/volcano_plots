'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
Performs hypergeometric analysis for genetic perturbation screens
Input: 1. Input file with column headers
       2. Chip File
       3. Option to include/exclude singletons in output file;Default: Singletons included
       4. Option to specify n% for mean p-val, LFC calculation
       5. Option of number of control guides to be included in a random set
       6. Option of minimum number of perturbations for a gene to be included in the volcano plot
       7. Option of maximum number of perturbations for a gene to be included in the volcano plot
       8. Option of number of genes to be labeled on the plot
       9. Option to include or exclude No_Site controls in plot
       10. Option to include or exclude One_Non-Gene_Site controls in plot
'''
import pandas as pd
import numpy as np
from scipy import stats
from math import log10
import csv, argparse, sys, re, traceback, os
from datetime import datetime
from decimal import Decimal
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
mpl.use('Agg')
mpl.rc('pdf',fonttype=42)
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['axes.unicode_minus']=False
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import warnings

warnings.filterwarnings("ignore")
log = open('error.log','w')

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, help='File containing sgRNA sequence and score(s) columns with headers')
    parser.add_argument('--chip-file', type=str, help='Chip file')
    parser.add_argument('--sing-pert', type=str, default='Y', help='Y to include genes with single perturbations in output file, N to exclude; Default: Y')
    parser.add_argument('--fraction', type=int, default=100, help='Fraction of perturbations to be included in mean LFC, pval calculation, 1 - 100; Default: 100')
    parser.add_argument('--ctrls-num', type=int, default=4,help='Number of control guides to be included in a random set')
    parser.add_argument('--min-pert', type=int, default=3, help='Minimum number of perturbations for a gene to be included in the volcano plot')
    parser.add_argument('--max-pert', type=int, default=8, help='Maximum number of perturbations for a gene to be included in the volcano plot')
    parser.add_argument('--label-num', type=int, default=3, help='Number of genes to be labeled on the plot')
    parser.add_argument('--disp-nosite', type=str, default='Y', choices=['Y','N'],help='Y to include No_Site controls in plot')
    parser.add_argument('--disp-onesite', type=str, default='Y', choices=['Y','N'],help='Y to include One_Non-Gene_Site controls in plot')
    return parser

'''
Sorts a data frame on the column specified and re-indexes
Argument: Dataframe, column to sort on, direction of sorting
Return Value: Sorted and re-index dataframe
'''
def sort_reindex(df, col, direction):
    if direction == 'P' or direction == 'p':
        df = df.sort_values(by=col, ascending=False)
        df.index = range(0, len(df))
        #df['Rank'] = map(lambda x:x+1, range(len(df)))
    elif direction == 'N' or direction == 'n':
        df = df.sort_values(by=col, ascending=True)
        df.index = range(0, len(df))
    else:
        print('Please enter a relevant direction; P for positive and N for negative')
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


def get_random_sets(df, num, label):
    random_assign = []
    for x in range(1,(len(df)/num)+1):
        random_assign.extend([x]*num)
    rem_df = len(df) - len(random_assign)
    random_assign.extend(range(random_assign[-1] + 1, random_assign[-1] + 1 + rem_df))
    df.index = range(0,len(df))
    df['random_assign'] = random_assign
    new_df = pd.DataFrame()
    for i in list(set(random_assign)):
        ra_df = df[df['random_assign'] == i]
        ra_df['Gene Symbol'] = label+'_'+str(i)
        #if len(ra_df) == num:
        new_df = new_df.append(ra_df)
    new_df = new_df[['Guide Sequence','Gene Symbol']]
    return new_df


def generate_chip(chip,num):
    print('Generating temp chip file...')
    nosite_ctrls = chip[chip['Gene Symbol'].str.startswith('NO_SITE')]
    nongene_ctrls = chip[chip['Gene Symbol'].str.startswith('ONE_NON-GENE_SITE')]
    new_chip = chip[(~chip['Gene Symbol'].str.startswith('NO_SITE')) & (~chip['Gene Symbol'].str.startswith('ONE_NON-GENE_SITE'))]
    new_chip = new_chip[['Guide Sequence','Gene Symbol']]
    if len(nosite_ctrls) != 0:
        nosite_random = get_random_sets(nosite_ctrls, num, 'NO_SITE')
        if len(nosite_random) > 0:
            new_chip = new_chip.append(nosite_random)
    if len(nongene_ctrls) != 0:
        nongene_random = get_random_sets(nongene_ctrls, num, 'ONE_NON-GENE_SITE')
        if len(nongene_random) > 0:
            new_chip = new_chip.append(nongene_random)
    return new_chip

def plot_volcano(outputfile, min_pert, max_pert,label_num, c, disp_nosite, disp_onesite):
    print('Generating volcano plot...')
    output_df = pd.read_table(outputfile)
    output_df = output_df[(output_df['Number of perturbations']>=min_pert)&(output_df['Number of perturbations']<=max_pert)]
    nosite_ctrls = output_df[output_df['Gene Symbol'].str.contains('NO_SITE')]
    nongene_ctrls = output_df[output_df['Gene Symbol'].str.contains('ONE_NON-GENE_SITE')]
    fig = plt.figure()
    ax = plt.subplot2grid((1,3),(0,0), colspan=2)
    ax1 = plt.subplot2grid((1,3),(0,2))
    ax.scatter(output_df['Average LFC'],output_df['Average -log(p-values)'],color='black',zorder=1,s=30)
    legend_elements = []
    if disp_onesite == 'Y':
        ax.scatter(nongene_ctrls['Average LFC'], nongene_ctrls['Average -log(p-values)'], color='lightgrey', zorder=3,s=30)
        legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label='Controls'))
        legend_elements.append(Line2D([],[], marker='o', color='lightgrey', label='ONE_NON-GENE_SITE',linestyle='None'))
    if disp_nosite == 'Y':
        ax.scatter(nosite_ctrls['Average LFC'],nosite_ctrls['Average -log(p-values)'],color='dimgrey',zorder=2,s=30)
        if len(legend_elements) == 0:
            legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label='Controls'))
        legend_elements.append(Line2D([0], [0], marker='o', color='dimgrey', label='NO_SITE',linestyle='None'))
    ax.set_title(c,fontname="Helvetica",fontsize=12)
    ax.set_xlabel('Average fold change(log2)',fontname="Helvetica",fontsize=12)
    ax.set_ylabel('Average p-value(-log10)',fontname="Helvetica",fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for tick in ax.get_xticklabels():
        tick.set_fontname('Helvetica')
        tick.set_fontsize(12)
    for tick in ax.get_yticklabels():
        tick.set_fontname('Helvetica')
        tick.set_fontsize(12)
    if label_num != 0:
        output_df = output_df.sort(columns='Average LFC',ascending=False)
        color_list = plt.cm.Set1(np.linspace(0, 1, 12))
        color_index = 0
        top_hits = output_df.head(label_num)
        legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0))
        legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label='Top enrichment hits'))
        for i,r in top_hits.iterrows():
            ax.scatter(r['Average LFC'],r['Average -log(p-values)'],c=color_list[color_index],s=50)
            legend_elements.append(Line2D([],[], marker='o', color=color_list[color_index], label=r['Gene Symbol'], linestyle='None', markeredgecolor='black'))
            color_index+=1
        bottom_hits = output_df.tail(label_num)
        legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0))
        legend_elements.append(Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label='Top depletion hits'))
        for i,r in bottom_hits.iterrows():
            ax.scatter(r['Average LFC'],r['Average -log(p-values)'],c=color_list[color_index],s=50)
            legend_elements.append(Line2D([0], [0], marker='o', color=color_list[color_index], label=r['Gene Symbol'], linestyle='None', markeredgecolor='black'))
            color_index+=1
    ax1.legend(handles=legend_elements, loc='center',frameon=False, numpoints = 1, fontsize=10, handletextpad=0.0)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.tick_params(axis='both', which='both', labelbottom=False, labelleft=False, bottom=False, left=False, top=False, right=False)
    plt.tight_layout()
    fig.savefig(os.path.dirname(os.path.abspath(outputfile))+'/'+c+'.pdf')
    return 1


def input_check(input_df):
    for i, x in enumerate(input_df.columns[1:]):
        if input_df.dtypes.values[i + 1] != 'float64':
            sys.exit('Only the first column of the input file should contain non-numeric values; Please check your input file.')
    if input_df.isnull().values.any():
        sys.exit('Please remove empty spaces/NaNs from your input file.')
    return


def read_args(args):
    inputfile = args.input_file
    input_df = pd.read_table(inputfile)
    cols = list(input_df.columns)[1:]
    cols = [re.sub('[\s+\.]', '_', x) for x in cols]
    cols.insert(0, 'Guide Sequence')
    input_df.columns = cols
    input_check(input_df)
    ref = pd.read_table(args.chip_file)
    ref_colnames = list(ref.columns)
    ref_colnames[0:2] = ['Guide Sequence', 'Gene Symbol']
    ref.columns = ref_colnames
    include = args.sing_pert
    frac = args.fraction
    num = args.ctrls_num
    min_pert = args.min_pert
    max_pert = args.max_pert
    label_num = args.label_num
    disp_nosite = args.disp_nosite
    disp_onesite = args.disp_onesite
    return inputfile, input_df, ref, include, frac, num, min_pert, max_pert, label_num, disp_nosite, disp_onesite


if __name__ == '__main__':
    try:
        args = get_parser().parse_args()
        inputfile, input_df, ref, include, frac, num, min_pert, max_pert, label_num, disp_nosite, disp_onesite = read_args(args)
        cols_iter = list(input_df.columns)[1:]
        inputname = inputfile.split('/')[-1]
        ref = generate_chip(ref,num)
        o_folder = 'Volcano-plots_'+inputname[:-4]+'_'+str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
        if not os.path.exists(o_folder):
            os.makedirs(o_folder)
        ref.to_csv(o_folder+'/temp_chip_file.txt',sep='\t',index=False)
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
                print('Analyzing '+c)
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
            val = plot_volcano(outputfile, min_pert, max_pert, label_num, c, disp_nosite, disp_onesite)
    except Exception as e:
        traceback.print_exc(file=log)
        sys.exit('Error processing request.')
