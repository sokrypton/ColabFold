import numpy as np
from scipy.stats import zscore
import pandas as pd
import matplotlib.pyplot as plt
import os
import jax
import argparse
import linecache
import sys
import PROCESS
import ANALYZE
import PLOT



def main():
    parser = argparse.ArgumentParser(description="Visualize attention heads.")
    parser.add_argument('--attention1', help='Path to the attention heads folder for the first protein.')
    parser.add_argument('--protein1', help='Name of the first protein.')
    parser.add_argument('--sequence1', help='Path to the sequence file (.a3m or .fasta) for the first protein.')
    parser.add_argument('--attention2', help='Path to the attention heads folder for the second protein.', default=None)
    parser.add_argument('--protein2', help='Name of the second protein.', default=None)
    parser.add_argument('--sequence2', help='Path to the sequence file (.a3m or .fasta) for the second protein.', default=None)
    parser.add_argument('--alignment', help='Path to the alignment file.', default=None)
    parser.add_argument('--save_path', help='Path to save output')

    args = parser.parse_args()

    sequence1 = PROCESS.read_sequence_file(args.sequence1)
    n1 = PROCESS.get_n(args.attention1)
    attention_spectrum1 = PROCESS.get_attention(args.attention1, n1)
    average_attention1 = PROCESS.average(attention_spectrum1)
    min_max_attention1 = PROCESS.min_max(average_attention1)
    zscores1 = zscore(min_max_attention1)
    important1 = ANALYZE.find_important(min_max_attention1, zscores1)

    # Only 1 protein then just show average attention
    if args.protein2 is None and args.attention2 is None:
        if args.save_path.endswith('/'):
            new_dir = f'{args.save_path}{args.protein1}/'
        else:
            new_dir = f'{args.save_path}/{args.protein1}/'
        os.makedirs(new_dir, exist_ok=True)
        ANALYZE.find_highest_attention(average_attention1, sequence1, new_dir, args.protein1)
        PLOT.plot_attention(min_max_attention1, important1, args.protein1, new_dir, sequence=sequence1)
    # An alignment isn't necessary for 2 proteins because they could already be the same length
    elif args.protein2 and args.attention2:
        sequence2 = PROCESS.read_sequence_file(args.sequence2)
        if (len(sequence1) != len(sequence2)) and args.alignment is None:
            print('Sequences are differeng lengths. Must input MSA')
            sys.exit()

        n2 = PROCESS.get_n(args.attention2) 
        attention_spectrum2 = PROCESS.get_attention(args.attention2, n2)
        average_attention2 = PROCESS.average(attention_spectrum2)
        if args.save_path.endswith('/'):
            new_dir = f'{args.save_path}{args.protein1}_{args.protein2}/'
        else:
            new_dir = f'{args.save_path}/{args.protein1}_{args.protein2}/'
        os.makedirs(new_dir, exist_ok=True)


        if len(sequence1) == len(sequence2):
            min_max_attention2 = PROCESS.min_max(average_attention2)
            zscores2 = zscore(min_max_attention2)
            important2 = ANALYZE.find_important(min_max_attention2, zscores2)
            blosum1, blosum2 = ANALYZE.blosum_scores(sequence1, sequence2)
            diff1, diff2 = ANALYZE.calculate_differences(blosum1, blosum2, min_max_attention1, min_max_attention2, 
                                                         [], [], False)
            ANALYZE.find_highest_attention(average_attention1, sequence1, new_dir, args.protein1)
            ANALYZE.find_highest_attention(average_attention2, sequence2, new_dir, args.protein2)
            PLOT.plot_attention(min_max_attention1, important1, args.protein1, new_dir, sequence=sequence1)
            PLOT.plot_attention(min_max_attention2, important2, args.protein2, new_dir, sequence=sequence2)
            PLOT.plot_difference(diff1, args.protein1, new_dir, sequence=sequence1)
            PLOT.plot_difference(diff2, args.protein2, new_dir, sequence=sequence2)
            # PLOT.plot_comparison(min_max_attention1, min_max_attention2, important1, important2, 
            #                 args.protein1, args.protein2, [], [], new_dir, False, sequence1=sequence1, sequence2=sequence2)
        elif args.alignment and args.protein2:
            aligned_sequence1, aligned_sequence2 = PROCESS.read_alignment(args.protein1, args.protein2, args.alignment)
            aligned_average_attention1, aligned_average_attention2, gaps1, gaps2 = PROCESS.align_attention(average_attention1, 
                                                                                                           average_attention2, aligned_sequence1, aligned_sequence2)
            aligned_min_max_attention1 = PROCESS.min_max(aligned_average_attention1)
            aligned_min_max_attention2 = PROCESS.min_max(aligned_average_attention2)
            aligned_zscores1 = zscore(aligned_min_max_attention1)
            aligned_zscores2 = zscore(aligned_min_max_attention1)
            aligned_important1 = ANALYZE.find_important(aligned_min_max_attention1, aligned_zscores1)
            aligned_important2 = ANALYZE.find_important(aligned_min_max_attention2, aligned_zscores2)
            aligned_blosum1, aligned_blosum2 = ANALYZE.blosum_scores(aligned_sequence1, aligned_sequence2)
            aligned_diff1, aligned_diff2 = ANALYZE.calculate_differences(aligned_blosum1, aligned_blosum2, aligned_min_max_attention1, aligned_min_max_attention2, 
                                                         gaps1, gaps2, True)
            ANALYZE.find_highest_attention(aligned_average_attention1, aligned_sequence1, new_dir, args.protein1)
            ANALYZE.find_highest_attention(aligned_average_attention2, aligned_sequence2, new_dir, args.protein2)
            PLOT.plot_attention(aligned_min_max_attention1, aligned_important1, args.protein1, new_dir, sequence=aligned_sequence1)
            PLOT.plot_attention(aligned_min_max_attention2, aligned_important2, args.protein2, new_dir, sequence=aligned_sequence2)
            PLOT.plot_difference(aligned_diff1, args.protein1, new_dir, sequence=aligned_sequence1)
            PLOT.plot_difference(aligned_diff2, args.protein2, new_dir, sequence=aligned_sequence2)
            # PLOT.plot_comparison(aligned_min_max_attention1, aligned_min_max_attention2, aligned_important1, aligned_important2, 
            #                 args.protein1, args.protein2, gaps1, gaps2, new_dir, True, sequence1=aligned_sequence1, sequence2=aligned_sequence2)

        # If an alignment is given then there needs to be 2 proteins
        elif args.alignment and (args.attention2 is None or args.protein2 is None):
            parser.error('--alignment requires 2 proteins')




if __name__ == '__main__':
    main()

