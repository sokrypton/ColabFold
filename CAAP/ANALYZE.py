import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

def find_highest_attention(attention, sequence, new_dir, protein):
    # Get sorted indices (descending order)
    descending_indices = (-attention).argsort()
    # Sorted attention, sequence, and original indices
    sorted_attention = attention[descending_indices]
    sorted_residue_numbers = descending_indices + 1  # 1-based indexing
    sorted_sequence = np.array(list(sequence))[descending_indices]

    # Build DataFrame
    df = pd.DataFrame({
        'Rank': np.arange(1, len(sequence)+1),  # 1-based rank
        'Residue number': sorted_residue_numbers,
        'Amino acid': sorted_sequence,
        'Attention score': sorted_attention  
    })

    df.to_csv(f'{new_dir}{protein}_residue_ranking.csv', index=False)


# Finds which positions have positive zscores meaning they're important
def find_important(attention, zscores):
    # This code finds the top 5
    important = np.zeros_like(attention)

    # Get indices of top_n attention values
    top_indices = np.argsort(attention)[-5:]
    important[top_indices] = attention[top_indices]

    return important
    # This code returns the ones that have z-scores above 1.5 or something
    # mean = np.mean(attention)
    # important = np.zeros([len(attention)])
    # for a, (i2, z) in zip(attention, enumerate(zscores)):
    #     if z > 1.5:
    #         important[i2] = a
    # return important


# Gets the blosum scores for each residue and then multiplies by the average attention
# to return 
def blosum_scores(sequence1, sequence2):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62') 
    scores1 = np.zeros([len(sequence1)])
    scores2 = np.zeros([len(sequence2)])

    for i, (s1, s2) in enumerate(zip(sequence1, sequence2)):
        if s1 == '-' and s2 == '-':
            s1_temp = '*'
            s2_temp = '*'
            scores1[i] = aligner.substitution_matrix[s1_temp, s2_temp] 
            scores2[i] = aligner.substitution_matrix[s1_temp, s2_temp] 
        elif s1 == '-':
            s1_temp = '*'
            scores1[i] = aligner.substitution_matrix[s1_temp, s2] 
        elif s2 == '-':
            s2_temp = '*'
            scores2[i] = aligner.substitution_matrix[s1, s2_temp] 
        else:
            scores1[i] = aligner.substitution_matrix[s1, s2] 
            scores2[i] = aligner.substitution_matrix[s1, s2] 

    return scores1, scores2


def calculate_differences(scores1, scores2, attention1, attention2, gaps1, gaps2, bool_alignment):
    # If there's no alignment then just take the differences
    if bool_alignment:
        # Poositive difference for protein 1, meaning important for 1
        difference1 = attention1 - attention2
        important_diff1 = np.where(difference1 > 0, difference1, 0)
        # After taking the difference, setting all that are gaps to 0 so it doesn't skew distribution
        # Set the positions where the gaps are in 2 equal to 0 for 1
        for g in gaps2:
            important_diff1[g] = 0
        # Important for 2
        difference2 = attention2 - attention1
        important_diff2 = np.where(difference2 > 0, difference2, 0)
        # Where the gaps are in 1, set them equal to 0 in 2
        for g in gaps1:
            important_diff2[g] = 0

    else:
        # Poositive difference for protein 1, meaning important for 1
        difference1 = attention1 - attention2
        important_diff1 = np.where(difference1 > 0, difference1, 0)

        # Important for 2
        difference2 = attention2 - attention1
        important_diff2 = np.where(difference2 > 0, difference2, 0)
    
    # Multiplying the differences by the blosum scores
    important_diff_blosum1 = [diff * s for diff, s in zip(important_diff1, scores1)]
    important_diff_blosum2 = [diff * s for diff, s in zip(important_diff2, scores2)]
    return important_diff_blosum1, important_diff_blosum2
