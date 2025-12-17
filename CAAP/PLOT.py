import matplotlib.pyplot as plt
import numpy as np 

def create_custom_xticks(x, sequence, interval=5):
    # Use sequence[i-1]because x is 1-indexed
    labels = [f'{sequence[i-1]}\n{i}' if i % interval == 0 else sequence[i-1] for i in x]
    return labels


# Average attention
def plot_attention(attention, important, protein, new_dir, sequence=None):
    x = np.arange(1, attention.size + 1 ) # 1-indexed positions

    plt.figure(figsize=(8, 6))
    plt.bar(x, attention, color='gray', zorder=1)
    #plt.bar(x, important, color='#CA9823', zorder=2)
    plt.bar(x, important, color='#AE0639', zorder=2)


    # Display amino acid and every 5th residue number
    if sequence:
        plt.xticks(x, create_custom_xticks(x, sequence))

    plt.xlabel('Amino acids')
    plt.ylabel('Average attention')
    plt.title(f'{protein}')

    plt.savefig(f'{new_dir}{protein}_average.png', dpi=600)
    #plt.savefig(f'{new_dir}Figure2_{protein}_average.png', dpi=600)
    plt.show()

def plot_difference(difference, protein, new_dir, sequence=None):
    x = np.arange(1, len(difference) + 1 ) # 1-indexed positions

    plt.figure(figsize=(8, 6))
    plt.bar(x, difference, color='gray')

    # Display amino acid and every 5th residue number
    if sequence:
        plt.xticks(x, create_custom_xticks(x, sequence))

    plt.xlabel('Amino acids')
    plt.ylabel('Positive difference')
    plt.title(f'{protein}')

    plt.savefig(f'{new_dir}{protein}_positive_difference.png', dpi=600)
    plt.show()

# # Plots average attention and the difference 
# def plot_comparison(attention1, attention2, important1, important2, protein1, protein2, 
#                     new_dir, sequence1=None, sequence2=None):
#     x1 = np.arange(1, attention1.size + 1 ) # 1-indexed positions
#     x2 = np.arange(1, attention2.size + 1 )

    
#     # Gets the max y values in order to scale y axis 
#     max_y_important = max(np.max(important1), np.max(important2))


#     plt.figure(figsize=(8, 6))

#     # Plot average attention for protein 1
#     plt.bar(x1, attention1, color='gray', width=0.4)
#     plt.xlabel('Residues')
#     plt.ylabel('Min-Max Normalized Average Attention')
#     plt.ylim(0, 1) 
#     plt.title(f'{protein1}')
#     if sequence1:
#         plt.xticks(x1, create_custom_xticks(x1, sequence1))
#     plt.savefig(f'{new_dir}{protein1}_aligned_average_attention.png', dpi=600)
#     plt.close()

#     # Plot important for protein 1
#     plt.bar(x1, important_diff1, width=0.4, color='b')
#     plt.bar(x1, important1, width=0.4, color='black', alpha=0.3)
#     plt.xlabel('Residues')
#     plt.ylabel('Min-Max Normalized Average Attention Positive Difference')
#     plt.ylim(0, max_y_important)  
#     plt.title(f'Important to {protein1}')
#     if sequence1:
#         plt.xticks(x1, create_custom_xticks(x1, sequence1))
#     plt.savefig(f'{new_dir}{protein1}_{protein2}_difference.png', dpi=600)
#     plt.close()

#     # Plot original normalized dataset 2
#     plt.bar(x2, attention2, width=0.4, color='gray')
#     plt.xlabel('Residues')
#     plt.ylabel('Min-Max Normalized Average Attention')
#     plt.ylim(0, 1) 
#     plt.title(f'{protein2}')
#     if sequence2:
#         plt.xticks(x2, create_custom_xticks(x2, sequence2))
#     plt.savefig(f'{new_dir}{protein2}_aligned_average_attention.png', dpi=600)
#     plt.close()
        
#     # Plot positive differences for data2 - data1
#     plt.bar(x2, important_diff2, width=0.4, color='r')
#     plt.bar(x2, important2, width=0.4, color='black', alpha=0.3)
#     plt.xlabel('Residues')
#     plt.ylabel('Min-Max Normalized Average Attention Positive Difference')
#     plt.ylim(0, max_y_important)  
#     plt.title(f'Important to {protein2}')
#     if sequence2:
#         plt.xticks(x2, create_custom_xticks(x2, sequence2))
#     plt.savefig(f'{new_dir}{protein1}_{protein2}_difference.png', dpi=600)
#     plt.close()



#     plt.show()
