import os
import jax
import numpy as np
import linecache

def get_n(folder_path):
    shape_counts = {}

    for fname in os.listdir(folder_path):
        if fname.endswith('.npy'):
            arr = np.load(os.path.join(folder_path, fname))
            shape = arr.shape
            shape_counts[shape] = shape_counts.get(shape, 0) + 1

    print("Unique shapes found:")
    for shape, count in shape_counts.items():
        print(f"{shape}: {count} files")
    shapes = list(shape_counts.keys())
    return shapes[0][0]

def get_attention(folder_path, n):
    # Extract file roots and sort numerically
    file_list = [fname for fname in os.listdir(folder_path) if fname.endswith('.npy')]
    file_roots = sorted(file_list, key=lambda x: int(x.split('_')[-1].split('.')[0]))
    #file_roots = sorted(file_list, key=lambda x: int(x.split('.')[0]))

    heads_n_4_n_n = []

    for fname in file_roots:
        arr1 = np.load(os.path.join(folder_path, fname))
        arr1 = arr1.view(dtype=np.float16) #key thing to display as numbers
        arr1 = jax.nn.softmax(arr1)
        if arr1.shape == (n, 4, n, n):
            heads_n_4_n_n.append(arr1)

    attention_spectrum = []
    for i in range(len(heads_n_4_n_n)):
        # aggregate of all the heads
        attn = np.sum(heads_n_4_n_n[i], axis=(0, 1, 2))
        attention_spectrum.append(attn)
    attention_spectrum = np.array(attention_spectrum)

    return attention_spectrum

def average(attention_spectrum):
    average_attention = np.mean(attention_spectrum, axis=0)
    return average_attention


def min_max(data):
    min_max = np.zeros([len(data)])
    min = data.min()
    max = data.max()
    
    for i, value in enumerate(data):
        min_max[i] = (value-min)/(max-min)
    return min_max


# Gets sequence 
def read_sequence_file(sequence_file):
    with open(sequence_file, 'r') as f:
        content = f.readlines()
        start_end = []
        found = 0
        sequence = ''
        for linenum, line in enumerate(content):
            # Finding first instance meaning next line is target sequence
            if line.startswith('>'):
                start_end.append(linenum)
                found +=1
            # Searching for next > meaning the rest of the MSA will be ignored if MSA is used
            if found == 2:
                break
        # if only 1 > then not an MSA
        if found == 1:
            # When using linecache it's not 0 indexed
            # +2 in order to start with the sequence line
            start_line = start_end[0]+2
            # Keep going until end of file in case the sequence got broken up into multiple lines
            while linecache.getline(sequence_file, start_line):
                sequence += linecache.getline(sequence_file, start_line).strip()
                start_line += 1
        # if 2 > found then an MSA was used 
        elif found == 2:
            start_line = start_end[0]
            end_line = start_end[1]
            # appending all lines to sequence after stripping
            # in case sequence got broken up into multiple lines
            for line in content[start_line+1:end_line]:
                sequence += line.strip()
    return sequence

# Reads the alignment to get aligned sequences with gaps
def read_alignment(protein1, protein2, alignment):
    with open(alignment, 'r') as f:
        content = f.readlines()
        for linenum, line in enumerate(content):
            if protein1.casefold() in line.casefold():
                aligned_sequence1 = content[linenum+1].strip()
            elif protein2.casefold() in line.casefold():
                aligned_sequence2 = content[linenum+1].strip()
    return aligned_sequence1, aligned_sequence2

# Aligns the attention values so that if there are gaps it's 0
def align_attention(attention1, attention2, sequence1, sequence2):
    # Taking the ratio between the protein lengths if alignment is necessary because attention on longer proteins will have 
    # lower absolute value
    og_seq1_len = len(attention1)
    og_seq2_len = len(attention2)

    seq_ratio = og_seq1_len/og_seq2_len

    aligned_attention1 = np.zeros([len(sequence1)])
    aligned_attention2 = np.zeros([len(sequence2)])

    # To keep track of the indices of the gaps
    gaps1 = []
    gaps2 = []
    # Index of non gaps that correspond with attention values
    i1_vals = 0
    i2_vals = 0
    for (i1, resi1), (i2, resi2) in zip(enumerate(sequence1), enumerate(sequence2)):
        if resi1 == '-':
            aligned_attention1[i1] = 0
            gaps1.append(i1)
            #print(resi1, aligned_attention1_avg_norm[i1])
        elif resi1.isalpha():
            aligned_attention1[i1] = attention1[i1_vals]
            #print(resi1, aligned_attention1_avg_norm[i1])
            i1_vals += 1

        if resi2 == '-':
            aligned_attention2[i2] = 0
            gaps2.append(i2)
            #print(resi2, aligned_attention2_avg_norm[i2])
        elif resi2.isalpha():
            aligned_attention2[i2] = attention2[i2_vals]
            #print(resi2, aligned_attention2_avg_norm[i2])
            i2_vals += 1
    
    print('gaps1: ', gaps1)
    print('gaps2: ', gaps2)

    # Because I'm always doing seq 1/seq 2
    # If seq 1 is longer then the ratio will be > 1
    # Need to multiple seq 1 by ratio as its longer and needs more attention
    # If seq 1 is shorter the ratio is < 1
    # Multiple seq 1 by the ratio as its shorter and needs less attention
    aligned_attention1 = aligned_attention1 * seq_ratio


    return aligned_attention1, aligned_attention2, gaps1, gaps2


