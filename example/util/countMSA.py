import json
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('--input', '-i', required=True)

args = argparser.parse_args()

with open(args.input) as f:
    data = json.load(f)
    seqs = data['sequences']
    for seq in seqs:
        prot = seq['protein']
        id = prot['id']
        length = len(prot['sequence'])
        unpaired = prot['unpairedMsa']
        unpaired_lines = len([x for x in unpaired.split('\n') if x.startswith('>')])
        unpaired_lengths = [len(x) for x in unpaired.split('\n') if not x.startswith('>')]

        paired = prot['pairedMsa']
        paired_lines = len([x for x in paired.split('\n') if x.startswith('>')])
        paired_lengths = [len(x) for x in paired.split('\n') if not x.startswith('>')]
        print(f'{id}: {length} residues,{unpaired_lines} unpaired, {paired_lines} paired')
        print(f'Unpaired lengths: {unpaired_lengths}')
        print(f'Paired lengths: {paired_lengths}')
        