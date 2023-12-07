import numpy as np
import sys
tmp_params = np.load(sys.argv[1], allow_pickle=True)['arr_0'].flat[0]
params = {}
for k, v in tmp_params.items():
    for i, j in v.items():
        new_key = k.replace('deepfold', 'alphafold').replace('alphafold_batch/', '') + '//' + i
        params[new_key] = j
np.savez(sys.argv[2], **params)