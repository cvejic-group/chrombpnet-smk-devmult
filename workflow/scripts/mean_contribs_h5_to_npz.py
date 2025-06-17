import numpy as np
import deepdish

shap_dict = {}
projected_shap_dict = {}

for i, h5 in enumerate(snakemake.input['h5']):
    shap_dict[i] = deepdish.io.load(h5, '/shap/seq')
    projected_shap_dict[i] = deepdish.io.load(h5, '/projected_shap/seq')

mean_shap = np.mean(np.array([shap_dict[fold] for fold in shap_dict]), axis=0)
mean_projected_shap = np.mean(np.array([projected_shap_dict[fold] for fold in projected_shap_dict]), axis=0)

onehot = deepdish.io.load(snakemake.input['h5'][0], '/raw/seq')

d = {
    'raw': {'seq': onehot},
    'shap': {'seq': mean_shap},
    'projected_shap': {'seq': mean_projected_shap}
}

deepdish.io.save(snakemake.output['h5'], d, compression='blosc')

np.savez(snakemake.output['shap'], mean_shap)
np.savez(snakemake.output['onehot'], onehot)
