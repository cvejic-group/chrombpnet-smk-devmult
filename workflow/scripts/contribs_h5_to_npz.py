import numpy as np
import deepdish

shap = deepdish.io.load(snakemake.input['h5'], '/shap/seq')
onehot = deepdish.io.load(snakemake.input['h5'], '/raw/seq')

np.savez(snakemake.output['shap'], shap)
np.savez(snakemake.output['onehot'], onehot)
