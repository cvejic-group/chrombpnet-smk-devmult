[![DOI](https://zenodo.org/badge/1003459024.svg)](https://doi.org/10.5281/zenodo.20024504)

# chrombpnet-smk
This is a Snakemake pipeline for running ChromBPNet, TF-MoDiScO, and Fi-NeMo on pseudobulked snATAC-seq data.

## Overview
The pipeline uses the arguments in `config/config.yaml` to define what files are going to be used in TF footprinting.
The pipeline primarily runs using a pretrained ATAC bias model, however, some rules are implemented to train *de novo* bias. In practice, we found that the pretrained bias model did not lead to any substantial residual Tn5 footprints in the no bias model, and therefore was sufficient. The Snakefile `MAIN` section 
defines what outputs will be generated.

### Inputs
- `samples` - dictionary of deduplicated BAM files
- `peak_bed` - dictionary of BED files to be used in model training
- `peak_gene_link` (optional) - dictionary of BEDPE files defining peak-gene links

### Outputs
- standard ChromBPNet, TF-MoDiScO, and Fi-NeMo outputs
- additional BED, BEDPE, and BigWig files representing the putative footprints and contribution scores

## Running the Pipeline
The pipeline has been previously run with Snakemake v8.27 in SLURM-based mode.
It is expected to run with Conda+Mamba (Miniforge installation).
Note that rules that use GPU include a `resources:` entry:

```
slurm_extra="'--gres=gpu:h100:1'"
```

to request a GPU in SLURM. Such rules may need adjusting if not using SLURM.

### Install Snakemake from YAML
YAML representations of the Conda environment with Snakemake v8.27 are provided in minimal (".min") 
and full (".full") forms for reuse and reproduction, respectively. The minimal environment could be
created with

```bash
conda env create -n smk_8_27 -f workflow/envs/smk_8_27.min.yaml
```

### SLURM Execution
To run on a SLURM configuration, first configure a SLURM profile for Snakemake.
We have used

**~/.config/snakemake/slurm_basic/config.v8+.yaml**
```yaml
executor: slurm
use-conda: true
jobs: 10000
default-resources:
  mem_mb_per_cpu: 4096
```

which can be used with:

```bash
snakemake --profile slurm_basic -s workflow/Snakefile
```

### Generating Individual Outputs
As a Snakemake pipeline, one can generate individual files *ad hoc* by specifying them in the command.
For example, one could request a BED track of seqlets for predicting accessibility in the "MEMP" sample with

```bash
snakemake --profile slurm_basic \
  -s workflow/Snakefile \
  results/chrombpnet_nobias/pretrained_bias/MEMP/mean/modisco/MEMP_mean.counts_scores.top_tf.bed.bgz
```
