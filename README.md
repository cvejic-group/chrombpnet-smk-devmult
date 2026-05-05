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

## Pipeline Setup
In practice, the pipeline used a mixture of Snakemake-managed Conda environments and user-managed Conda environments.
The user-managed Conda environments are documented in the `workflow/envs` folder, with the **".full"** file suffix,
and capture all software in those environments. These can only be replicated on Linux systems:

```bash
## chrombpnet_0_1_7
conda env create -n chrombpnet_0_1_7 -f workflow/envs/chrombpnet_0_1_7.full.yaml

## finemo_0_25
conda env create -n finemo_0_25 -f workflow/envs/finemo_0_25.full.yaml

## modiscolite_2_2_1
conda env create -n modiscolite_2_2_1 -f workflow/envs/modiscolite_2_2_1.full.yaml
```

Furthermore, the **modiscolite_2_2_1** environment included patching to the `modiscolite` module, which is provided in `workflow/envs/modiscolite_2.1.1.patch`. To apply this patch, it should be run from the `lib/python3.10/site-packages/` directory. E.g.,

```bash
## move to the env directory
pushd ${CONDA_PREFIX}/envs/modiscolite_2_2_1/lib/python3.10/site-packages
  ## check first
  git apply --stat --check /path/to/chrombpnet-smk-devmult/workflow/envs/modiscolite_2.1.1.patch
  # modiscolite/io.py         |    4 ++--
  # modiscolite/bed_writer.py |    4 ++--
  # 2 files changed, 4 insertions(+), 4 deletions(-)

  ## if looks like above, then apply
  git apply /path/to/chrombpnet-smk-devmult/workflow/envs/modiscolite_2.1.1.patch

## change back to original directory
popd
```

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
