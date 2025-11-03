message("Loading libraries...")
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(SummarizedExperiment)
})

if (interactive()) { ## fake data for testing
  message("Mocking Snakemake object...")
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric",
                                             wildcards="list"))
  snakemake <- Snakemake(
    input=list(rdss=c(
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_HSC_HSC.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_HSC_IM-B.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_HSC_Large-PreB.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_IM-B_HSC.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_IM-B_IM-B.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_IM-B_Large-PreB.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_Large-PreB_HSC.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_Large-PreB_IM-B.Rds",
      "/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_Large-PreB_Large-PreB.Rds"
    )),
    output=list(rds="results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/se_seqlets_celltypes.Rds"),
    params=list(),
    threads=1L
  )
}

message("Inspecting Snakemake object...")
print(snakemake)

message("Defining helper functions...")
order_pattern_str <- function (x) {
  m <- str_match(x, "^[^.]+\\.(pos|neg)_pattern_(\\d+)\\.[^.]+$")
  type_pattern  <- factor(m[,2], levels=c("pos", "neg"))
  idx_pattern  <- as.numeric(m[,3])

  x[order(type_pattern, idx_pattern)]
}

message("Preparing matrix blocks and indices...")
blocks_row <- sort(snakemake@input$rdss) %>%
  split(str_replace_all(., "(^.*intersect_|_[^_]+.Rds)", ""))

idx_cols <- names(blocks_row)
idx_rows <- lapply(blocks_row, function (files)  {
  readRDS(files[[1]]) %>%
    rownames() %>%
    order_pattern_str()
})

message("Loading data...")
mat_pattern_sample <- lapply(names(blocks_row), function (idx_celltype) {
  files <- blocks_row[[idx_celltype]]
  names_cols <- str_replace_all(files, "(^.*intersect_[^_]+_|.Rds)", "")
  lapply(files, function(f) {
    readRDS(f)[idx_rows[[idx_celltype]],]
  }) %>%
    do.call(what=cbind) %>%
    `colnames<-`(names_cols) %>%
    `[`(,idx_cols)
}) %>% do.call(what=rbind)

message("Constructing SummarizedExperiment object...")
df_rowdata <- tibble(pattern_id=rownames(mat_pattern_sample),
                     n_seqlets=rowMaxs(mat_pattern_sample)) %>%
  mutate(sample_id=str_extract(pattern_id, "^[^.]+"),
         top_tf=str_extract(pattern_id, "[^.]+$")) %>%
  DataFrame(row.names=.$pattern_id)
df_coldata <- tibble(sample_id=idx_cols,
                     n_seqlet_patterns=sapply(idx_cols, function (idx) {
                       sum(str_detect(rownames(mat_pattern_sample), str_c("^", idx, "\\.")))
                     }),
                     n_seqlets=sapply(idx_cols, function (idx) {
                       idx_sample <- str_detect(rownames(mat_pattern_sample), str_c("^", idx, "\\."))
                       sum(mat_pattern_sample[idx_sample, idx])
                     })) %>%
  DataFrame(row.names=.$sample_id)
se_pattern_sample <- SummarizedExperiment(assays=list(counts=mat_pattern_sample),
                                           rowData=df_rowdata,
                                           colData=df_coldata)

message("Exporting result...")
saveRDS(se_pattern_sample, snakemake@output$rds)

message("Done.")
