message("Loading libraries...")
suppressPackageStartupMessages({
  library(magrittr)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
  library(stringr)
  library(plyranges)
})

if (interactive()) { ## fake data for testing
  message("Mocking Snakemake object...")
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric",
                                             wildcards="list"))
  snakemake <- Snakemake(
    input=list(bedpe="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/Large-PreB/mean/modisco/Large-PreB_mean.counts_scores.top_tf.bedpe.gz",
               bed="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/IM-B/mean/modisco/IM-B_mean.counts_scores.sorted.bed.gz"),
    output=list(rds="/work/aaa/projects/chrombpnet-devmult/pipeline/results/overlaps/chrombpnet_nobias/pretrained_bias/celltypes/counts/intersect_Large-PreB_IM-B.Rds"),
    params=list(sample1="Large-PreB", sample2="IM-B"),
    threads=1L
  )
}

message("Inspecting Snakemake object...")
print(snakemake)

message("Defining helper functions...")
read_bedpe_seqlet <- function (file) {
  read_tsv(file,
           col_names=c("seqnames", "start", "end",
                       "name", "strand", "pattern_id"),
           col_types="cii---c-f--c") %>%
    as_granges()
}

message("Loading data...")
gr_seqlets1 <- read_bedpe_seqlet(snakemake@input$bedpe)
gr_seqlets2 <- read_bed(snakemake@input$bed)

message("Counting overlaps...")
cts_seqlet_celltype <- gr_seqlets1 %>%
  mutate(pattern_group=str_c(snakemake@params$sample1, ".", str_remove(pattern_id, ".[0-9]+$"), ".", name)) %>%
  mutate(pattern_group=str_remove(pattern_group, "patterns.")) %>%
  split(f=.$pattern_group) %>%
  lapply(function(gr) {
    filter_by_overlaps(gr, gr_seqlets2, minoverlap=20) %>% length()
  }) %>% {
    matrix(unlist(.), ncol=1, dimnames=list(names(.), snakemake@params$sample2))
  }

message("Exporting result...")
saveRDS(cts_seqlet_celltype, snakemake@output$rds)

message("Done.")
