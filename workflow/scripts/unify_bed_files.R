message("Loading libraries...")
library(tidyverse)
library(magrittr)
library(plyranges)
library(BiocParallel)

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(bed_agg="/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC/HSC.no_blacklist.bed",
               beds_pcw=c(
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW5/HSC_PCW5.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW6/HSC_PCW6.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW7/HSC_PCW7.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW8/HSC_PCW8.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW9/HSC_PCW9.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW10/HSC_PCW10.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW11/HSC_PCW11.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW12/HSC_PCW12.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW13/HSC_PCW13.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW14/HSC_PCW14.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW15/HSC_PCW15.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW16/HSC_PCW16.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW17/HSC_PCW17.no_blacklist.bed",
                 "/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/HSC_PCW18/HSC_PCW18.no_blacklist.bed"
               )),
    output=list(bed="/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/unified/HSC.unified.bed"),
    params=list(sample="HSC"),
    threads=8L)
}

message("Inspecting `snakemake` object...")
print(snakemake)

message("Initializing threading...")
register(MulticoreParam(snakemake@threads))

message("Loading aggregated BED peaks...")
gr_agg <- read_narrowpeaks(snakemake@input$bed_agg)

message("Loading individual BED peaks...")
grl_pcw <- lapply(snakemake@input$beds_pcw, read_narrowpeaks) %>%
  as("GRangesList")

message("Filtering to unique peaks...")
gr_unique <- lapply(grl_pcw, filter_by_non_overlaps, y=gr_agg, minoverlap=250L) %>%
  as("GRangesList") %>% unlist()

gr_summits <- resize(gr_unique, width=1, fix='center')
gr_clusters <- reduce(gr_summits, min.gapwidth=100)



hits <- gr_clusters %>%
  findOverlaps(gr_summits, ignore.strand=TRUE) %>%
  { split(subjectHits(.), queryHits(.)) }

idx_summits <- bplapply(hits, function(idxs) {
  if (length(idxs) == 1) return(idxs)

  gr_sub <- gr_summits[idxs]
  pos_summit <- start(gr_sub)
  pos_median <- median(pos_summit)

  dists <- abs(pos_summit - pos_median)
  closest <- which(dists == min(dists))

  if (length(closest) == 1) {
    return(idxs[closest])
  } else {
    scores <- mcols(gr_sub)$score[closest]
    return(idxs[closest[which.max(scores)]])
  }
})

gr_cluster_peaks <- gr_summits[unlist(idx_summits, use.names=FALSE)] %>%
  resize(width=500L, fix='center')

gr_unified <- c(gr_agg, gr_cluster_peaks)

write_narrowpeaks(gr_unified, snakemake@output$bed)

# gr_summits[unlist(hits, use.names=FALSE)] %>%
#   resize(width=500L, fix='center') %>%
#   write_narrowpeaks(file="scratch/bed/HSC.unified.raw_clusters.bed")
