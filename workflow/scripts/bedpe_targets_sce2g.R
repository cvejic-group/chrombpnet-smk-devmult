message("Loading libraries...")
library(tidyverse)
library(magrittr)
library(plyranges)

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(bedpe_seqlet="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/pretrained_bias/IM-B/mean/modisco/IM-B_mean.counts_scores.top_tf.bedpe.gz",
               bedpe_sce2g="/work/DevM_analysis/data/E2G/scE2G/IM-B_scE2G.enh2gene.bedpe"),
    output=list(bedpe="/work/aaa/projects/chrombpnet-devmult/pipeline/results/seqlet-target/scE2G/IM-B_seqlet_scE2G.counts_scores.bedpe.gz",
                csv="/work/aaa/projects/chrombpnet-devmult/pipeline/results/seqlet-target/scE2G/IM-B_seqlet_scE2G.counts_scores.csv.gz"),
    params=list(sample="IM-B"),
    threads=1L)
}

message("Printing snakemake object...")
print(snakemake)

message("Reading seqlets BEDPE...")
df_seqlet <- read_tsv(snakemake@input$bedpe_seqlet,
                      col_names=c("chrom1", "start1", "end1",
                                  "chrom2", "start2", "end2",
                                  "name", "score",
                                  "strand1", "strand2",
                                  "target_score", "modisco_id")) %>%
  mutate(idx_seqlet=row_number())

message("Reading P2G BEDPE...")
df_sce2g <- read_tsv(snakemake@input$bedpe_sce2g,
                     col_names=c("chrom1", "start1", "end1",
                                 "chrom2", "start2", "end2",
                                 "name", "score",
                                 "strand1", "strand2",
                                 "enhancer_type", "ensembl_id")) %>%
  mutate(idx_enhancer=row_number())

message("Converting seqlet peaks to GRanges...")
gr_peak <- df_seqlet %>%
  dplyr::select(chrom2, start2, end2, idx_seqlet) %>%
  dplyr::rename(seqnames=chrom2, start=start2, end=end2) %>%
  as_granges()

message("Converting enhancer peaks to GRanges...")
gr_enhancer <- df_sce2g %>%
  dplyr::select(chrom1, start1, end1, idx_enhancer) %>%
  dplyr::rename(seqnames=chrom1, start=start1, end=end1) %>%
  as_granges()

message("Intersecting seqlet and enhancer peaks...")
df_overlap <- join_overlap_inner(gr_peak, gr_enhancer) %>%
  as_tibble %>%
  dplyr::select(idx_seqlet, idx_enhancer)

message("Generating seqlet-target tibble...")
df_combined <- df_overlap %>%
  left_join(df_seqlet, by="idx_seqlet") %>%
  left_join(df_sce2g, by="idx_enhancer", suffix = c("_seqlet", "_target")) %>%
  mutate(target_name=str_extract(name_target, "[^:]+$"),
         pattern_name=str_remove(modisco_id, fixed("patterns.")) %>%
           str_remove("[0-9]+$") %>% str_c(name_seqlet)) %>%
  mutate(name=str_c(name_seqlet, ":", target_name)) %>%
  transmute(
    chrom1=chrom1_seqlet, start1=start1_seqlet, end1=end1_seqlet,
    chrom2=chrom2_target, start2=start2_target, end2=end2_target,
    name=name, score=score_target,
    strand1=strand1_seqlet, strand2=strand2_target,
    seqlet_score=score_seqlet, peak_score=target_score, sce2g_score=score_target,
    modisco_id=modisco_id, pattern_name=pattern_name,
    tf_top=name_seqlet, sce2g_id=name_target, target_name=target_name,
    enhancer_type=enhancer_type, ensembl_id=ensembl_id
  )

message("Summarizing links by modisco pattern...")
df_summary <- df_combined %>%
  group_by(pattern_name, tf_top) %>%
  arrange(target_name) %>%
  summarize(n_targets=length(unique(ensembl_id)),
            target_names=str_c(unique(target_name), collapse=";"),
            target_ids=str_c(unique(ensembl_id), collapse=";"),
            mean_sce2g_score=mean(sce2g_score),
            median_sce2g_score=median(sce2g_score),
            mean_seqlet_score=mean(seqlet_score),
            median_seqlet_score=median(seqlet_score),
            mean_peak_score=mean(peak_score),
            median_peak_score=median(peak_score),
            .groups='drop')

message("Adding info on non-linked seqlet counts...")
df_summary <- df_seqlet %>%
  mutate(has_link=idx_seqlet %in% df_overlap$idx_seqlet,
         pattern_name=str_remove(modisco_id, fixed("patterns.")) %>%
           str_remove("[0-9]+$") %>% str_c(name)) %>%
  group_by(pattern_name) %>%
  summarize(n_seqlets=dplyr::n(),
            n_linked=sum(has_link),
            frac_linked=mean(has_link)) %>%
  full_join(df_summary, by="pattern_name") %>%
  arrange(-n_seqlets)

message("Writing seqlet-target BEDPE...")
write_tsv(df_combined, snakemake@output$bedpe, col_names=FALSE)
message("Writing modisco pattern stats...")
write_csv(df_summary, snakemake@output$csv)

message("Done.")
