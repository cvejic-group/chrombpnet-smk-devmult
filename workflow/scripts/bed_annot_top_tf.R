library(tidyverse)
library(magrittr)
library(plyranges)
library(jsonlite)

df_hm <- stream_in(file(snakemake@params$json))

df_map <- tibble(
  hocomoco_id=df_hm$name,
  gene_symbol=df_hm$masterlist_info$species$HUMAN$gene_symbol,
  uniprot_id=df_hm$masterlist_info$species$HUMAN$uniprot_id,
  uniprot_ac=df_hm$masterlist_info$species$HUMAN$uniprot_ac,
  tf_name=df_hm$masterlist_info$tf,
  tf_superclass=df_hm$masterlist_info$tfclass_superclass,
  tf_class=df_hm$masterlist_info$tfclass_class,
  tf_family=df_hm$masterlist_info$tfclass_family,
  tf_subfamily=df_hm$masterlist_info$tfclass_subfamily,
  entrez=df_hm$masterlist_info$species$HUMAN$entrez
)

df_motifs <- read_csv(snakemake@input$csv)

df_tfs <- df_motifs %>%
  pivot_longer(cols=starts_with("match") | starts_with("qval"),
               names_to=c(".value", "match_rank"),
               names_pattern="(match|qval)(\\d+)") %>%
  filter(!is.nan(qval) | match_rank == 0) %>%
  left_join(df_map, by=c("match"="hocomoco_id"))

df_top <- df_tfs %>%
  group_by(cluster, pattern) %>%
  mutate(qratio=qval/qval[match_rank == 0]) %>%
  summarize(tf_names_top=tf_name[match_rank == 0],
            tf_names_2=str_c(sort(tf_name[qratio < 2]), collapse="|"),
            .groups='drop')

pat2tf_hsc <- df_top %>% pull(name=pattern, var=tf_names_top)


read_bed(snakemake@input$bed) %>%
  mutate(pattern=str_extract(name, "^(neg|pos)_patterns\\.[^.]+"),
         tf=pat2tf_hsc[pattern],
         modisco_id=name,
         name=tf) %>%
  write_bed(str_remove(snakemake@output$bed,  ".bgz$"), index=TRUE)
