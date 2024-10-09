# Load deltas


# Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(289)


export_dir <- "presentations/2024-03_rerun/"

neuron_properties <- read_csv("data/neuron_properties.csv") |>
  mutate(Strata = factor(Strata))








#~ Load ----

#> slow, skip to "Load precomputed"

#~ deltapsi ----
data_dir <- "data/2024-03-04_outs/deltapsi/"

files_dpsi <- list.files(data_dir,
                         pattern = "\\.tsv$",
                         full.names = FALSE)

dpsidta <- map_dfr(files_dpsi,
                  ~read_tsv(file.path(data_dir, .x),
                            col_names = c("gene_id", "lsv_id", "lsv_type",
                                          "dpsi", "p20", "p05", "psiA", "psiB",
                                          "nb_sj", "nb_exons",
                                          "sj_coords", "ir_coords"),
                            col_types = cols(
                              gene_id = col_character(),
                              lsv_id = col_character(),
                              lsv_type = col_character(),
                              dpsi = col_character(),
                              p20 = col_character(),
                              p05 = col_character(),
                              psiA = col_character(),
                              psiB = col_character(),
                              nb_sj = col_double(),
                              nb_exons = col_double(),
                              sj_coords = col_character(),
                              ir_coords = col_character()
                            ),
                            skip = 1,
                            na = "na",
                            progress = FALSE) |>
                    add_column(neurA = str_split(.x,"[-\\.]")[[1]][1],
                               neurB = str_split(.x,"[-\\.]")[[1]][2]),
                  .progress = TRUE)

# neurs_here <- neuron_properties |>
#   filter(include == "yes") |>
#   pull(Neuron_type)
# 
# dpsidta <- dpsidta |>
#   filter(neurA %in% neurs_here,
#          neurB %in% neurs_here)

# separate data in E(PSI) and Std(PSI) fields
# We get one row per junction
dpsi <- dpsidta |>
  separate_rows(dpsi, p20, p05, psiA, psiB,
                sep = ";") |>
  group_by(neurA, neurB, lsv_id) |>
  mutate(junction_id = row_number()) |>
  ungroup()  |>
  mutate(across(c(dpsi,p20, p05, psiA, psiB), as.double),
         junction_id = factor(junction_id))

# qs::qsave(dpsi, "intermediates/240304_dpsi.qs")





#~ Load PSI ----
data_dir <- "data/2024-03-04_outs/psi/"

files_psi <- list.files(data_dir,
                         pattern = "\\.tsv$",
                         full.names = FALSE)

# xx <- read_tsv(file.path(data_dir, files_psi[[2]]),na = "na")

psidta <- map_dfr(files_psi,
                   ~read_tsv(file.path(data_dir, .x),
                             col_names = c("gene_id", "lsv_id", "lsv_type",
                                           "mean_psi_per_lsv_junction", "stdev_psi_per_lsv_junction",
                                           "nb_sj", "nb_exons",
                                           "sj_coords", "ir_coords"),
                             col_types = cols(
                               gene_id = col_character(),
                               lsv_id = col_character(),
                               lsv_type = col_character(),
                               mean_psi_per_lsv_junction = col_character(),
                               stdev_psi_per_lsv_junction = col_character(),
                               nb_sj = col_double(),
                               nb_exons = col_double(),
                               sj_coords = col_character(),
                               ir_coords = col_character()
                             ),
                             skip = 1,
                             na = "na",
                             progress = FALSE) |>
                     add_column(neuron = str_split_1(.x, "\\.")[[1]]),
                  .progress = TRUE)

psi <- psidta |>
  separate_rows(mean_psi_per_lsv_junction, stdev_psi_per_lsv_junction,
                sep = ";") |>
  group_by(neuron, lsv_id) |>
  mutate(junction_id = row_number()) |>
  ungroup()  |>
  mutate(across(c(mean_psi_per_lsv_junction,stdev_psi_per_lsv_junction), as.double),
         junction_id = factor(junction_id))
# qs::qsave(psi, "intermediates/240304_psi.qs")





# >>>>> Load precomputed <<<<< ----

psi <- qs::qread("intermediates/240304_psi.qs")
dpsi <- qs::qread("intermediates/240304_dpsi.qs")
all_neurs_sequenced <- unique(c(dpsi$neurA, dpsi$neurB))

stopifnot(identical(sort(all_neurs_sequenced), sort(unique(psi$neuron))))





#~ genes expressed ----

# use CeNGEN scRNA-Seq, threshold 3
gene_expression <- 1L * (cengenDataSC::cengen_sc_3_bulk > 0)

neurs_with_known_expr <- colnames(gene_expression)

stopifnot(all(all_neurs_sequenced %in% neurs_with_known_expr ))




#~ filter dpsi based on expression ----
gene_expr_tib <- gene_expression |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expressed")

# filter to keep only genes that are expressed in this neur type
dpsi_filt <- dpsi |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron_id")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron_id")) |>
  filter(expressed.x > 0, expressed.y > 0) |>
  select(-expressed.x, -expressed.y)

psi_filt <- psi |>
  left_join(gene_expr_tib,
            by = c("gene_id", neuron = "neuron_id")) |>
  filter(expressed > 0) |>
  select(-expressed)



# DEGs from Alec's integrations
nb_degs_per_pair <- read.delim("data/2021-11-30_alec_integration/Total_integrated_DEGS_pairwise_113021.tsv") |>
  as_tibble() |>
  rowwise() |>
  mutate(pair = c_across(starts_with("cell_")) |> sort() |> paste0(collapse = "-"))

neurs_with_deg_data <- unique(union(nb_degs_per_pair$cell_A,
                                    nb_degs_per_pair$cell_B))


genes_degs <- readRDS("data/genes/integrated_significant_genes_pairwise_bsn9_012722.rds") |>
  map(as_tibble, rownames = "gene_id") |>
  imap_dfr(~add_column(.x, pair = .y)) |>
  mutate(pair = str_remove_all(pair, "[()]"))

strsplit(genes_degs$pair, "-") |>
  unlist() |>
  unique() |>
  sort() |>
  all.equal(neurs_with_deg_data)





# DS genes  ----


# ds_genes <- dpsi |>
#   filter(p20 >.50 & p05 < .05) |>
#   pull(gene_id) |>
#   unique()


dsf_genes <- dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique()







# Nb of genes with DS ----

# nb genes tested

# length(unique(dpsi$gene_id))
length(unique(dpsi_filt$gene_id))

# total DS
dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique() |>
  length()






# Compare literature ----

bib_all <- read_tsv("data/biblio/bib_ds_genes.tsv")$gene_id



genes_measurable <- dpsi_filt$gene_id |> unique()

plot(eulerr::euler(list(literature = bib_all,
                        `this work` = genes_measurable)),
     quantities = TRUE,
     main = "Number of detectable genes")


bib_all <- intersect(bib_all, genes_measurable)




(litt_plot <- plot(eulerr::euler(list(literature = bib_all,
                                      `this work` = dsf_genes)),
                   quantities = TRUE,
                   main = "Number of dAS genes"))

# ggsave("compare_literature_filt.pdf", path = export_dir, plot = litt_plot,
#        width = 10, height = 6, units = "cm")

#> Fig. 5A

# --> export source data

# writexl::write_xlsx(list(literature = data.frame(bib_all),
#                          this_work = data.frame(dsf_genes)),
#                     file.path(export_dir, "fig5A_sourceData.xlsx"))



# Gene type ----

#~ GO with background ----


# Restrict to genes expressed in more than 2 neurons in our dataset

nb_neurs_where_gene_expr <- rowSums(gene_expression[,all_neurs_sequenced])
genes_in_our_neur_sample <- names(nb_neurs_where_gene_expr)[nb_neurs_where_gene_expr > 2]


# writeLines(genes_in_our_neur_sample, "intermediates/240304_background_genes.txt")

writeClipboard(dsf_genes[dsf_genes %in% genes_in_our_neur_sample])
# -> use Wormbase enrichment analysis

#> Fig. S1A


#~ Gene families ----

bind_rows(wormDatasets::gene_families |>
            select(gene_id, gene_name, family),
          wormDatasets::list_rbp_tamburino2013 |>
              select(gene_id, gene_name) |>
              mutate(family = "RBP")) |>
  filter(gene_id %in% genes_measurable) |>
  mutate(family = recode(
    family,
    DEG_ENaC = "DEG/ENaC channel",
    downstream_GPCR = "GPCR signaling",
    GPCR = "GPCR signaling",
    neuropept_metabo = "neuropeptide signaling",
    nt_degradation = "neuropeptide signaling",
    nt_synthesis = "neuropeptide signaling",
    nt_transporter = "neuropeptide signaling",
    potassium_channel = "potassium channel",
    ribosome = "ribosome subunit",
    RBP = "RNA-binding protein",
    synaptic_ves = "synaptic vesicle",
    TF = "transcription factor",
    trp_channel = "TRP channel"
  )) |>
  mutate(is_ds = gene_id %in% dsf_genes) |>
  summarize(nb_ds = sum(is_ds),
            nb_tot = n(),
            prop_ds = nb_ds / nb_tot,
            .by = "family") |>
  mutate(p_higher = phyper(nb_ds, nb_tot,
                    sum(nb_tot) - nb_ds,sum(nb_ds), lower.tail = FALSE),
         p_lower = phyper(nb_ds, nb_tot,
                          sum(nb_tot) - nb_ds,sum(nb_ds), lower.tail = TRUE),
         padj_higher = p.adjust(p_higher, "fdr", n = 2*length(p_higher)),
         padj_lower = p.adjust(p_lower, "fdr", n = 2*length(p_higher)),
         sig_higher = case_when(
           padj_higher >= .05 ~ "",
           padj_higher < .001 ~ "***",
           padj_higher < .01 ~ "**",
           padj_higher < .05 ~ "*"
         ),
         sig_lower = case_when(
           padj_lower >= .05 ~ "",
           padj_lower < .001 ~ "###",
           padj_lower < .01 ~ "##",
           padj_lower < .05 ~ "#"
         )) |>
  arrange(desc(prop_ds)) |> mutate(family = fct_inorder(family)) |>
  ggplot() +
  theme_classic() +
  xlab(NULL) + ylab("Proportion dAS genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        
        legend.position = "none") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0,1)) +
  geom_col(aes(x = family, y = prop_ds)) +
  geom_text(aes(x = family, y = 0,
                label = nb_tot),
            nudge_y = .01,
            hjust = 0,
            angle = 90,
            color = 'white') +
  geom_text(aes(x = family, y = prop_ds, label = sig_higher),
            nudge_y = .03) +
  geom_text(aes(x = family, y = prop_ds, label = sig_lower),
            nudge_y = .09)

# ggsave("ds_by_family.pdf", path = export_dir,
#        width = 10, height = 12, units = "cm")

#> Fig. 5B

# --> source data

# bind_rows(wormDatasets::gene_families |>
#             select(gene_id, gene_name, family),
#           wormDatasets::list_rbp_tamburino2013 |>
#             select(gene_id, gene_name) |>
#             mutate(family = "RBP")) |>
#   filter(gene_id %in% genes_measurable) |>
#   mutate(family = recode(
#     family,
#     DEG_ENaC = "DEG/ENaC channel",
#     downstream_GPCR = "GPCR signaling",
#     GPCR = "GPCR signaling",
#     neuropept_metabo = "neuropeptide signaling",
#     nt_degradation = "neuropeptide signaling",
#     nt_synthesis = "neuropeptide signaling",
#     nt_transporter = "neuropeptide signaling",
#     potassium_channel = "potassium channel",
#     ribosome = "ribosome subunit",
#     RBP = "RNA-binding protein",
#     synaptic_ves = "synaptic vesicle",
#     TF = "transcription factor",
#     trp_channel = "TRP channel"
#   )) |>
#   mutate(is_ds = gene_id %in% dsf_genes) |>
#   writexl::write_xlsx(file.path(export_dir, "sourceData_families.xlsx"))


# bind_rows(wormDatasets::gene_families |>
#             select(gene_id, gene_name, family),
#           wormDatasets::list_rbp_tamburino2013 |>
#             select(gene_id, gene_name) |>
#             mutate(family = "RBP")) |>
#   filter(gene_id %in% genes_measurable) |>
#   mutate(family = recode(
#     family,
#     DEG_ENaC = "DEG/ENaC channel",
#     downstream_GPCR = "GPCR signaling",
#     GPCR = "GPCR signaling",
#     neuropept_metabo = "neuropeptide signaling",
#     nt_degradation = "neuropeptide signaling",
#     nt_synthesis = "neuropeptide signaling",
#     nt_transporter = "neuropeptide signaling",
#     potassium_channel = "potassium channel",
#     ribosome = "ribosome subunit",
#     RBP = "RNA-binding protein",
#     synaptic_ves = "synaptic vesicle",
#     TF = "transcription factor",
#     trp_channel = "TRP channel"
#   )) |>
#   mutate(is_ds = gene_id %in% dsf_genes) |>
#   summarize(nb_ds = sum(is_ds),
#             nb_tot = n(),
#             prop_ds = nb_ds / nb_tot,
#             .by = "family") |>
#   mutate(p_higher = phyper(nb_ds, nb_tot,
#                            sum(nb_tot) - nb_ds,sum(nb_ds), lower.tail = FALSE),
#          p_lower = phyper(nb_ds, nb_tot,
#                           sum(nb_tot) - nb_ds,sum(nb_ds), lower.tail = TRUE),
#          padj_higher = p.adjust(p_higher, "fdr", n = 2*length(p_higher)),
#          padj_lower = p.adjust(p_lower, "fdr", n = 2*length(p_higher)),
#          sig_higher = case_when(
#            padj_higher >= .05 ~ "",
#            padj_higher < .001 ~ "***",
#            padj_higher < .01 ~ "**",
#            padj_higher < .05 ~ "*"
#          ),
#          sig_lower = case_when(
#            padj_lower >= .05 ~ "",
#            padj_lower < .001 ~ "###",
#            padj_lower < .01 ~ "##",
#            padj_lower < .05 ~ "#"
#          )) |>
#   arrange(desc(prop_ds)) |> mutate(family = fct_inorder(family)) |>
#   writexl::write_xlsx(file.path(export_dir, "sourceData_families_byfamily.xlsx"))



# Heatmap ----


annot_df <- neuron_properties |>
  select(-include) |>
  column_to_rownames("Neuron_type")


gene_expr_tib <- gene_expression |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron",
               values_to = "expressed")

# find nb of coexpressed genes which are DS
coexpr_ds <- dpsi_filt |>
  filter(neurA %in% neurs_with_known_expr,
         neurB %in% neurs_with_known_expr) |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  summarize(has_ds = any(ds),
            .by = c(neurA, neurB, gene_id)) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, has_ds, coexpressed) |>
  group_by(neurA, neurB) |>
  summarize(nb_ds = sum(has_ds),
            nb_coexpr = sum(coexpressed),
            nb_both = sum(has_ds & coexpressed),
            nb_total = n(),
            .groups = "drop")

# make matrix of proportion coexpr genes that are DS
hm_coexpr <- coexpr_ds |>
  mutate(prop_ds = 100*nb_both/nb_coexpr) |>
  select(neurA,neurB,prop_ds)


ds_mat <- hm_coexpr |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(hm_coexpr) |>
  pivot_wider(names_from = neurB,
              values_from = prop_ds) |>
  arrange(neurA) |>
  (\(.x) select(.x, order(colnames(.x))))() |>
  column_to_rownames("neurA") |> 
  as.matrix()


# Show neurons clustered by prop of DS
dist_mat_prop <- ds_mat/max(ds_mat, na.rm = TRUE)
hc <- hclust(as.dist(dist_mat_prop), method = "ward.D")
# plot(hc)

pheatmap::pheatmap(ds_mat,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name =
                                                                       "Blues"))(151),
                   scale = "none",
                   breaks = (0:150)/5,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # cutree_rows = 4,
                   # cutree_cols = 4,
                   # main = "Proportion of DS",
                   legend_breaks = c(0,5,10,15,20, 25, 30),
                   legend_labels = c("0%","5%","10%","15%","20%","25%","30%"),
                   # annotation_row = annot_df,
                   annotation_col = annot_df,
                   border_color = NA,
                   # filename = file.path(export_dir, "heatmap_ds_big.pdf"),
                   width = 14,
                   height = 11
)

dev.off()

#> Fig. 5C


# --> source data
# writexl::write_xlsx(ds_mat |> as.data.frame() |> rownames_to_column("neuron"),
#                     file.path(export_dir, "sourceData_heatmap.xlsx"))





# Show neurons clustered by nb of DS
nb_ds_mat <-  coexpr_ds |>
  select(neurA = neurB, neurB = neurA, nb_both) |>
  bind_rows(coexpr_ds |> select(neurA, neurB, nb_both)) |>
  pivot_wider(names_from = neurB,
              values_from = nb_both) |>
  arrange(neurA) |>
  (\(.x) select(.x, order(colnames(.x))))() |>
  column_to_rownames("neurA") |>
  as.matrix()


dist_mat_nb <- nb_ds_mat/max(nb_ds_mat, na.rm = TRUE)
heatmap(dist_mat_nb, Rowv = NA, Colv = NA)
plot(hclust(as.dist(dist_mat_nb), method = "complete"))





# Show neurons clustered by similarity of profile

ds_mat |> 
  dist(method = "canberra") |>
  hclust(method = "complete") |>
  plot()

nb_ds_mat |> 
  dist(method = "canberra") |>
  hclust(method = "complete") |>
  plot()


# Pairs of similar neurons ----

source("R/utils.R")


plot_hc_clustered(hc, h = .2, cex = .7,
                  xlab = "Neuron type", sub = NA, main = "Neuron similarity (by proportion of genes dAS)")
abline(h = .2, lty = 'dashed', col = 'grey80')

# pdf(file.path(export_dir, "dendrogram_similarity.pdf"),
#     width = 9, height = 4)
# # 
# plot_hc_clustered(hc, h = .2, cex = .9,
#                   xlab = "Neuron type", sub = NA, main = "Neuron similarity (by proportion of genes dAS)")
# abline(h = .2, lty = 'dashed', col = 'grey80')
# dev.off()


#> Fig. S1B

# --> source data: same as heatmap






# # Event broadness ----
# Probably wrong approach: we would prefer to work on PSI, not dPSI
# as it's more interesting to say "this spliceform is un/common" than to
# look at a number of pairs (could be 1 single neuron that differs)
# 
# lsv_dpsi <- dpsi |>
#   mutate(is_ds_jct = (p20 >.50 & p05 < .05)) |>
#   group_by(lsv_id, neurA, neurB) |>
#   summarize(is_ds_lsv = any(is_ds_jct),
#             .groups = 'drop')
# 
# 
# lsv_ds <- lsv_dpsi |>
#   summarize(nb_times_ds = sum(is_ds_lsv),
#             n = length(is_ds_lsv),
#             .by = "lsv_id") |>
#   mutate(prop_ds = nb_times_ds/n)
# 
# hist(lsv_ds$n, breaks = 100)
# hist(lsv_ds$prop_ds, breaks = 100)
# hist(lsv_ds$nb_times_ds, breaks = 100)
# 
# lsv_ds |>
#   mutate(category = case_when(
#     n <= 2 ~ "unmeasured",
#     nb_times_ds == 0 ~ "not_ds",
#     nb_times_ds <= 3 ~ "rare",
#     .default = "common"
#   ))
# 
# 
# dpsi |>
#   filter(lsv_id == "WBGene00000006:s:2573999-2574267") |>
#   View()
# 
# 
# 
# dpsi |>
#   mutate(is_ds = (p20 >.50 & p05 < .05)) |>
#   group_by(neurA, neurB) |>
#   summarize(nb_times_ds = sum(is_ds),
#             n = length(is_ds))




# Dominant spliceform ----

#> no used in paper


#~ by lsv ----

dominant_jct_pattern <- function(jct_id, threshold = .7){
  
  tab <- sort(table(jct_id), decreasing = TRUE)
  
  stopifnot(sum(tab) == length(jct_id))
  prop_tab <- tab/sum(tab)
  
  if(prop_tab[[1]] >= threshold){
    return("single dominant form")
    
  } else if(sum(prop_tab[1:2]) >= threshold){
    return("two main forms")
    
  } else if(sum(prop_tab[1:3]) >= threshold){
    return("three main forms")
    
  } else{
    return("complex")
  }
  
}




dominant_jct <- psi |>
  filter(mean_psi_per_lsv_junction > .3) |>
  summarize(jct_pattern = dominant_jct_pattern(junction_id) |>
              factor(levels = c("single dominant form", "two main forms", "three main forms", "complex")),
            gene_id = gene_id[[1]],
            .by = c(lsv_id))





dominant_jct$jct_pattern |> table()

gg_dominant_jct <- dominant_jct |>
  mutate(jct_pattern = if_else(jct_pattern == "three main forms", "complex", jct_pattern) |>
           factor(levels = c("single dominant form", "two main forms", "complex"))) |>
  mutate(jct_pattern = recode_factor(
    jct_pattern,
    `single dominant form` = "single\ndominant\nform",
    `two main forms` = "two\nmain\nforms"
  )) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = jct_pattern)) +
  # scale_y_log10() +
  ylab("Number of LSVs") +
  xlab("Spliceform expression pattern")

# ggsave("dominant_spliceform_by_lsv.pdf", path = export_dir,
#        width = 9, height = 12, units = "cm")


#~ by gene ----

dominant_gene_pattern <- function(jct_pat){
  
  jct_pat <- as.integer(jct_pat)
  
  nb_tot <- length(jct_pat)
  nb_single <- sum(jct_pat == 1L)
  nb_two <- sum(jct_pat == 2L)
  nb_three <- sum(jct_pat == 3L)
  nb_more <- sum(jct_pat > 3L)
  
  if(nb_single == nb_tot){
    return("single dominant form")
    
  } else if(nb_two == 1 && (nb_two + nb_single) == nb_tot){
    return("two main forms")
  } else if(nb_three == 1 && (nb_two + nb_single + nb_three) == nb_tot){
    return("three main forms")
  } else{
    return("complex")
  }
}


dominant_jct_gene <- dominant_jct |>
  summarize(gene_pattern = dominant_gene_pattern(jct_pattern) |>
              factor(levels = c("single dominant form", "two main forms", "three main forms", "complex")),
            .by = gene_id)



dominant_jct_gene$gene_pattern |> table()

gg_dominant_jct_gene <- dominant_jct_gene |>
  mutate(gene_pattern = if_else(gene_pattern == "three main forms", "complex", gene_pattern) |>
           factor(levels = c("single dominant form", "two main forms", "complex"))) |>
  mutate(gene_pattern = recode_factor(
    gene_pattern,
    `single dominant form` = "single\ndominant\nform",
    `two main forms` = "two\nmain\nforms"
  )) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = gene_pattern)) +
  # scale_y_log10() +
  ylab("Number of genes") +
  xlab("Spliceform expression pattern")

# ggsave("dominant_spliceform_by_gene.pdf", path = export_dir,
#        width = 8, height = 12, units = "cm")

patchwork::wrap_plots(gg_dominant_jct, gg_dominant_jct_gene)

ggsave("dominant_spliceforms.pdf", path = export_dir,
       width = 16, height = 10, units = "cm")








# DE vs DS ----

#~ neuron level ----




#~ By neur classes ----
nb_ds_genes_by_pair <- dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_DS_genes = n(),
            .groups = 'drop') |>
  filter(neurA %in% neurs_with_deg_data,
         neurB %in% neurs_with_deg_data) |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |> sort() |> paste0(collapse = "-"))


de_ds <- full_join(nb_degs_per_pair,
                   nb_ds_genes_by_pair,
                   by = "pair")

ggplot(de_ds, aes(x = total_integrated_DEGs, y = nb_DS_genes)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_text_repel(aes(label = pair),
                           max.overlaps = 7, nudge_x = .1) +
  scale_x_log10() + scale_y_log10() +
  xlab("Number of DE genes (log)") +
  ylab("Number of DAS genes (log)")

# ggsave("DS_vs_DE.pdf", path = export_dir,
#        width = 9, height = 7, units = "cm")

#> Fig. 5D

# ---> source data

# writexl::write_xlsx(de_ds,
#                     file.path(export_dir, "sourceData_das_vs_deg.xlsx"))




mod <- lm(log(nb_DS_genes)~log(total_integrated_DEGs), data = de_ds)
mod <- lm(nb_DS_genes~total_integrated_DEGs, data = de_ds)

plot(fitted.values(mod), residuals(mod))
summary(mod)
qqnorm(residuals(mod)); qqline(residuals(mod))




#~ Top DEGs vs top DSGs ----


top_degs <- genes_degs |>
  filter(p.hmp < .05) |>
  group_by(pair) |>
  slice_max(order_by = abs(Average_logFC),
            n = 100) |>
  summarize(top_DEGs = list(as.character(gene_id)))


# check on scRNA-Seq app:
# top_degs$top_DEGs[[150]] |> clipr::write_clip()
# top_degs$pair[[150]]

top_dsgs <- dpsi |>
  filter(neurA %in% neurs_with_deg_data,
         neurB %in% neurs_with_deg_data) |>
  filter(p20 >.50 & p05 < .05) |>
  summarize(max_dpsi = max(abs(dpsi)),
            .by = c(neurA, neurB, gene_id)) |>
  group_by(neurA, neurB) |>
  slice_max(order_by = max_dpsi,
            n = 100) |>
  summarize(top_DSGs = list(as.character(gene_id)),
            .groups = 'drop') |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |>
           sort() |>
           paste0(collapse = "-")) |>
  select(-neurA, -neurB)

sapply(top_degs$top_DEGs, length) |> hist()
sapply(top_dsgs$top_DSGs, length) |> hist()
table(lengths(top_dsgs$top_DSGs) == 100)

tops <- left_join(top_degs, top_dsgs,
          by = "pair") |>
  filter(lengths(top_DEGs) == 100,
         lengths(top_DSGs) == 100) |>
  rowwise() |>
  mutate(overlap = length(intersect(top_DEGs, top_DSGs)),
         nb_degs = length(top_DEGs),
         nb_dsgs = length(top_DSGs))

table(tops$overlap)


empty_eulerr_degs_dsgs <- eulerr::euler(list(`Top 100 DE genes` = 1:5,
                   `Top 100 dAS genes` = 5:9)) |>
  plot()

# ggsave("empty_eulerr_degs_dsgs.pdf", plot = empty_eulerr_degs_dsgs,
#        path = export_dir,
#        width = 7, height = 5, units = "cm")

ggplot(tops) +
  theme_classic() +
  geom_bar(aes(x = overlap)) +
  scale_x_continuous(breaks = 0:7) +
  xlab("Number of overlapping genes") +
  ylab("Number of neuron pairs")

# ggsave("hist_degs_dsgs.pdf", path = export_dir,
#        width = 9, height = 6, units = "cm")

#> Fig. S1C

# ---> source data

# writexl::write_xlsx(tops,
#                     file.path(export_dir, "sourceData_topDEGs_DAS.xlsx"))




#~ Only RBPs ----

# not in paper

rbps <- wormDatasets::list_rbp_tamburino2013 |>
  filter(gene_id %in% rownames(cengenDataSC::cengen_proportion_bulk)) |>
  pull(gene_id)

# genes_degs <- readRDS("data/genes/integrated_significant_genes_pairwise_bsn9_012722.rds") |>
#   map(as_tibble, rownames = "gene_id") %>%
#   map2_dfr(., names(.), ~add_column(.x, pair = .y))

nb_rbps_deg <- genes_degs |>
  filter(gene_id %in% rbps) |>
  group_by(pair) |>
  summarize(nb_degs = n()) |>
  mutate(pair = str_remove_all(pair, "[()]"))

stopifnot(sum(nb_rbps_deg$pair %in% nb_ds_genes_by_pair$pair) == sum(nb_ds_genes_by_pair$pair %in% nb_rbps_deg$pair))
stopifnot(sum(nb_rbps_deg$pair %in% nb_ds_genes_by_pair$pair) == nrow(nb_rbps_deg))
stopifnot(sum(nb_rbps_deg$pair %in% nb_ds_genes_by_pair$pair) == nrow(nb_ds_genes_by_pair))

de_ds_rbps <- full_join(nb_rbps_deg, nb_ds_genes_by_pair, by = "pair")


ggplot(de_ds_rbps, aes(x = nb_degs, y = nb_DS_genes)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  # ggrepel::geom_text_repel(aes(label = pair)) +
  scale_x_log10() + scale_y_log10() +
  xlab("Number of DE RBPs (log)") +
  ylab("Number of DAS genes (log)")
# ggsave("DS_vs_DE_rbps.pdf", path = export_dir,
#        width = 4, height = 4, units = "in")

mod <- lm(log(nb_DS_genes)~log(nb_degs), data = de_ds_rbps)

plot(fitted.values(mod), residuals(mod))
summary(mod)
qqnorm(residuals(mod));qqline(residuals(mod))





# Subsample neurons ----

# How many neurons do we need to sequence to detect that many DS genes?
# Subsamples



#~ single replicate ----
nb_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
}





#~ 10 replicates per subsample ----

#~~ DAS ----
sub_nb_genes_ds_reps <- expand_grid(nb_neurs = 1:length(all_neurs_sequenced),
                                    rep = 1:10) |>
  mutate(nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub,
                               .progress = TRUE))

# qs::qsave(sub_nb_genes_ds_reps, "intermediates/240502_sub_nb_genes_ds_reps.qs")
sub_nb_genes_ds_reps <- qs::qread("intermediates/240502_sub_nb_genes_ds_reps.qs")



#~~ same with gene expr ----

nb_expr_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  gene_expression[,neurset] |> matrixStats::rowAnys() |> sum()
}

sub_nb_genes_expr_reps <- expand_grid(nb_neurs = 2:length(all_neurs_sequenced),
                                      rep = 1:10) |>
  mutate(nb_expr_genes = map_int(nb_neurs, nb_expr_genes_sub,
                                 .progress = TRUE))


ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs))+
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DAS")

ggplot(sub_nb_genes_expr_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_expr_genes, group = nb_neurs))+
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DEG")






#~ Fit on log-linear ----

#~~ DAS ----
mod_lm <- lm(nb_ds_genes ~ log(nb_neurs),
             data = sub_nb_genes_ds_reps)
plot(fitted.values(mod_lm), residuals(mod_lm))
plot(fitted.values(mod_lm)^2, residuals(mod_lm))
summary(mod_lm)
qqnorm(residuals(mod_lm)); qqline(residuals(mod_lm))



sub_nb_genes_ds_reps$fit_lm <- predict(mod_lm)

projection <- tibble(
  nb_neurs = 1:119,
  prediction = predict(mod_lm, newdata = data.frame(nb_neurs=nb_neurs))
)

ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons sampled") +
  ylab("Number of dAS genes detected") +
  scale_x_continuous(breaks = c(0,20,40,60,80,100,119)) +
  scale_y_continuous(limits = c(0,NA),
                     labels = scales::label_comma(),
                     breaks = c(1000,2000,3000,3192)) +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=prediction), color = "red3", linetype = "dashed",
            data = projection) +
  geom_vline(aes(xintercept = 119), linetype = "dotted", color = "grey") +
  geom_hline(aes(yintercept = predict(mod_lm,
                                      newdata = data.frame(nb_neurs = 119))),
             linetype = "dotted", color = "grey")
# ggsave("nb_genes_das.pdf", path = export_dir,
#        width = 11, height = 7, units = "cm")

predict(mod_lm, newdata = data.frame(nb_neurs = 119))

# > Fig 5E

# ---> source data
# sub_nb_genes_ds_reps |>
#   writexl::write_xlsx(file.path(export_dir, "sourceData_subsample_DAS.xlsx"))





#~~ expr ----

gmod_lm <- lm(nb_expr_genes ~ log(nb_neurs),
              data = sub_nb_genes_expr_reps)
plot(fitted.values(gmod_lm), residuals(gmod_lm))
summary(gmod_lm)
predict(gmod_lm, newdata = data.frame(nb_neurs = 119))
sub_nb_genes_expr_reps$gmod_lm <- predict(gmod_lm, newdata = data.frame(nb_neurs=sub_nb_genes_expr_reps$nb_neurs))

ggplot(sub_nb_genes_expr_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_expr_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=gmod_lm), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  # geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as expressed")



# Prop
predict(mod_lm, newdata = data.frame(nb_neurs = 119))/predict(gmod_lm, newdata = data.frame(nb_neurs = 119))




#~ micment ----

#~~ DAS ----


dmod_micment <- nls(nb_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                   data = sub_nb_genes_ds_reps,
                   start = list(Gmax = 3000, b = 10))
summary(dmod_micment)
plot(fitted.values(dmod_micment), residuals(dmod_micment))
predict(dmod_micment, newdata = data.frame(nb_neurs = 119))


sub_nb_genes_ds_reps <- cbind(sub_nb_genes_ds_reps, fit=predict(dmod_micment))

ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  geom_hline(aes(yintercept = coef(dmod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DS")


#~~ expr ----



gmod_micment <- nls(nb_expr_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                    data = sub_nb_genes_expr_reps,
                    start = list(Gmax = 8000, b = 10))
summary(gmod_micment)
plot(fitted.values(gmod_micment), residuals(gmod_micment))
predict(gmod_micment, newdata = data.frame(nb_neurs = 119))

# Prop
predict(dmod_micment, newdata = data.frame(nb_neurs = 119))/predict(gmod_micment, newdata = data.frame(nb_neurs = 119))













