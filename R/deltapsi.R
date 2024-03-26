# Load deltas


# Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)


export_dir <- "presentations/2022-12_rerun/"

neuron_properties <- read_csv("data/neuron_properties.csv")





#~ Load ----
data_dir <- "data/2022-03-23_outs/deltapsi/"

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
                            na = "na") |>
                    add_column(neurA = str_split(.x,"[-\\.]")[[1]][1],
                               neurB = str_split(.x,"[-\\.]")[[1]][2]))

neurs_here <- neuron_properties |>
  filter(include == "yes") |>
  pull(Neuron_type)

dpsidta <- dpsidta |>
  filter(neurA %in% neurs_here,
         neurB %in% neurs_here)

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

# qs::qsave(dpsi, "intermediates/221214_dpsi_220323.qs")









dpsi <- qs::qread("intermediates/221214_dpsi_220323.qs")
all_neurs_sequenced <- unique(c(dpsi$neurA, dpsi$neurB))



#~ Alec genes expressed ----
# Use Alec's integrated GeTMMs to determine what genes are expressed in each neuron type
# Update 2022: use more recent data, already binary
gene_expr_int <- read.delim("data/genes/bsn9_subtracted_integrated_binarized_expression_withVDDD_FDR_0.1_092022.tsv")

neurs_integrated <- colnames(gene_expr_int)





ds_genes <- dpsi |>
       filter(p20 >.50 & p05 < .05) |>
       pull(gene_id) |> unique()

neur_genes <- wormDatasets::genes_by_pattern


table(ds_genes %in% neur_genes$present_in_neurons)

ds_genes |>
  intersect(neur_genes$present_in_neurons) |>
  head() |>
  i2s(gids)

ds_genes |>
  intersect(neur_genes$nonneuronal) |>
  head() |>
  i2s(gids)





# Nb of genes with DS ----

# nb genes tested
length(unique(dpsi$gene_id))

# total
dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique() |>
  length()

#~ By neur classes ----
nb_signif_genes_by_test <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_DS_genes = n())

nb_signif_genes_by_test |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_genes) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = neurA, y=neurB, fill = nb_DS_genes))





# Compare literature ----

bib_all <- read_tsv("data/biblio/bib_ds_genes.tsv")$gene_id


# Restrict to genes expressed in at least 2 neurons in our dataset
# expr_sc <- cengenDataSC::cengen_sc_3_bulk > 0


nb_neurs_where_gene_expr <- rowSums(gene_expr_int[,all_neurs_sequenced])
genes_in_our_neur_sample <- names(nb_neurs_where_gene_expr)[nb_neurs_where_gene_expr > 2]

table(bib_all %in% genes_in_our_neur_sample)
bib_all <- intersect(bib_all, genes_in_our_neur_sample)


all_signif_genes <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |>
  unique()

length(all_signif_genes)

(litt_plot <- plot(eulerr::euler(list(literature = bib_all,
                                      `this work` = all_signif_genes)),
                   quantities = TRUE))

# ggsave("compare_literature.png", path = export_dir, plot = litt_plot,
#        width = 4, height = 3, units = "in")




# GO with background ----

writeLines(genes_in_our_neur_sample, "intermediates/221215_background_genes.txt")
writeClipboard(all_signif_genes[all_signif_genes %in% genes_in_our_neur_sample])
# -> use Wormbase enrichment analysis



# Heatmap ----


gene_expr_tib <- gene_expr_int |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron",
               values_to = "expressed")

# find nb of coexpressed genes which are DS
coexpr_ds <- dpsi |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  select(gene_id, ds, neurA, neurB) |>
  filter(neurA %in% neurs_integrated,
         neurB %in% neurs_integrated) |>
  distinct() |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, ds, coexpressed) |>
  group_by(neurA, neurB) |>
  summarize(nb_ds = sum(ds),
            nb_coexpr = sum(coexpressed),
            nb_both = sum(ds & coexpressed),
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


pheatmap::pheatmap(ds_mat,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "Blues"))(100),
                   scale = "none",
                   breaks = (0:100)/5,
                   # cutree_rows = 2,
                   # cutree_cols = 2,
                   main = "Proportion of coexpressed genes DS",
                   # filename = file.path(export_dir, "heatmap_ds.png"),
                   width = 8,
                   height = 7
)


# Show neurons clustered by prop of DS
dist_mat_prop <- ds_mat/max(ds_mat, na.rm = TRUE)
heatmap(dist_mat_prop, Rowv = NA, Colv = NA)
plot(hclust(as.dist(dist_mat_prop), method = "complete"))

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





# Cluster analysis...
# need to redefine clusters first

table(cutree(hc2, k=2))
high_ds <- hc2$labels[cutree(hc2, k=2) == 1]
low_ds <- hc2$labels[cutree(hc2, k=2) == 2]

high_ds_mat <- ds_mat[high_ds, high_ds]
high_ds_mat[lower.tri(high_ds_mat)] <- NA
xx <- as.numeric(high_ds_mat)
xx <- xx[!is.na(xx)]
range(xx)
median(xx)


low_ds_mat <- ds_mat[low_ds, low_ds]
low_ds_mat[lower.tri(low_ds_mat)] <- NA
xx <- as.numeric(low_ds_mat)
xx <- xx[!is.na(xx)]
range(xx)
median(xx)




# DE vs DS ----

#~ neuron level ----
neurs_integrated_noD <- neurs_integrated |> setdiff(c("VD","DD"))

degs <- read.delim("data/2021-11-30_alec_integration/Total_integrated_DEGS_pairwise_113021.tsv") |>
  as_tibble() |>
  rowwise() |>
  mutate(pair = c_across(starts_with("cell_")) |> sort() |> paste0(collapse = "-"))

dsgs <- nb_signif_genes_by_test |>
  ungroup() |>
  filter(neurA %in% neurs_integrated_noD,
         neurB %in% neurs_integrated_noD) |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |> sort() |> paste0(collapse = "-"))


de_ds <- full_join(degs, dsgs, by = "pair")

ggplot(de_ds, aes(x = total_integrated_DEGs, y = nb_DS_genes)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_text_repel(aes(label = pair)) +
  scale_x_log10() + scale_y_log10() +
  xlab("Number of DE genes (log)") +
  ylab("Number of DAS genes (log)")
# ggsave("DS_vs_DE.pdf", path = export_dir,
#        width = 4, height = 4, units = "in")

mod <- lm(log(nb_DS_genes)~log(total_integrated_DEGs), data = de_ds)

plot(fitted.values(mod), residuals(mod))
summary(mod)
qqnorm(residuals(mod)); qqline(residuals(mod))



#~ Only RBPs ----

sf_walhout <- read_csv("data/biblio/list_rbp_walhout.csv",
                       skip=1,
                       col_types = cols(
                         `Gene Name` = col_character(),
                         `WBID (WS219)` = col_character(),
                         ORF = col_character(),
                         RBD = col_character(),
                         Group = col_double(),
                         Source = col_character(),
                         GO = col_character(),
                         UniProtKB = col_character()))

rbps <- sf_walhout$`WBID (WS219)`[sf_walhout$`WBID (WS219)` %in% rownames(cengenDataSC::cengen_proportion_bulk)]


genes_degs <- readRDS("data/genes/integrated_significant_genes_pairwise_bsn9_012722.rds") |>
  map(as_tibble, rownames = "gene_id") %>%
  map2_dfr(., names(.), ~add_column(.x, pair = .y))

nb_rbps_deg <- genes_degs |>
  filter(gene_id %in% rbps) |>
  group_by(pair) |>
  summarize(nb_degs = n()) |>
  mutate(pair = str_remove_all(pair, "[()]"))

stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == sum(dsgs$pair %in% nb_rbps_deg$pair))
stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == nrow(nb_rbps_deg))
stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == nrow(dsgs))

de_ds_rbps <- full_join(nb_rbps_deg, dsgs, by = "pair")


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

sub_nb_genes_ds <- tibble(nb_neurs = 1:length(all_neurs_sequenced),
                          nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub))



# Fit Michaelis-Menten model
mod_micment <- nls(nb_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                   data = sub_nb_genes_ds,
                   start = list(Gmax = 3000, b = 10))
summary(mod_micment)
plot(fitted.values(mod_micment), residuals(mod_micment))
qqnorm(residuals(mod_micment))
qqline(residuals(mod_micment))


sub_nb_genes_ds <- cbind(sub_nb_genes_ds, fit=predict(mod_micment))

ggplot(sub_nb_genes_ds) +
  theme_classic() +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DS")
# ggsave("nb_genes_ds.pdf", path = export_dir,
#        width = 5, height = 5/1.375)





#~ 10 replicates per subsample ----

#~~ DAS ----
sub_nb_genes_ds_reps <- expand_grid(nb_neurs = 1:length(all_neurs_sequenced),
                                    rep = 1:10) |>
  mutate(nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub))



#~~ expr ----

nb_expr_genes_sub <- function(n){
  neurset <- sample(neurs_integrated, n)
  
  gene_expr_int[,neurset] |> matrixStats::rowAnys() |> sum()
}

sub_nb_genes_expr_reps <- expand_grid(nb_neurs = 2:length(neurs_integrated),
                                      rep = 1:10) |>
  mutate(nb_expr_genes = map_int(nb_neurs, nb_expr_genes_sub))


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
summary(mod_lm)



sub_nb_genes_ds_reps$fit_lm <- predict(mod_lm)

ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=fit_lm), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DAS") +
  scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA), labels = scales::label_comma())
# ggsave("nb_genes_das.png", path = export_dir,
#        width = 5, height = 5/1.375)

predict(mod_lm, newdata = data.frame(nb_neurs = 119))



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
  ylab("Number of genes detected as DS")



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
# ggsave("nb_genes_ds_reps.pdf", path = export_dir,
#        width = 5.5, height = 3.5)


#~~ expr ----



gmod_micment <- nls(nb_expr_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                    data = sub_nb_genes_expr_reps,
                    start = list(Gmax = 8000, b = 10))
summary(gmod_micment)
plot(fitted.values(gmod_micment), residuals(gmod_micment))
predict(gmod_micment, newdata = data.frame(nb_neurs = 119))

# Prop
predict(dmod_micment, newdata = data.frame(nb_neurs = 119))/predict(gmod_micment, newdata = data.frame(nb_neurs = 119))







#~ By proportion of expressed genes ----
# Say a gene is expressed if threshold 2 from sc detects it in at least 2 of the sampled classes

prop_ds_genes_sub <- function(n){
  neurset <- sample(neurs_integrated, n)
  
  nb_neurs_where_gene_expr <- rowSums(gene_expr_int[,neurset])
  nb_genes_expr_in_neurset <- sum(nb_neurs_where_gene_expr > 1)
  
  
  nb_genes_ds_in_neurset <- dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
  
  nb_genes_ds_in_neurset/nb_genes_expr_in_neurset
}


gene_coexpr_by_neur_pair <- dpsi |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  select(gene_id, ds, neurA, neurB) |>
  filter(neurA %in% neurs_integrated,
         neurB %in% neurs_integrated) |>
  distinct() |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, ds, coexpressed) 


prop_ds_genes_sub <- function(n){
  neurset <- sample(neurs_integrated, n)
  
  gene_coexpr_by_neur_pair |>
    filter(neurA %in% neurset,
           neurB %in% neurset,
           coexpressed) |>
    group_by(neurA,neurB) |>
    summarize(prop_ds_in_pair = mean(ds),
              .groups = "drop") |>
    summarize(mean_prop_ds = mean(prop_ds_in_pair)) |>
    pull(mean_prop_ds)
}

sub_prop_reps <- expand_grid(nb_neurs = 3:length(neurs_integrated),
                        rep = 1:10) |>
  mutate(prop_ds_genes = map_dbl(nb_neurs, prop_ds_genes_sub))




ggplot(sub_prop_reps) +
  theme_classic() +
  # geom_boxplot(aes(x = nb_neurs, y = prop_ds_genes, group = nb_neurs))+
  # geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=prop_ds_genes)) +
  # geom_hline(aes(yintercept = coef(mod_prop_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  # scale_y_continuous(limits = c(0,.5), labels = \(xx) scales::percent(xx, accuracy = 1)) +
  xlab("Number of neurons sampled") +
  ylab("Mean proportion of coexpressed genes detected as DS")






#~ older version, using sc ----

# expr_sc <- cengenDataSC::cengen_sc_3_bulk > 0


prop_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  nb_neurs_where_gene_expr <- rowSums(gene_expr_int[,neurset])
  nb_genes_expr_in_neurset <- sum(nb_neurs_where_gene_expr > 1)
  
  
  nb_genes_ds_in_neurset <- dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
  
  nb_genes_ds_in_neurset/nb_genes_expr_in_neurset
}


sub_prop_genes_ds <- tibble(nb_neurs = 3:length(all_neurs_sequenced),
                            prop_ds_genes = map_dbl(nb_neurs, prop_ds_genes_sub))


# Fit Michaelis-Menten model
mod_micment <- nls(prop_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                   data = sub_prop_genes_ds,
                   start = list(Gmax = 3000, b = 10))

sub_prop_genes_ds <- cbind(sub_prop_genes_ds, fit=predict(mod_micment))

ggplot(sub_prop_genes_ds) +
  theme_classic() +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=prop_ds_genes)) +
  geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  scale_y_continuous(limits = c(0,.3), labels = \(xx) scales::percent(xx, accuracy = 1)) +
  xlab("Number of neurons sampled") +
  ylab("Proportion of genes detected as DS")
# ggsave("prop_genes_ds.pdf", path = export_dir,
#        width = 5, height = 5/1.375)









# ~~~~~~~~~~~~~~~~~~~~~~~~ ----


# Long exons with DS ----

# Consider only genes witrh DS
genes_ds <- dpsi |>
  filter(neurA %in% all_neurs_sequenced,
         neurB %in% all_neurs_sequenced) |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique()


# of those, how many have an exon with > 500 bp

exons <- wb_load_exon_coords(281)

big_ex_ds <- exons |>
  filter(gene_id %in% genes_ds) |>
  mutate(exon_length = end-start+1) |>
  filter(exon_length > 1000) |>
  mutate(gene_name = i2s(gene_id, gids))




#~ old Heatmap normd by nb_genes_coexpressedping genes ----


hm_mat <- nb_signif_genes_by_test |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(nb_signif_genes_by_test) |>
  pivot_wider(names_from = neurB,
              values_from = nb_DS_genes) |>
  column_to_rownames("neurA") |> 
  as.matrix()

hm_mat <- hm_mat[sort(rownames(hm_mat)), sort(colnames(hm_mat))]

heatmap_annot <- neuron_properties |>
  column_to_rownames("Neuron_type")

hc <- hclust(dist(hm_mat, method = "canberra"), method = "complete")

nb_genes_coexpressed <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)] >0) %*% (cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)]>0)
all.equal(colnames(hm_mat), colnames(nb_genes_coexpressed))
all.equal(rownames(hm_mat), rownames(nb_genes_coexpressed))



# To set the order (roughly)
weights <- set_names(c("OLL","IL1","VD","DD","PHA","AWC","AVH","DA","PVC","I5","AVM","VC","AIN","OLQ",
                       "RIM","RMD","ASER","ASI","AFD","BAG","ASEL","ASG","AWA","NSM","CAN","ADL",
                       "ASK","IL2","RIA","RIC","AVK","SMD","AVA","AVE","AVG","AWB","PVD","PVM","RIS",
                       "VB","AIY"),
                     exp(1:41)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}

hc2 <- hm_callback(hc)

pheatmap::pheatmap((hm_mat/nb_genes_coexpressed),
                   scale = "none",
                   cluster_rows = hc2,
                   cluster_cols = hc2,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   main = "Normalized by number of nb_genes_coexpressedping genes",
                   # filename = file.path(export_dir, "heatmap_ds_bynb_genes_coexpressed.pdf"),
                   # width = 10,
                   # height = 9
)

table(cutree(hc2, k=2))




#~ other Heatmaps ----



pheatmap::pheatmap(hm_mat, scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   # filename = file.path(export_dir, "211111_heatmap_ds.pdf"),
                   # width = 10,
                   # height = 9
                   )






pheatmap::pheatmap((hm_mat/nb_genes_coexpressed),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   main = "Normalized by number of nb_genes_coexpressedping genes",
                   # filename = file.path(export_dir, "heatmap_ds_bynb_genes_coexpressed_no_annot.pdf"),
                   # width = 10,
                   # height = 9
)


#~ Heatmap normalized by neurons ----

hc <- hclust(dist(hm_mat, method = "canberra"), method = "complete")
hm_mat_by_row <- hm_mat/colSums(hm_mat, na.rm = TRUE)
hm_mat_by_col <- t(hm_mat/colSums(hm_mat, na.rm = TRUE))

pheatmap::pheatmap(hm_mat_by_row, scale = "none",
                   cluster_cols = hc,
                   cluster_rows = hc,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by row",
                   # filename = file.path(export_dir, "211111_heatmap_ds_byrow.pdf"),
                   # width = 10,
                   # height = 9
)
# 
# pheatmap::pheatmap(hm_mat_by_col, scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(t(scale(t(hm_mat), center = FALSE)), scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(hm_mat, scale = "row",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)


#~ Normalization comparison ----

#~~ bargraph ----

tot_genes_per_neur <- nb_signif_genes_by_test |>
  pivot_longer(cols = c("neurA","neurB")) |>
  group_by(value) |>
  summarize(tot_neur = sum(nb_DS_genes),
            .groups = "drop")

nb_signif_genes_by_test |>
  ungroup() |>
  left_join(tot_genes_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_genes_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_DS_genes/tot_neurA,
         nb_by_B = nb_DS_genes/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  pivot_longer(-pair, names_to = "norm_by", values_to = "norm_nb_DS_genes") |>
  ggplot() +
  geom_col(aes(x = pair, y = norm_nb_DS_genes, fill = norm_by))

#~~ scatterplot ----

nb_signif_genes_by_test |>
  ungroup() |>
  left_join(tot_genes_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_genes_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_DS_genes/tot_neurA,
         nb_by_B = nb_DS_genes/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  ggplot(aes(x = nb_by_A, y = nb_by_B)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Number of DS genes for neuron pair normalized by first neuron") +
  ylab("Number of DS genes for neuron pair normalized by second neuron") +
  ggrepel::geom_text_repel(aes(label = pair))

# ggsave("nb_ds_genes_norm_by_both_neurs.pdf", path = export_dir)
# ggsave("nb_ds_genes_norm_by_both_neurs_labels.pdf", path = export_dir)

# How many tests this DS gene appears in?
nb_tests_per_ds_gene <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n())

hist(nb_tests_per_ds_gene$nb_tests, breaks = 100)

nb_tests_per_ds_gene |>
    filter(nb_tests > 600) |>
  mutate(gene_name = i2s(gene_id, gids))

dpsi |>
  filter(p20 >.50 & p05 < .05,
         neurA == "AIY" | neurB == "AIY") |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n()) |>
  pull(nb_tests) |>
  hist(breaks = 100)


xx <- dpsi |>
  filter(neurA == "AIY" | neurB == "AIY") |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  mutate(neur = if_else(neurA == "AIY", neurB, neurA)) |>
  select(-neurA, -neurB) |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = mean(nb_ds > 0))

hist(xx$prop_tests_ds, breaks = 50)



xx <- dpsi |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = sum(nb_ds > 0)/nrow(select(dpsi, neurA, neurB) |> distinct()))

hist(xx$prop_tests_ds, breaks = 50)

xx |>
  filter(prop_tests_ds > .8) |>
  mutate(gene_name = i2s(gene_id, gids))


# Nb of LSVs with DS ----


dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(lsv_id) |> unique() |>
  length()
#> 4611 / 10,232

#~ By neur classes ----
signif_lsvs <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(lsv_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_ds_lsvs = n(),
            .groups = "drop")

signif_lsvs |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_lsvs) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = neurA, y=neurB, fill = nb_ds_lsvs))



#~ Heatmap ----
hm_mat_lsvs <- signif_lsvs |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_lsvs) |>
  pivot_wider(names_from = neurB,
              values_from = nb_ds_lsvs) |>
  column_to_rownames("neurA") |> 
  as.matrix()

hm_mat_lsvs <- hm_mat_lsvs[sort(rownames(hm_mat_lsvs)), sort(colnames(hm_mat_lsvs))]


pheatmap::pheatmap(hm_mat_lsvs, scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   # filename = file.path(export_dir, "211111_heatmap_lsvs.pdf"),
                   # width = 10,
                   # height = 9
)



#~ Heatmap normd by nb_genes_coexpressedping genes ----


nb_genes_coexpressed <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]) %*% cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]
all.equal(colnames(hm_mat_lsvs), colnames(nb_genes_coexpressed))
all.equal(rownames(hm_mat_lsvs), rownames(nb_genes_coexpressed))



# To set the order (roughly)
weights <- set_names(c("OLL","OLQ","PHA","AWC","IL1","M4","VD","DD","AVH","DA","PVC","I5","AVM","VC","AIN","RIM","RMD","ASER","ASI","AFD","BAG","ADF","ASK","IL2","RIA","RIC","AVK","SMD","ASEL","ASG","AWA","NSM","CAN","ADL","AVA","AVE","AVG","AWB","PVD","PVM","RIS","VB","AIY"),
                     exp(1:43)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}


pheatmap::pheatmap((hm_mat_lsvs/nb_genes_coexpressed),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by number of nb_genes_coexpressedping LSVs",
                   # filename = file.path(export_dir, "heatmap_lsvs_bynb_genes_coexpressed.pdf"),
                   # width = 10,
                   # height = 9
)




#~ Heatmap normalized by neurons ----

hc <- hclust(dist(hm_mat_lsvs, method = "canberra"), method = "complete")
hm_mat_by_row <- hm_mat_lsvs/colSums(hm_mat_lsvs, na.rm = TRUE)
hm_mat_by_col <- t(hm_mat_lsvs/colSums(hm_mat_lsvs, na.rm = TRUE))

pheatmap::pheatmap(hm_mat_by_row, scale = "none",
                   cluster_cols = hc,
                   cluster_rows = hc,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by row",
                   # filename = file.path(export_dir, "211111_heatmap_lsvs_byrow.pdf"),
                   # width = 10,
                   # height = 9
)
# 
# pheatmap::pheatmap(hm_mat_by_col, scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(t(scale(t(hm_mat), center = FALSE)), scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(hm_mat, scale = "row",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)



#~ Normalization comparison ----
#~~ bargraph ----

tot_lsvs_per_neur <- signif_lsvs |>
  pivot_longer(cols = c("neurA","neurB")) |>
  group_by(value) |>
  summarize(tot_neur = sum(nb_ds_lsvs),
            .groups = "drop")

signif_lsvs |>
  ungroup() |>
  left_join(tot_lsvs_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_lsvs_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_ds_lsvs/tot_neurA,
         nb_by_B = nb_ds_lsvs/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  pivot_longer(-pair, names_to = "norm_by", values_to = "norm_nb_DS_lsvs") |>
  ggplot() +
  geom_col(aes(x = pair, y = norm_nb_DS_lsvs, fill = norm_by))


#~~ scatterplot ----

signif_lsvs |>
  ungroup() |>
  left_join(tot_lsvs_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_lsvs_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_ds_lsvs/tot_neurA,
         nb_by_B = nb_ds_lsvs/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  ggplot(aes(x = nb_by_A, y = nb_by_B)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Number of DS LSVs for neuron pair normalized by first neuron") +
  ylab("Number of DS LSVs for neuron pair normalized by second neuron") +
  ggrepel::geom_text_repel(aes(label = pair))

# ggsave("nb_ds_lsvs_norm_by_both_neurs.pdf", path = export_dir)
# ggsave("nb_ds_lsvs_norm_by_both_neurs_labels.pdf", path = export_dir)


# How many tests this DS gene appears in?
nb_tests_per_ds_gene <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n())

hist(nb_tests_per_ds_gene$nb_tests, breaks = 100)

nb_tests_per_ds_gene |>
  filter(nb_tests > 600) |>
  mutate(gene_name = i2s(gene_id, gids))

dpsi |>
  filter(p20 >.50 & p05 < .05,
         neurA == "AIY" | neurB == "AIY") |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n()) |>
  pull(nb_tests) |>
  hist(breaks = 100)


xx <- dpsi |>
  filter(neurA == "AIY" | neurB == "AIY") |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  mutate(neur = if_else(neurA == "AIY", neurB, neurA)) |>
  select(-neurA, -neurB) |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = mean(nb_ds > 0))

hist(xx$prop_tests_ds, breaks = 50)



xx <- dpsi |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = sum(nb_ds > 0)/nrow(select(dpsi, neurA, neurB) |> distinct()))

hist(xx$prop_tests_ds, breaks = 50)

xx |>
  filter(prop_tests_ds > .8) |>
  mutate(gene_name = i2s(gene_id, gids))





# ~~~~~~~~~~~  ----














