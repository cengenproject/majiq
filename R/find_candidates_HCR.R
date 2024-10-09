# Load deltaPSI to find events with variation, and filter to keep only genes
# that are not too broadly expressed.


# Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)



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

# Identify broadly expressed genes
expression_breadth <- rowMeans(gene_expr_int)
hist(expression_breadth, breaks = 50)

table(expression_breadth < .2)

restricted_pattern <- names(expression_breadth)[expression_breadth < .2]

# differentially spliced
ds_genes <- dpsi |>
  filter(p20 >.50 & p05 < .1 & abs(dpsi) > 0.6 ) |>
  pull(gene_id) |> unique()



# Keep only DS event if it happens btw neurons where gene is expressed
dpsi_expr <- dpsi |>
  filter(p20 >.50 & p05 < .1 & abs(dpsi) > 0.6 ) |>
  mutate(exprA = map2_lgl(gene_id, neurA, \(.g, .n) gene_expr_int[.g, .n]),
         exprB = map2_lgl(gene_id, neurB, \(.g, .n) gene_expr_int[.g, .n]))

ds_in_expr <- dpsi_expr |>
  filter(exprA & exprB) |>
  pull(gene_id) |> unique()




neur_genes <- wormDatasets::genes_by_pattern


table(ds_genes %in% neur_genes$nonneuronal)


# Selections ----

list(`Differentially Spliced` = ds_genes,
     `DS in neurons where expressed` = ds_in_expr) |>
  eulerr::euler() |>
  plot(quantities = TRUE)




list(`Restricted pattern` = restricted_pattern,
     `DS in neurons where expressed` = ds_in_expr,
     `Non-neuronal expression` = neur_genes$nonneuronal) |>
  eulerr::euler() |>
  plot(quantities = TRUE)


# Get list
candidates <- ds_in_expr |>
  setdiff(neur_genes$nonneuronal) |>
  intersect(restricted_pattern)

length(candidates)
head(candidates)
head(i2s(candidates, gids))
writeClipboard(i2s(candidates, gids))
writeClipboard(candidates)



# To Excel
xx <- s2i(clipr::read_clip(), gids)
writeClipboard(xx)

xx <- clipr::read_clip()

length(xx)
table(candidates %in% xx)
table(xx %in% candidates)

clipr::write_clip(xx %in% candidates)


# previous list
candidates2 <- dpsi |>
  filter(p20 >.50 & p05 < .05 & abs(dpsi) > 0.6 ) |>
  pull(gene_id) |> unique() |>
  setdiff(neur_genes$nonneuronal) |>
  intersect(restricted_pattern)


list(new = candidates,
     old = candidates2) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

list(excel = xx,
     old_candidates = candidates2) |>
  eulerr::euler() |>
  plot(quantities = TRUE)

