# Load quantifications from MAJIQ and use them for general descriptions


## Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(277)

# Load and format data ----
path_psi <- "data/psi/"
files_psi <- list.files(path_psi, pattern = "*.tsv")

# Load file, one row per LSV
mjqdta <- map_dfr(files_psi,
                  ~read_tsv(file.path(path_psi, .x),
                            col_types = cols(
                              `Gene ID` = col_character(),
                              `LSV ID` = col_character(),
                              `LSV Type` = col_character(),
                              `E(PSI) per LSV junction` = col_character(),
                              `StDev(E(PSI)) per LSV junction` = col_character(),
                              A5SS = col_logical(),
                              A3SS = col_logical(),
                              ES = col_logical(),
                              `Num. Junctions` = col_double(),
                              `Num. Exons` = col_double(),
                              `Junctions coords` = col_character(),
                              `IR coords` = col_character()
                            ),
                            na = "na") |>
                    add_column(neuron = str_split(.x,"\\.")[[1]][1]))

# separate data in E(PSI) and Std(PSI) fields
# We get one row per junction
mjq <- mjqdta |>
  separate_rows(`E(PSI) per LSV junction`,
                `StDev(E(PSI)) per LSV junction`,
                sep = ";") |>
  group_by(neuron, `LSV ID`) |>
  mutate(junction_id = row_number()) |>
  ungroup() |>
  dplyr::select(neuron,
                gene_id = `Gene ID`,
                lsv_id = `LSV ID`,
                junction_id,
                psi = `E(PSI) per LSV junction`,
                sd = `StDev(E(PSI)) per LSV junction`,
                jct_coords = `Junctions coords`,
                ir_coords = `IR coords`) |>
  mutate(psi = as.double(psi),
         sd = as.double(sd),
         junction_id = factor(junction_id))





# Global LSV characterizations ----

nb_jctions_per_lsv <- mjq |>
  group_by(gene_id, neuron, lsv_id) |>
  summarize(nb_jctions = n()) |> 
  group_by(gene_id, lsv_id) |>
  summarize(max_nb_jctions = max(nb_jctions))

ggplot(nb_jctions_per_lsv) +
  theme_classic() +
  geom_histogram(aes(x=max_nb_jctions),bins=60, color = "white") +
  scale_x_continuous(trans = "log1p", breaks = c(0,1:10,12,15,20,30,40,50,100)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,10,100,1000,10000)) +
  xlab("Number of junctions per LSV")

# A look at the LSVs with a lot of junctions
nb_jctions_per_lsv |>
  filter(max_nb_jctions >= 10) |>
  group_by(gene_id) |>
  summarize(nb_jctions = paste0(max_nb_jctions, collapse=", "),
            max_max = max(max_nb_jctions)) |>
  arrange(desc(max_max)) |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  select(gene_name,
         `Number of junctions per LSV` = nb_jctions) |>
  clipr::write_clip()



nb_lsv_per_gene <- mjq |>
  group_by(neuron, gene_id)|>
  summarize(nb_lsv = n()) |>
  group_by(gene_id) |>
  summarize(max_nb_lsv = max(nb_lsv))

ggplot(nb_lsv_per_gene) +
  theme_classic() +
  geom_histogram(aes(x=max_nb_lsv),bins=60, color = "white") +
  scale_x_continuous(trans = "log1p", breaks = c(0,1:10,12,15,20,30,40,50,100,200)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,10,100,1000,2000)) +
  xlab("Number of LSV per gene")

# Look at genes with many LSVs
nb_lsv_per_gene |>
  filter(max_nb_lsv >= 30) |>
  arrange(desc(max_nb_lsv)) |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  clipr::write_clip()


## Call differential splicing ----


# Examples to develop method

mjq |>
  filter(gene_id == s2i("unc-40", gids),
         lsv_id == "WBGene00006776:s:5682284-5682373") |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mjq |>
  filter(gene_id == s2i("unc-40", gids)) |>
  group_by(lsv_id, junction_id) |>
  mutate(is_DS = max(psi-sd) > min(psi+sd) |  min(psi+sd) < max(psi-sd)) |> 
  pull(is_DS)


mjq |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549") |>
  group_by(lsv_id, junction_id) |>
  mutate(is_DS = max(psi-sd) > .25 & min(psi+sd) < .25 |
                 min(psi+sd) < .75 & max(psi-sd) > .75) |> 
  select(1,3,4,5,6,9)


mjq |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549") |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




mjq |>
  filter(gene_id == s2i("unc-40", gids)) |> 
  pull(lsv_id) |> unique()



all_ds |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549")



# Compute DS ----

all_ds <- mjq |>
  group_by(lsv_id, junction_id) |>
  mutate(is_DS = max(psi-sd) > .25 & min(psi+sd) < .25 |
           min(psi+sd) < .75 & max(psi-sd) > .75) |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_DS = any(is_DS),
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_DS = sum(is_DS),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(has_DS = any(nb_DS > 0),
            .groups = "drop")


plot(all_ds$nb_jnctions,all_ds$nb_DS)

all_ds |>
  group_by(lsv_id) |>
  summarize(has_DS = nb_DS > 0) |>
  pull(has_DS) |>
  table()
