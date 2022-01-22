library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

export_dir <- "presentations/2021-11_for_preprint_december/"


# Load deltaPSI
dpsi <- readRDS("intermediates/211130_dpsi.rds")

# Load modulizer
modul_summ <- read_tsv("data/2021-11-30_outs/voila_modulizer/summary.tsv",comment = "#",
               col_types = cols(
                 module_id = col_character(),
                 gene_id = col_character(),
                 gene_name = col_character(),
                 seqid = col_character(),
                 strand = col_character(),
                 lsv_id = col_character(),
                 cassette = col_double(),
                 tandem_cassette = col_double(),
                 alt3 = col_double(),
                 alt5 = col_double(),
                 putative_alt3 = col_double(),
                 putative_alt5 = col_double(),
                 alt3_5 = col_double(),
                 mxe = col_double(),
                 ir = col_double(),
                 ale = col_double(),
                 afe = col_double(),
                 putative_ale = col_double(),
                 putative_afe = col_double(),
                 orphan_junc = col_double(),
                 other = col_double(),
                 constitutive_junc = col_double(),
                 persistent_ir = col_double(),
                 multi_exon_spanning = col_double(),
                 complex = col_logical(),
                 denovo_juncs = col_double(),
                 denovo_introns = col_double(),
                 num_events = col_double(),
                 module_event_combination = col_character()
               ))



# Extract lsv with signif dpsi
signif_lsv <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(lsv_id) |> unique()


# map_lgl(signif_lsv,
#        ~ any(str_detect(modul_summ$lsv_id, .x))) |>
#   table()

match_mod_dpsi <- map_int(signif_lsv,
        ~ which(str_detect(modul_summ$lsv_id, .x)))

length(match_mod_dpsi)
head(match_mod_dpsi)

dmod <- modul_summ |>
  slice(match_mod_dpsi) |>
  distinct()

# saveRDS(dmod, "intermediates/220118_dmod.rds")

head(sort(table(dmod$module_event_combination), decreasing = TRUE),10)

event_types <- list(`Alt. 3' ss`      = dmod$module_id[!is.na(dmod$alt3)],
                    `Alt. 5' ss`      = dmod$module_id[!is.na(dmod$alt5)],
                    `Alt. last exon`  = dmod$module_id[!is.na(dmod$ale)],
                    `Alt. first exon` = dmod$module_id[!is.na(dmod$afe)],
                    `Multiple exons`  = dmod$module_id[!is.na(dmod$mxe) | !is.na(dmod$multi_exon_spanning)],
                    Cassette          = dmod$module_id[!is.na(dmod$cassette)],
                    `Intron Retention`= dmod$module_id[!is.na(dmod$ir)])

p <- UpSetR::upset(UpSetR::fromList(event_types),
              nsets = 7,nintersects = 20,
              sets.x.label = "Number of modules",
              text.scale = 1.5)

p

pdf(file.path(export_dir, "upset_module_type.pdf"),
    width = 8, height = 5)
p
dev.off()



# Lump together the "complex" events that don't have a single type.

tot_nb_modules <- nrow(dmod)

dmod2 <- dmod |>
  select(module_id, cassette,ir,mxe,multi_exon_spanning,afe,ale,alt3,alt5,module_event_combination) |>
  mutate(across(-c(module_id, module_event_combination),
                as.character))

for(categ in c("cassette","ir","mxe", "multi_exon_spanning","afe","ale","alt3","alt5")){
  dmod2[!is.na(dmod2[[categ]]), categ] <- categ
}

dmod2[is.na(dmod2)] <- ""
dmod2$any_multi <- if_else(dmod2$mxe != "" | dmod2$multi_exon_spanning != "",
                           "multi",
                           "")
dmod2 <- dmod2 |>
  select(-mxe, -multi_exon_spanning)

dmod2$nb_cat <- apply(dmod2 |>
                        select(-module_id, -module_event_combination),
                      MARGIN = 1,
                      \(x) sum(x != ""))

dmod2$category <- if_else(dmod2$nb_cat > 1,
                          "complex",
                          NA_character_)
dmod2$category[is.na(dmod2$category)] <- apply(dmod2[is.na(dmod2$category),] |>
                                                 select(-module_id, -module_event_combination, -nb_cat, -category),
                                               MARGIN = 1,
                                               \(.x) paste0(.x, collapse = ""))

dmod2 <- dmod2[dmod2$nb_cat != 0,]

dmod2 |>
  mutate(category = recode_factor(category,
                                  alt3 = "Alt. 3' ss",
                                  alt5 = "Alt. 5' ss",
                                  ir   = "Intron retention",
                                  cassette = "Cassette exon",
                                  afe = "Alt. first exon",
                                  ale = "Alt. last exon",
                                  multi = "Multi-exon",
                                  complex = "Complex")) |>
  select(module_id, category) |>
  group_by(category) |>
  summarize(nb_modules = n(),
            .groups = "drop") |>
  mutate(percent = round(100*nb_modules/tot_nb_modules),
         percent = paste0(percent, "%")) |>
  ggplot(aes(x = category, y = nb_modules, label = percent)) +
  theme_classic() +
  geom_col() +
  geom_text(nudge_y = 25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Number of modules") + xlab(NULL)

# ggsave("DS_module_category.pdf", path = export_dir,
#        width = 5, height = 4, units = "in")
