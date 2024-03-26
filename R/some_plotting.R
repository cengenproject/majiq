
## Initializations ----
# library(tidyverse)
# library(wbData)
suppressPackageStartupMessages(library(ggbio))

gids <- wbData::wb_load_gene_ids("WS277")


scale_color_randomized <- function(...){
  discrete_scale("color",
                 "hue_rand",
                 palette = function(n) sample(scales::hue_pal(...)(n)))
}

txdb_WS277 <- GenomicFeatures::makeTxDbFromGFF(wb_get_gtf_path(277),
                                               dataSource = "Wormbase",
                                               organism = "Caenorhabditis elegans")
all_exons <- GenomicFeatures::exonsBy(txdb_WS277, by = "gene")


## fctions to plot gene and LSV positions ----

gene_as_gr <- function(my_gene){
  check_collisions(my_gene)
  
  pmjq %>%
    filter(gene_id == my_gene) %>%
    dplyr::select(sh_lsv_id, chr, strand, ends_with("coords")) %>%
    distinct() %>%
    separate(jct_coords, ";",
             into = paste0("Cj_",1:50),
             fill="right") %>%
    separate(ir_coords, ";",
             into = paste0("Ci_",51:70),
             fill="right") %>%
    pivot_longer(starts_with("C", ignore.case = FALSE),
                 names_to="junction_id") %>%
    filter(nchar(value) > 0) %>%
    separate(value,
             into = c("start", "end"),
             "-") %>%
    GenomicRanges::GRanges()
}

plot_single_lsv <- function(gr, my_gene, my_lsv){
  gr <- gr[gr$sh_lsv_id == my_lsv,]
  tracks(
    ggplot() +
      geom_rect(unique(all_exons[[my_gene]])) +
      theme_minimal() +
      ggtitle(my_gene),
    {no_ir <- ggplot() +
      geom_arch(gr[startsWith(gr$junction_id, "Cj"),],
                aes(color = sh_lsv_id, height = 10)) +
      ggrepel::geom_text_repel(data = as_tibble(gr[startsWith(gr$junction_id, "Cj"),]),
                               mapping = aes(x=(end+start)/2,y=10,label = junction_id),
                               ylim = c(10,20),nudge_y = 1) +
      expand_limits(y=15) +
      theme_minimal() +
      theme(legend.position="none")
    
    if(length(gr[startsWith(gr$junction_id, "Ci"),])>0){
      no_ir +
        geom_rect(gr[startsWith(gr$junction_id, "Ci"),],
                  aes(fill = sh_lsv_id)) +
        geom_text(data = as_tibble(gr[startsWith(gr$junction_id, "Ci"),]),
                  mapping = aes(x=(end+start)/2,y=1.1,label = sh_lsv_id),
                  size=3)
    } else{
      no_ir
    }
    }
  )
}



plot_gene_LSVs <- function(gr, my_gene){
  tracks(
    ggplot() +
      geom_rect(unique(all_exons[[my_gene]])) +
      theme_minimal() +
      ggtitle(my_gene),
    ggplot() +
      geom_arch(gr[startsWith(gr$junction_id, "Cj"),],
                aes(color = sh_lsv_id, height = 10)) +
      geom_rect(gr[startsWith(gr$junction_id, "Ci"),],
                aes(fill = sh_lsv_id)) +
      ggrepel::geom_text_repel(data = as_tibble(gr[startsWith(gr$junction_id, "Cj"),]),
                               mapping = aes(x=(end+start)/2,y=10,label = sh_lsv_id),
                               ylim = c(10,20),nudge_y = 1) +
      {if(any(startsWith(gr$junction_id, "Ci"))){
        geom_text(data = as_tibble(gr[startsWith(gr$junction_id, "Ci"),]),
                  mapping = aes(x=(end+start)/2,y=1.1,label = sh_lsv_id),
                  size=3)
      } else{
        NULL
      }} +
      expand_limits(y=15) +
      theme_minimal() +
      theme(legend.position="none")
  )
}

check_collisions <- function(my_gene){
  stopifnot(all.equal(n_distinct(pmjq[pmjq$gene_id == my_gene, "lsv_id"]),
                      n_distinct(pmjq[pmjq$gene_id == my_gene, "sh_lsv_id"])))
}





## Prep data ----

pmjq <- mjq |>
  mutate(chr = AnnotationDbi::select(txdb_WS277, gene_id,"CDSCHROM","GENEID")$CDSCHROM,
         strand = AnnotationDbi::select(txdb_WS277, gene_id,"CDSSTRAND","GENEID")$CDSSTRAND)

# create short LSV ID (as visual guide, rare collisions OK)
set.seed(654)
xx <- data.frame(row.names = unique(pmjq$lsv_id),
                 short_id = map_chr(1:length(unique(pmjq$lsv_id)),
                                    ~paste0(sample(LETTERS, 3, replace = TRUE),collapse="")))
pmjq$sh_lsv_id <- xx[pmjq$lsv_id,]






# Plot with ggbio ----


# plot individual genes ----

my_gene <- "WBGene00012389"

my_gene <- s2i("casy-1", gids)
my_gene <- s2i("ric-7", gids)
my_gene <- s2i("mbk-2", gids)
my_gene <- s2i("skn-1", gids)
my_gene <- s2i("unc-40", gids)
my_gene <- s2i("unc-64", gids)
my_gene <- s2i("sax-3", gids)
my_gene <- s2i("sca-1", gids)

# plotSpliceGraph(sg_feat, geneName = my_gene)


# Plot all LSV quantifs for 1 gene
gr <- gene_as_gr(my_gene)

plot_gene_LSVs(gr, my_gene)


pmjq %>%
  filter(gene_id == my_gene) %>% #, neuron %in% c("M4","BAG")
  ggplot(aes(x=neuron, y =psi, ymin=psi-sd,ymax=psi+sd, color=neuron)) +
  geom_point() +
  geom_errorbar() +
  facet_grid(vars(junction_id), vars(sh_lsv_id)) +
  theme_bw() +
  scale_color_randomized() +
  ggtitle(i2s(my_gene, gids)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))



# Plot sigle LSV quantif
my_lsv <- "EGH"

plot_single_lsv(gr, my_gene, my_lsv)


pmjq %>%
  filter(gene_id == my_gene & sh_lsv_id == my_lsv) %>%   # , neuron %in% c("AWA","AWB","ASG","AVG")
  ggplot(aes(x=neuron, y =psi, ymin=psi-sd,ymax=psi+sd, color=neuron)) +
  geom_point() +
  geom_errorbar() +
  facet_wrap(~junction_id) +
  scale_color_randomized() +
  ggtitle(paste0(i2s(my_gene, gids), ", ", my_lsv)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# As violin (same as VOILA)
pmjq %>%
  filter(gene_id == my_gene,
         sh_lsv_id == my_lsv,
         # neuron %in% c("CAN","ASEL","ASG","AVG"),
         junction_id == 2) %>%
  mutate(alpha = psi*(psi*(1-psi)/sd^2 -1),
         beta = (1-psi)*(psi*(1-psi)/sd^2 -1),
         points = map2(alpha, beta, ~rbeta(1000, .x,.y))) %>%
  unnest(points) %>%
  ggplot(aes(x=neuron)) +
  geom_violin(aes(y=points, color = neuron), width=3) +
  geom_boxplot(aes(y=points), width = .1, fill="black", outlier.alpha = 0) +
  geom_point(aes(y = psi), size = .8,color="white") +
  # geom_errorbar(aes(ymin=psi-sd,ymax=psi+sd), width = .1) +
  ylim(c(0,1)) +
  theme_classic() +
  ylab("PSI")

# As violin (same as VOILA)
pmjq %>%
  filter(gene_id == my_gene,
         sh_lsv_id == my_lsv,
         # neuron %in% c("CAN","ASEL","ASG","AVG"),
         junction_id %in% c(1,2,3)) %>%
  mutate(alpha = psi*(psi*(1-psi)/sd^2 -1),
         beta = (1-psi)*(psi*(1-psi)/sd^2 -1),
         points = map2(alpha, beta, ~rbeta(1000, .x,.y))) %>%
  unnest(points) %>%
  ggplot(aes(x=neuron)) +
  geom_violin(aes(y=points, color = neuron), width=4) +
  geom_boxplot(aes(y=points), width = .2, fill="black", outlier.alpha = 0) +
  geom_point(aes(y = psi), size = .8,color="white") +
  facet_grid(rows = vars(junction_id)) +
  # geom_errorbar(aes(ymin=psi-sd,ymax=psi+sd), width = .1) +
  ylim(c(0,1)) +
  theme_bw() +
  ylab("PSI") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Check LSV in original data
(mjq %>%
    filter(gene_id == my_gene & sh_lsv_id == my_lsv) %>% pull(lsv_id))[1]

mjqdta %>%
  filter(`Gene ID` == my_gene,
         `LSV ID` == "WBGene00000912:s:10751038-10751188") %>%
  dplyr::slice(1) %>%
  as.data.frame()





## ~~~~~ END ~~~~~ ----
# ~~~~~~~~~~~~~~~ ----


# path_psi <- "data/210215_majiq/psi/"
# files_psi <- list.files(path_psi, pattern = "*.tsv")
# 
# mjqdta <- map_dfr(files_psi,
#                   ~read_tsv(file.path(path_psi, .x),
#                             col_types = cols(
#                               `Gene ID` = col_character(),
#                               `LSV ID` = col_character(),
#                               `LSV Type` = col_character(),
#                               `E(PSI) per LSV junction` = col_character(),
#                               `StDev(E(PSI)) per LSV junction` = col_character(),
#                               A5SS = col_logical(),
#                               A3SS = col_logical(),
#                               ES = col_logical(),
#                               `Num. Junctions` = col_double(),
#                               `Num. Exons` = col_double(),
#                               `Junctions coords` = col_character(),
#                               `IR coords` = col_character()
#                             ),
#                             na = "na") %>%
#                     add_column(neuron = str_split(.x,"\\.")[[1]][1]))
# 
# # separate data in E(PSI) and Std(PSI) fields
# mjq <- mjqdta %>%
#   separate(`E(PSI) per LSV junction`,
#            into = paste0("e_junction_",1:100),
#            sep = ";", fill = "right") %>%
#   separate(`StDev(E(PSI)) per LSV junction`,
#            into = paste0("s_junction_",1:100),
#            sep = ";", fill = "right") %>%
#   pivot_longer(contains("_junction_"),
#                names_to = c(".value", "junction"),
#                names_pattern = "(e|s)_junction_(\\d+)") %>%
#   filter(!is.na(e)) %>%
#   dplyr::select(neuron,
#                 gene_id = `Gene ID`,
#                 lsv_id = `LSV ID`,
#                 junction_id = junction,
#                 psi = e,
#                 sd = s,
#                 jct_coords = `Junctions coords`,
#                 ir_coords = `IR coords`) %>%
#   mutate(psi = as.double(psi),
#          sd = as.double(sd),
#          chr = AnnotationDbi::select(txdb_WS277, gene_id,"CDSCHROM","GENEID")$CDSCHROM,
#          strand = AnnotationDbi::select(txdb_WS277, gene_id,"CDSSTRAND","GENEID")$CDSSTRAND)









# Plot with SGSeq ----

# Prepare global
library(SGSeq)

sg_feat <- wb_get_gtf_path(277) %>%
  importTranscripts() %>%
  convertToTxFeatures() %>%
  convertToSGFeatures()


# Plot gene

my_gene <- "WBGene00012389"

my_gene <- s2i("casy-1", gids)
my_gene <- s2i("ric-7", gids)
my_gene <- s2i("mbk-2", gids)

plotSpliceGraph(sg_feat, geneName = my_gene, main = i2s(my_gene, gids))
plotSpliceGraph(sg_feat, geneName = my_gene, toscale = "gene")


# Prepare LSVs

gr <- gene_as_gr(my_gene)

unique(all_exons[[my_gene]])


gr0 <- gr[startsWith(gr$junction_id, "Cj"),]

gr0$sh_j_id <- paste(gr0$sh_lsv_id, gr0$junction_id, sep=":")

gr[startsWith(gr$junction_id, "Ci"),]
plotSpliceGraph(sg_feat, geneName = my_gene, main = i2s(my_gene, gids),
                ranges = split(gr0, gr0$sh_j_id), ranges_ypos = c(0.1,0.7))


aa(sg_feat)
plotSpliceGraph(sg_feat, geneName = my_gene)

ggrepel::geom_text_repel(data = as_tibble(gr[startsWith(gr$junction_id, "Cj"),]),
                         mapping = aes(x=(end+start)/2,y=10,label = sh_lsv_id),
                         ylim = c(10,20),nudge_y = 1) +
  geom_text(data = as_tibble(gr[startsWith(gr$junction_id, "Ci"),]),
            mapping = aes(x=(end+start)/2,y=1.1,label = sh_lsv_id),
            size=3)



gr0 <- gr[gr$sh_lsv_id == my_lsv,]
ggplot() +
  geom_rect(unique(all_exons[[my_gene]])) +
  theme_minimal() +
  ggtitle(my_gene)


no_ir <- ggplot() +
  geom_arch(gr[startsWith(gr$junction_id, "Cj"),],
            aes(color = sh_lsv_id, height = 10)) +
  ggrepel::geom_text_repel(data = as_tibble(gr[startsWith(gr$junction_id, "Cj"),]),
                           mapping = aes(x=(end+start)/2,y=10,label = junction_id),
                           ylim = c(10,20),nudge_y = 1) +
  expand_limits(y=15) +
  theme_minimal() +
  theme(legend.position="none")








# END MAIN SCRIPT ----


##

# full isoforms structures
p_lsv <- ggplot() +
  geom_arch(gr[startsWith(gr$name, "Cj"),],
            aes(color = sh_lsv_id)) +
  geom_rect(gr[startsWith(gr$name, "Ci"),],
            aes(fill = sh_lsv_id))

tracks(transcripts = p_genomic,
       LSVs = p_lsv,
       heights = c(4,1))



#
library(SplicingGraphs)
options(ucscChromosomeNames=FALSE)
# sg <- SplicingGraphs(txdb_WS277)
# saveRDS(sg, "data/SplicingGraphs_WS277.rds")
sg <- readRDS("data/SplicingGraphs_WS277.rds")

plot(sg[[my_gene]])

# ----

# Various previous tries
mjqdta %>%
  filter(`Gene ID` == "WBGene00000092") %>%
  # filter(`LSV ID` == "WBGene00000092:s:5303474-5303600") %>%
  View()

autoplot(txdb_WS277, which = gene_tr)


grlist <- mjq %>%
  filter(gene_id == "WBGene00000092") %>%
  dplyr::select(lsv_id, chr, ends_with("coords")) %>%
  distinct() %>%
  separate(jct_coords, ";",
           into = paste0("Cj_",1:50),
           fill="right") %>%
  separate(ir_coords, ";",
           into = paste0("Ci_",51:70),
           fill="right") %>%
  pivot_longer(starts_with("C", ignore.case = FALSE)) %>%
  filter(nchar(value) > 0) %>%
  separate(value,
           into = c("start", "end"),
           "-") %>%
  group_by(lsv_id) %>%
  group_split() %>%
  IRanges::DataFrameList() %>%
  as("GRangesList")

grlist
autoplot(grlist, geom="link")

ggplot() +
  geom_alignment(grlist, type="model")

GenomicRanges::GRangesList(.$value)

aa <- GenomicRanges::GRangesList(yy$value, lsv_id = yy$lsv_id)


tracks(transcripts = autoplot(txdb_WS277, which = gene_tr),
       LSVs = autoplot(grlist, geom="splice"))

# map(xx$jct_coords, ~ strsplit(.x, ";"))
# yy <- xx$jct_coords[1]
GenomicRanges::GRanges(strsplit(yy, ";")[[1]])
GenomicRanges::GRanges()  

# Plot LSVs positions for 1 gene
gene_tr <- GenomicFeatures::transcriptsBy(txdb_WS277, by="gene", use.names=FALSE)[[my_gene]]


ggplot() +
  geom_rect(unique(all_exons[[my_gene]])) +
  geom_arch(gr[startsWith(gr$name, "Cj"),],
            aes(color = sh_lsv_id)) +
  geom_rect(gr[startsWith(gr$name, "Ci"),],
            aes(fill = sh_lsv_id))

# In genomic coordinates with all full transcripts
p_genomic <- autoplot(txdb_WS277, which = gene_tr)
p_lsv <- ggplot() +
  geom_arch(gr[startsWith(gr$name, "Cj"),],
            aes(color = sh_lsv_id)) +
  geom_rect(gr[startsWith(gr$name, "Ci"),],
            aes(fill = sh_lsv_id))

tracks(transcripts = p_genomic,
       LSVs = p_lsv,
       heights = c(2,1))


# heatmap ----


xx <- mjqdta %>%
  dplyr::select(id = `LSV ID`, neuron, PSI_all = `E(PSI) per LSV junction`) %>%
  mutate(PSI = map_dbl(PSI_all, ~as.double(strsplit(.x, ";")[[1]][1])),
         neuron = map_chr(neuron, ~strsplit(.x, "\\.")[[1]][1])) %>%
  dplyr::select(-PSI_all) %>%
  group_by(id, neuron) %>%
  summarize(PSI = mean(PSI)) %>%
  pivot_wider(names_from = "neuron", values_from = "PSI")

yy <- as.matrix(xx[,-1])
rownames(yy) <- xx$id

pheatmap::pheatmap(yy,cluster_rows = FALSE, cluster_cols = FALSE)


dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
gplots::heatmap.2(yy, na.color="grey", distfun=dist_no_na,scale="col",trace="none")
