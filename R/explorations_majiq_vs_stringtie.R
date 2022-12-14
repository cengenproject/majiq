library(tidyverse)
library(wbData)
gids <- wb_load_gene_ids(277)

# from majiq_global_splicing.R
quantif_genes <- unique(mjqdta$`Gene ID`)

# from explore_build_database.R
built_genes <- unique(xx$gene_id)

plot(eulerr::euler(list(quantif = quantif_genes, built = built_genes)))

# genes with biological role in neuron and likely alt spliced by visual examination

my_candidates <- readxl::read_excel("../candid_based_on_bio.xlsx") |>
  filter(interest > 1) |>
  pull(gene)

my_candidates

plot(eulerr::euler(list(quantif = quantif_genes, built = built_genes, candidates = s2i(my_candidates, gids))))



not_candidates <- readxl::read_excel("../candid_based_on_bio.xlsx") |>
  filter(interest == 0) |>
  pull(gene)
not_candidates
plot(eulerr::euler(list(quantif = quantif_genes,
                        built = built_genes,
                        candidates = s2i(my_candidates, gids),
                        not_candidates = s2i(not_candidates, gids))))


UpSetR::upset(UpSetR::fromList(list(quantif = quantif_genes,
                                    built = built_genes,
                                    candidates = s2i(my_candidates, gids),
                                    not_candidates = s2i(not_candidates, gids))))

# Screen some genes manually
xx <- sample(built_genes[! built_genes %in% quantif_genes], 1)

mjqdta |> 
  filter(`Gene ID` == xx)

i2s(xx, gids); writeClipboard(i2s(xx, gids))

xx; writeClipboard(xx)







# Non neuronal genes ----
known_nonneur <- c("ceh-17","pal-1","ceh-22","ceh-51","ceh-5","ceh-49","eyg-1","pax-3","ceh-33","ceh-40","ceh-75","ceh-83","ceh-87","ceh-92","myo-1","myo-2","myo-3","hlh-1","unc-22","lev-11","fkh-6","fhod-1","lon-3","bli-1","ifb-2","act-5")
known_panneur <- c("nsf-1","rab-3","sng-1","snt-1","unc-10","unc-18","ehs-1","unc-11","unc-57","snn-1","unc-104","egl-3","egl-21","maco-1","rgef-1","ceh-30","ceh-21","ceh-38","ceh-39","ceh-41","ceh-48","ceh-74","ceh-82","ceh-88","ceh-44","ceh-58","ztf-11","ztf-18","ztf-26","ztf-3","ztf-4","ztf-7","ztf-9","zip-2","zip-4","ric-4","snb-1","unc-64","syd-2","ric-19","unc-31","unc-108","tbb-1")
known_broadneur <- c("lim-4","unc-4","ceh-43","mab-5","egl-5","ceh-54","ceh-32","ceh-79","mgl-3","lin-39","cog-1","ceh-23","ceh-34","irx-1","ceh-14","lim-7","unc-42","alr-1","lin-11","vab-3","inx-19","unc-86","gar-3","che-7","inx-18","ceh-20","unc-62","mgl-1","gar-1","gar-2","nsy-7","gbb-2","mgl-2","ceh-18","zfh-2","zag-1","inx-2","inx-7","inx-14","unc-7","gbb-1","inx-1","inx-10","unc-9")
any_truth <- read_csv("../../../cengen_10x/ROC/output/ground_truth.csv")
plot(eulerr::euler(list(quantif = quantif_genes,
                        nonneur = any_truth$`gene ID`[rowSums(any_truth[,4:131]) == 0],
                        panneur = any_truth$`gene ID`[rowSums(any_truth[,4:131]) >= 127],
                        broad_neur = any_truth$`gene ID`[rowSums(any_truth[,4:131]) > 10 &
                                                           rowSums(any_truth[,4:131]) < 127])))

i2s(intersect(quantif_genes, any_truth$`gene ID`[rowSums(any_truth[,4:131]) == 0]), gids)


