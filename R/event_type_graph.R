# Starting script to trace the splicing graph to determine whether LSV affects TSS or polyA

# Check build DB
library(DBI)


# Retrive splice graph
con <- dbConnect(RSQLite::SQLite(), "data/splicegraph.sql")

# junctions
res <- dbSendQuery(con, "SELECT * FROM junction
                              WHERE has_reads;")
sg_sj <- dbFetch(res)
dbClearResult(res)


head(sg_sj)


dbListTables(con)

# Junction ----
dbListFields(con, "junction")


dbClearResult(res)


head(xx)
nrow(xx)

length(unique(xx$gene_id))



# Junction reads ----
dbListFields(con, "junction_reads")


res <- dbSendQuery(con, "SELECT * FROM junction_reads
                              WHERE reads > 3;")
xx <- dbFetch(res)
dbClearResult(res)


head(xx)
nrow(xx)

length(unique(xx$junction_gene_id))

xx <- xx |>
  mutate(sj_coords = paste0(junction_gene_id,":",junction_start,"-",junction_end))

sg_sj <- sg_sj |>
  mutate(sj_coords = paste0(gene_id,":",start,"-",end))

jreads_uniq <- xx |>
  select(sj_coords) |>
  distinct()

j_uniq <- sg_sj |>
  select(sj_coords) |>
  distinct()

nrow(jreads_uniq)
nrow(j_uniq)


# Alt start ----
dbListFields(con, "alt_start")


res <- dbSendQuery(con, "SELECT * FROM alt_start;")
xx <- dbFetch(res)
dbClearResult(res)


head(xx)
nrow(xx)
table(duplicated(xx$gene_id))

length(unique(xx$junction_gene_id))





# dbDisconnect(con)