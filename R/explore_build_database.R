# Check build DB
library(DBI)

con <- dbConnect(RSQLite::SQLite(), "data/splicegraph.sql")


dbListTables(con)

# Junction ----
dbListFields(con, "junction")

res <- dbSendQuery(con, "SELECT gene_id FROM junction
                              WHERE has_reads AND NOT is_constitutive;")
xx <- dbFetch(res)
dbClearResult(res)


head(xx)
nrow(xx)

length(unique(xx$gene_id))



# Junction reads ----
dbListFields(con, "junction_reads")


res <- dbSendQuery(con, "SELECT junction_gene_id FROM junction_reads
                              WHERE reads > 3;")
xx <- dbFetch(res)
dbClearResult(res)


head(xx)
nrow(xx)

length(unique(xx$junction_gene_id))



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