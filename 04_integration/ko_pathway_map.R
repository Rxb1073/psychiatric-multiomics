## Build and cache the KO -> KEGG pathway map.
## Run ONCE, on a login node (needs internet). Everything downstream reads the cache.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(KEGGREST))

OUT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs/06_kegg_axis"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("querying KEGG for ko -> pathway links...\n")
lnk <- keggLink("pathway", "ko")                       # names = ko:K#####, values = path:ko#####
map <- data.frame(KO      = sub("^ko:", "", names(lnk)),
                  pathway = sub("^path:", "", as.character(lnk)),
                  stringsAsFactors = FALSE)
map <- map[grepl("^ko[0-9]{5}$", map$pathway), ]        # drop map#####/ko01100-style duplicates? keep all koNNNNN
map$suffix <- sub("^ko", "", map$pathway)               # the five digits shared with hsa

cat("querying pathway names...\n")
nm <- keggList("pathway", "ko")
names(nm) <- sub("^path:", "", names(nm))
map$pathway_name <- nm[map$pathway]

## KEGG BRITE/global maps are not metabolic pathways in the usual sense and
## dominate any enrichment; flag rather than drop, so the choice stays visible.
GLOBAL <- c("01100","01110","01120","01200","01210","01212","01230","01232","01250","01240","01220")
map$is_global <- map$suffix %in% GLOBAL

write.csv(map, file.path(OUT,"ko_pathway_map.csv"), row.names=FALSE)
cat(sprintf("\nwrote %s\n  %d KO-pathway links | %d unique KOs | %d unique pathways (%d global/overview)\n",
            file.path(OUT,"ko_pathway_map.csv"), nrow(map),
            length(unique(map$KO)), length(unique(map$pathway)), sum(!duplicated(map$pathway) & map$is_global)))
cat("\nfirst rows:\n"); print(head(map, 4))
quit(save="no")
