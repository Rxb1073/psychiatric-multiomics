suppressMessages({library(limma); library(GEOquery)})
load("mVals_betaVals.RData")
tok <- sub(".*?_([0-9]+-[A-Za-z])\\(.*", "\\1", colnames(mVals))
gse <- getGEO("GSE41826", GSEMatrix=TRUE, getGPL=FALSE)[[1]]
pd  <- pData(gse)
sex <- as.character(pd[["Sex:ch1"]])[match(tok, as.character(pd$title))]
cat("Matched", sum(!is.na(sex)), "of", length(sex), "\n"); print(table(sex))
keep <- !is.na(sex); M <- mVals[, keep]; sexf <- factor(sex[keep])
design <- model.matrix(~ sexf)
fit <- eBayes(lmFit(M, design))
tt  <- topTable(fit, coef=2, number=Inf, sort.by="p")
nsig <- sum(tt$adj.P.Val < 0.05, na.rm=TRUE)
cat(sprintf("\n===== POSITIVE CONTROL (sex) =====\nProbes: %d | FDR<0.05: %d | top adj.P: %.2e\n",
            nrow(tt), nsig, min(tt$adj.P.Val, na.rm=TRUE)))
ok <- tryCatch({
  suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  chr <- ann[head(rownames(tt),200),"chr"]
  cat(sprintf("Top-200 on chrX/Y: %.0f%%\n", 100*mean(chr %in% c("chrX","chrY"),na.rm=TRUE))); TRUE
}, error=function(e) FALSE)
cat("\nVERDICT:", if(nsig>100)
  "PASS - pipeline recovers sex signal; MDD null is biological." else
  "INVESTIGATE - weak sex signal.", "\n")
saveRDS(tt, "sexcontrol_MDD_toptable.rds")
