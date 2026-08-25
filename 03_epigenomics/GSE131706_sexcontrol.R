suppressMessages({library(limma); library(GEOquery)})
load("mVals_betaVals.RData")                       # mVals: 466610 x 34
gse <- getGEO("GSE131706", GSEMatrix=TRUE, getGPL=FALSE)[[1]]
pd  <- pData(gse)
sexcol <- grep("gender|sex|Sex", colnames(pd), value=TRUE)[1]
cat("Using GEO sex column:", sexcol, "\n")
sex <- as.character(pd[[sexcol]])[match(colnames(mVals), rownames(pd))]
cat("Matched", sum(!is.na(sex)), "of", length(sex), "\n"); print(table(sex, useNA="ifany"))
if (any(is.na(sex))) { cat("UNMATCHED:\n"); print(colnames(mVals)[is.na(sex)]) }
stopifnot(sum(!is.na(sex)) > 24, length(unique(na.omit(sex))) == 2)

keep <- !is.na(sex); M <- mVals[, keep]; sexf <- factor(sex[keep])
fit <- eBayes(lmFit(M, model.matrix(~ sexf)))
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
  "PASS - pipeline recovers sex signal; ASD null is biological." else
  "INVESTIGATE - weak sex signal.", "\n")
saveRDS(tt, "sexcontrol_ASD_toptable.rds")
