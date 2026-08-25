suppressMessages({library(limma); library(GEOquery)})
load("mVals_betaVals.RData")                       # mVals: 466610 x 34
pd <- pData(getGEO("GSE131706", GSEMatrix=TRUE, getGPL=FALSE)[[1]])

## parse "X years Y days" -> numeric years
raw <- as.character(pd[["age:ch1"]])[match(colnames(mVals), rownames(pd))]
yrs <- as.numeric(sub(".*?([0-9]+)\\s*years.*", "\\1", raw))
dys <- ifelse(grepl("days", raw), as.numeric(sub(".*?([0-9]+)\\s*days.*", "\\1", raw)), 0)
age <- yrs + dys/365.25
cat("Matched age for", sum(!is.na(age)), "of", length(age), "samples\n")
cat("Age range:", round(range(age, na.rm=TRUE),2), " | SD:", round(sd(age, na.rm=TRUE),2), "\n")
stopifnot(sum(!is.na(age)) > 24, sd(age, na.rm=TRUE) > 0)

keep <- !is.na(age); M <- mVals[, keep]; a <- age[keep]
fit <- eBayes(lmFit(M, model.matrix(~ a)))          # continuous age contrast
tt  <- topTable(fit, coef=2, number=Inf, sort.by="p")
nsig <- sum(tt$adj.P.Val < 0.05, na.rm=TRUE)
cat(sprintf("\n===== POSITIVE CONTROL (age) =====\nProbes: %d | FDR<0.05: %d | top adj.P: %.2e\n",
            nrow(tt), nsig, min(tt$adj.P.Val, na.rm=TRUE)))

## sanity: do we recover canonical age-associated CpGs? (ELOVL2, FHL2 are textbook)
clock <- c("cg16867657","cg21572722","cg24724428",  # ELOVL2 region
           "cg06639320","cg22454769")                # FHL2 region
hit <- clock[clock %in% rownames(tt)]
if(length(hit)) { cat("\nKnown age CpGs (ELOVL2/FHL2) adj.P.Val:\n"); print(round(tt[hit,"adj.P.Val"],4)) }

cat("\nVERDICT:", if(nsig>50)
  "PASS - pipeline recovers age signal; ASD null is biological." else
  "INVESTIGATE - weak age signal (note n=34 limits power even for age).", "\n")
saveRDS(tt, "agecontrol_ASD_toptable.rds")
