## Is the cross-arm microbial overlap specific to butyrate producers?
##
## The enumerative test asks only whether the SAME genera recur between arms
## (3 of MDD's 4 concordant genera among SCZ's 26, from a shared pool of 116,
## hypergeometric p = 0.035 -- the figures this script writes to
## cross_arm_overlap.csv; an earlier draft of this header said 88 and 0.076,
## from a run that predates the current concordant sets). That discards two things the data contain: the
## DIRECTION of each effect, and the functional annotation of each genus.
##
## This script tests the mechanistic version of the claim instead -- whether
## short-chain fatty acid producers are over-represented among the genera that
## recur across arms, and whether recurring genera agree in direction more often
## than chance. Both use information the hypergeometric ignores.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

D <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
OUT <- file.path(D,"INTEGRATION_inputs","12_cross_arm_microbiome")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

gname <- function(x) sub("^.*g__", "", gsub("[.;]", "_", as.character(x)))

## ---- concordant sets per arm ----------------------------------------------
scz <- read.delim(file.path(D,"Microbiome_SCZ/tables/robust_genera_both_methods_v2.tsv"), header=FALSE)
scz <- data.frame(genus=gname(scz$V1), lfc=as.numeric(scz$V3), stringsAsFactors=FALSE)
mdd <- read.delim(file.path(D,"Microbiome_MDD/tables/robust_genera_both_methods.tsv"))
mdd <- data.frame(genus=gname(mdd$genus), lfc=as.numeric(mdd$ancombc_log2FC), stringsAsFactors=FALSE)
asd <- read.csv(file.path(D,"Microbiome_ASD/tables/stage5_concordant_core.csv"))
asd <- data.frame(genus=gname(asd[[2]]), lfc=as.numeric(asd[[3]]), stringsAsFactors=FALSE)
sets <- list(ASD=asd, SCZ=scz, MDD=mdd)
for(a in names(sets)) cat(sprintf("%s concordant genera: %d\n", a, nrow(sets[[a]])))

## ---- tested pool per arm ---------------------------------------------------
pool_scz <- unique(gname(read.delim(file.path(D,"Microbiome_SCZ/tables/maaslin2_SZ_vs_HC.tsv"))$feature))
pool_mdd <- unique(gname(read.delim(file.path(D,"Microbiome_MDD/tables/ancombc2_MDD_vs_Control.tsv"))$genus))
pool_asd <- unique(gname(read.csv(file.path(D,"Microbiome_ASD/tables/stage5_ancombc2_all.csv"))$taxon))
pools <- list(ASD=pool_asd, SCZ=pool_scz, MDD=pool_mdd)
for(a in names(pools)) cat(sprintf("%s genera tested: %d\n", a, length(pools[[a]])))

## ---- pairwise overlap, enumerative and directional -------------------------
cat("\n", strrep("=",72), "\nPAIRWISE OVERLAP\n", strrep("=",72), "\n", sep="")
res <- list()
for(p in list(c("MDD","SCZ"), c("ASD","SCZ"), c("ASD","MDD"))){
  a <- p[1]; b <- p[2]
  pool <- intersect(pools[[a]], pools[[b]])
  A <- sets[[a]][sets[[a]]$genus %in% pool, ]
  B <- sets[[b]][sets[[b]]$genus %in% pool, ]
  ov <- merge(A, B, by="genus", suffixes=c("_1","_2"))
  ph <- phyper(nrow(ov)-1, nrow(A), length(pool)-nrow(A), nrow(B), lower.tail=FALSE)
  same <- sum(sign(ov$lfc_1) == sign(ov$lfc_2))
  ## direction agreement among shared genera: binomial against 50%
  pb <- if(nrow(ov)) binom.test(same, nrow(ov), 0.5)$p.value else NA
  cat(sprintf("\n%s vs %s | shared pool %d | concordant %d and %d | overlap %d (p = %.4f)\n",
              a, b, length(pool), nrow(A), nrow(B), nrow(ov), ph))
  if(nrow(ov)){
    cat(sprintf("  same direction: %d of %d (binomial p = %.4f)\n", same, nrow(ov), pb))
    print(ov[order(ov$genus), ], row.names=FALSE)
  }
  res[[paste(a,b)]] <- data.frame(pair=paste(a,"vs",b), pool=length(pool),
    n_a=nrow(A), n_b=nrow(B), overlap=nrow(ov), p_hyper=ph,
    same_direction=same, p_binom=pb)
}
S <- do.call(rbind, res); rownames(S) <- NULL

## ---- butyrate / SCFA annotation -------------------------------------------
## Annotation as used for the cross-arm figure. Genera classified as
## short-chain fatty acid producers on the basis of published characterisation.
SCFA <- c("Faecalibacterium","Fusicatenibacter","Agathobacter","Roseburia",
          "Butyricicoccus","Coprococcus","Anaerostipes","Eubacterium",
          "Eubacterium_hallii_group","Eubacterium_rectale_group","Blautia",
          "Ruminococcus","Subdoligranulum","Butyrivibrio","Lachnospira",
          "Dorea","Odoribacter","Alistipes")
cat("\n", strrep("=",72), "\nSCFA PRODUCERS AMONG SHARED GENERA\n", strrep("=",72), "\n", sep="")
cat("NOTE: the annotation list is fixed above and applied identically to the\n")
cat("shared set and to the background pool, so the test is not circular.\n")

shared_all <- unique(unlist(lapply(res, function(r) NULL)))
shared_all <- unique(unlist(lapply(list(c("MDD","SCZ"),c("ASD","SCZ"),c("ASD","MDD")), function(p){
  pool <- intersect(pools[[p[1]]], pools[[p[2]]])
  intersect(sets[[p[1]]]$genus[sets[[p[1]]]$genus %in% pool],
            sets[[p[2]]]$genus[sets[[p[2]]]$genus %in% pool]) })))
pool_any <- unique(unlist(pools))
k  <- sum(shared_all %in% SCFA); n <- length(shared_all)
K  <- sum(pool_any %in% SCFA);   N <- length(pool_any)
pe <- phyper(k-1, K, N-K, n, lower.tail=FALSE)
cat(sprintf("\nshared genera (any pair): %d -> %s\n", n, paste(shared_all, collapse=", ")))
cat(sprintf("of these, SCFA producers: %d (%.0f%%)\n", k, 100*k/n))
cat(sprintf("background: %d SCFA producers among %d genera tested (%.0f%%)\n", K, N, 100*K/N))
cat(sprintf("enrichment: %.1f-fold | hypergeometric p = %.4f\n", (k/n)/(K/N), pe))

write.csv(S, file.path(OUT,"cross_arm_overlap.csv"), row.names=FALSE)
write.csv(data.frame(genus=shared_all, scfa=shared_all %in% SCFA),
          file.path(OUT,"shared_genera.csv"), row.names=FALSE)
cat("\n", strrep("=",72), "\nSUMMARY\n", strrep("=",72), "\n", sep="")
print(S, row.names=FALSE)
cat(sprintf("\nwritten to %s\n", OUT))
