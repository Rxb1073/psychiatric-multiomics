.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
D  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
N2 <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2"
OUT<- file.path(D,"Microbiome_ASD","tables")
lt <- function(x){ x<-gsub('"','',as.character(x)); x<-sub(".*[;.]","",x); sub("^g__","",x) }
rd <- function(p,s="\t"){ if(!file.exists(p)){cat("MISSING:",p,"\n");return(NULL)}
  read.delim(p,sep=s,check.names=FALSE,stringsAsFactors=FALSE) }
grab <- function(df,arm,meth,g,e,pp,q){
  if(meth=="MaAsLin2" && "metadata" %in% names(df)) df <- df[df$metadata=="diagnosis" | df$metadata=="disease",]
  if(meth=="ANCOM-BC2") pp <- NA
  if(is.null(df)) return(NULL)
  miss <- setdiff(na.omit(c(g,e,pp,q)), names(df))
  if(length(miss)){ cat(sprintf("  [%s %s] missing: %s\n     available: %s\n",
    arm,meth,paste(miss,collapse=","),paste(names(df),collapse=", "))); return(NULL) }
  data.frame(arm=arm, method=meth, genus=lt(df[[g]]),
    effect=round(as.numeric(df[[e]]),4), p=if(is.na(pp)) NA_real_ else signif(as.numeric(df[[pp]]),3),
    q=signif(as.numeric(df[[q]]),3), stringsAsFactors=FALSE) }
files <- list(
  SCZ_ancom=file.path(D,"Microbiome_SCZ/tables/ancombc2_SZ_vs_HC.tsv"),
  SCZ_maas =file.path(D,"Microbiome_SCZ/tables/maaslin2_SZ_vs_HC.tsv"),
  MDD_ancom=file.path(D,"Microbiome_MDD/tables/ancombc2_MDD_vs_Control.tsv"),
  MDD_maas =file.path(N2,"MDD/microbiome/merged_genus/maaslin2_output/all_results.tsv"),
  ASD_ancom=file.path(D,"Microbiome_ASD/tables/stage5_ancombc2_all.csv"),
  ASD_maas =file.path(D,"Microbiome_ASD/tables/maaslin2_out/all_results.tsv"))
cat("=== column names ===\n")
dfs <- lapply(names(files), function(n){ f<-files[[n]]
  d <- if(grepl("\\.csv$",f)) rd(f,",") else rd(f,"\t")
  if(!is.null(d)) cat(sprintf("  %-10s %3d rows | %s\n",n,nrow(d),paste(names(d),collapse=", ")))
  d })
names(dfs) <- names(files)
cat("\n=== assembling ===\n")
parts <- list(
  grab(dfs$SCZ_ancom,"Schizophrenia","ANCOM-BC2","genus","log2FC_SZ_vs_HC","q_value","q_value"),
  grab(dfs$SCZ_maas ,"Schizophrenia","MaAsLin2","feature","coef","pval","qval"),
  grab(dfs$MDD_ancom,"Depression","ANCOM-BC2","genus","log2FC_MDD_vs_Control","q_value","q_value"),
  grab(dfs$MDD_maas ,"Depression","MaAsLin2","feature","coef","pval","qval"),
  grab(dfs$ASD_ancom,"Autism","ANCOM-BC2","taxon","lfc","q","q"),
  grab(dfs$ASD_maas ,"Autism","MaAsLin2","feature","coef","pval","qval"))
S <- do.call(rbind, Filter(Negate(is.null), parts))
if(is.null(S)){ cat("\nnothing assembled — check columns above\n") } else {
  S <- S[order(S$arm,S$method,S$q),]
  cat(sprintf("\nrows: %d\n",nrow(S))); print(table(S$arm,S$method))
  f <- file.path(OUT,"SupplementaryTableS4_differential_abundance_both_methods.csv")
  write.csv(S,f,row.names=FALSE); cat("\nwritten:",f,"\n")
  print(head(S,8), row.names=FALSE) }
