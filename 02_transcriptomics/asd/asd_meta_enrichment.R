#!/usr/bin/env Rscript
# Enrichment on the CONFOUND-FREE META-ANALYSIS ranked list.
# Primary: GSEA on genes ranked by signed -log10(fixed-effect p).
# Secondary: ORA on nominally significant (p<0.05) genes - reported
#            with the caveat that these do not survive FDR.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
logf<-function(...) cat(sprintf(...),"\n")
BASE<-"/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
IN  <-file.path(BASE,"POOLED_cortex","DE_ASD_child_metaanalysis_confoundfree.csv")
OUT <-file.path(BASE,"DISSERTATION_final","enrichment"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
stopifnot(file.exists(IN))

m<-read.csv(IN,stringsAsFactors=FALSE)
m$entrez<-suppressMessages(mapIds(org.Hs.eg.db,m$ens,"ENTREZID","ENSEMBL"))
m<-m[!is.na(m$entrez) & !duplicated(m$entrez) & is.finite(m$p_FE), ]
logf("genes with Entrez ID: %d",nrow(m))

# ---------- PRIMARY: GSEA on signed -log10(p) ----------
m$rankstat <- sign(m$beta_FE) * -log10(pmax(m$p_FE, .Machine$double.xmin))
gl <- sort(setNames(m$rankstat, m$entrez), decreasing=TRUE)
gl <- gl[is.finite(gl)]
logf("ranked list: %d genes (range %.2f to %.2f)",length(gl),min(gl),max(gl))

for(ont in c("BP","MF","CC")){
  g<-tryCatch(gseGO(geneList=gl,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont=ont,
                    pAdjustMethod="BH",pvalueCutoff=0.05,minGSSize=15,maxGSSize=500,
                    verbose=FALSE,seed=TRUE),error=function(e){message(e$message);NULL})
  n<-if(is.null(g)) 0 else nrow(as.data.frame(g))
  logf("GSEA GO:%s -> %d significant terms",ont,n)
  if(n>0){
    d<-as.data.frame(g); d$direction<-ifelse(d$NES>0,"up in ASD","down in ASD")
    write.csv(d,file.path(OUT,sprintf("META_GSEA_GO%s.csv",ont)),row.names=FALSE)
    cat(sprintf("   top: %s\n",paste(head(d$Description,6),collapse=" | ")))
    if(n>=3){
      p<-dotplot(g,showCategory=min(15,n),split=".sign")+facet_grid(.~.sign)+
         ggtitle(sprintf("Meta-analysis (confound-free) - GSEA GO:%s",ont))+
         theme(axis.text.y=element_text(size=8))
      ggsave(file.path(OUT,sprintf("META_GSEA_GO%s.png",ont)),p,width=9.5,height=6.5,dpi=200,type="cairo")
    }
  }
}
k<-tryCatch(gseKEGG(geneList=gl,organism="hsa",pAdjustMethod="BH",pvalueCutoff=0.05,
                    minGSSize=15,maxGSSize=500,verbose=FALSE,seed=TRUE),error=function(e)NULL)
nk<-if(is.null(k)) 0 else nrow(as.data.frame(k))
logf("GSEA KEGG -> %d significant pathways",nk)
if(nk>0){d<-as.data.frame(k);write.csv(d,file.path(OUT,"META_GSEA_KEGG.csv"),row.names=FALSE)
  cat(sprintf("   top: %s\n",paste(head(d$Description,6),collapse=" | ")))}

# ---------- SECONDARY: ORA on nominally significant genes ----------
cat("\n--- SECONDARY: ORA on nominal p<0.05 (NOT FDR-significant) ---\n")
uni<-m$entrez
nomD<-m$entrez[m$p_FE<0.05 & m$beta_FE<0]
nomU<-m$entrez[m$p_FE<0.05 & m$beta_FE>0]
logf("nominal p<0.05: down=%d up=%d",length(nomD),length(nomU))
for(s in list(list(g=nomD,nm="down"),list(g=nomU,nm="up"))){
  if(length(s$g)<10){logf("  [%s] too few (%d)",s$nm,length(s$g));next}
  for(ont in c("BP","MF","CC")){
    e<-tryCatch(enrichGO(gene=s$g,universe=uni,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont=ont,
                         pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.10,readable=TRUE),
                error=function(e)NULL)
    n<-if(is.null(e)) 0 else nrow(as.data.frame(e))
    logf("  ORA GO:%s [%s] -> %d terms",ont,s$nm,n)
    if(n>0){d<-as.data.frame(e)
      write.csv(d,file.path(OUT,sprintf("META_ORAnominal_GO%s_%s.csv",ont,s$nm)),row.names=FALSE)
      cat(sprintf("      top: %s\n",paste(head(d$Description,5),collapse=" | ")))}
  }
}
logf("\nOutput -> %s",OUT)
