#!/usr/bin/env Rscript
# ORA on the COMBINED nominally-significant gene set (p<0.05, both directions)
# from the confound-free meta-analysis. Exploratory only - these genes do not
# survive FDR correction, and ~half the list is expected by chance.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
set.seed(42)
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
logf<-function(...) cat(sprintf(...),"\n")
BASE<-"/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
IN  <-file.path(BASE,"POOLED_cortex","DE_ASD_child_metaanalysis_confoundfree.csv")
OUT <-file.path(BASE,"DISSERTATION_final","enrichment"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
stopifnot(file.exists(IN))

m<-read.csv(IN,stringsAsFactors=FALSE)
m$entrez<-suppressMessages(mapIds(org.Hs.eg.db,m$ens,"ENTREZID","ENSEMBL"))
m<-m[!is.na(m$entrez) & !duplicated(m$entrez) & is.finite(m$p_FE), ]
uni<-m$entrez
comb<-m$entrez[m$p_FE<0.05]
logf("universe=%d | combined nominal p<0.05 = %d (expected by chance ~%.0f)",
     length(uni),length(comb),0.05*length(uni))

for(ont in c("BP","MF","CC")){
  e<-tryCatch(enrichGO(gene=comb,universe=uni,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont=ont,
                       pAdjustMethod="BH",pvalueCutoff=0.05,qvalueCutoff=0.10,readable=TRUE),
              error=function(e){message(e$message);NULL})
  n<-if(is.null(e)) 0 else nrow(as.data.frame(e))
  logf("ORA GO:%s [combined] -> %d terms",ont,n)
  if(n>0){
    d<-as.data.frame(e)
    write.csv(d,file.path(OUT,sprintf("META_ORAcombined_GO%s.csv",ont)),row.names=FALSE)
    cat("   top 8:\n")
    print(head(d[,c("Description","GeneRatio","BgRatio","p.adjust","Count")],8),row.names=FALSE)
    if(n>=3){
      p<-dotplot(e,showCategory=min(15,n))+ggtitle(sprintf("Meta-analysis combined nominal set - GO:%s",ont))+
         theme(axis.text.y=element_text(size=8))
      ggsave(file.path(OUT,sprintf("META_ORAcombined_GO%s.png",ont)),p,width=8.5,height=6.5,dpi=200,type="cairo")
    }
  }
}
k<-tryCatch(enrichKEGG(gene=comb,universe=uni,organism="hsa",pAdjustMethod="BH",pvalueCutoff=0.05),
            error=function(e)NULL)
nk<-if(is.null(k)) 0 else nrow(as.data.frame(k))
logf("ORA KEGG [combined] -> %d pathways",nk)
if(nk>0){d<-as.data.frame(k);write.csv(d,file.path(OUT,"META_ORAcombined_KEGG.csv"),row.names=FALSE)
  cat("   top: ",paste(head(d$Description,8),collapse=" | "),"\n")}

# how much of the combined result is just the down set?
dn<-m$entrez[m$p_FE<0.05 & m$beta_FE<0]
logf("\nfor reference: down=%d (%.0f%% of combined), up=%d",
     length(dn),100*length(dn)/length(comb),length(comb)-length(dn))
