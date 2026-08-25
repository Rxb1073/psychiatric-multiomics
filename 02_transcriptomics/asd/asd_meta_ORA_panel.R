#!/usr/bin/env Rscript
# Three-way ORA comparison panel: ALL dysregulated / UP only / DOWN only
# from the confound-free meta-analysis (nominal p<0.05 - EXPLORATORY).
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

sets<-list(
  `All dysregulated` = m$entrez[m$p_FE<0.05],
  `Up in ASD`        = m$entrez[m$p_FE<0.05 & m$beta_FE>0],
  `Down in ASD`      = m$entrez[m$p_FE<0.05 & m$beta_FE<0])
logf("universe=%d | all=%d | up=%d | down=%d",length(uni),
     length(sets[[1]]),length(sets[[2]]),length(sets[[3]]))

for(ont in c("BP","MF","CC")){
  cc<-tryCatch(compareCluster(geneClusters=sets, fun="enrichGO", universe=uni,
                 OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                 pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10,
                 readable=TRUE), error=function(e){message(e$message);NULL})
  n<-if(is.null(cc)) 0 else nrow(as.data.frame(cc))
  logf("GO:%s -> %d enriched terms across the three sets",ont,n)
  if(n>0){
    d<-as.data.frame(cc)
    write.csv(d,file.path(OUT,sprintf("PANEL_ORA_GO%s.csv",ont)),row.names=FALSE)
    print(table(d$Cluster))
    p<-dotplot(cc, showCategory=12, by="count") +
       ggtitle(sprintf("Confound-free meta-analysis - GO:%s (exploratory, nominal p<0.05)",ont)) +
       theme(axis.text.y=element_text(size=8), axis.text.x=element_text(size=9),
             plot.title=element_text(size=11))
    ggsave(file.path(OUT,sprintf("PANEL_ORA_GO%s.png",ont)),p,width=9,height=7.5,dpi=200)
    ggsave(file.path(OUT,sprintf("PANEL_ORA_GO%s.pdf",ont)),p,width=9,height=7.5)
    logf("   wrote PANEL_ORA_GO%s.png/.pdf",ont)
  }
}

cc<-tryCatch(compareCluster(geneClusters=sets, fun="enrichKEGG", universe=uni,
               organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05),
             error=function(e){message(e$message);NULL})
nk<-if(is.null(cc)) 0 else nrow(as.data.frame(cc))
logf("KEGG -> %d enriched pathways across the three sets",nk)
if(nk>0){
  d<-as.data.frame(cc); write.csv(d,file.path(OUT,"PANEL_ORA_KEGG.csv"),row.names=FALSE)
  print(table(d$Cluster))
  p<-dotplot(cc, showCategory=12, by="count") +
     ggtitle("Confound-free meta-analysis - KEGG (exploratory, nominal p<0.05)") +
     theme(axis.text.y=element_text(size=8), plot.title=element_text(size=11))
  ggsave(file.path(OUT,"PANEL_ORA_KEGG.png"),p,width=9,height=7,dpi=200)
  ggsave(file.path(OUT,"PANEL_ORA_KEGG.pdf"),p,width=9,height=7)
}
logf("\nOutput -> %s",OUT)
