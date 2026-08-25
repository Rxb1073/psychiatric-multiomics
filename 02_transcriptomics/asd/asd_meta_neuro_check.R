#!/usr/bin/env Rscript
# Is there a neuronal/synaptic subset inside the 1,092 downregulated genes
# from the confound-free meta-analysis, or only the few candidate hits?
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
set.seed(42)
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)})
logf<-function(...) cat(sprintf(...),"\n")
BASE<-"/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
IN  <-file.path(BASE,"POOLED_cortex","DE_ASD_child_metaanalysis_confoundfree.csv")
OUT <-file.path(BASE,"DISSERTATION_final","enrichment")
m<-read.csv(IN,stringsAsFactors=FALSE)
m$entrez<-suppressMessages(mapIds(org.Hs.eg.db,m$ens,"ENTREZID","ENSEMBL"))
m<-m[!is.na(m$entrez)&!duplicated(m$entrez)&is.finite(m$p_FE),]
uni<-m$entrez; dn<-m[m$p_FE<0.05 & m$beta_FE<0,]
logf("universe=%d | down set=%d",length(uni),nrow(dn))

# GO terms of interest (neuronal / synaptic / interneuron)
gos<-c("GO:0007399"="nervous system development",
       "GO:0030182"="neuron differentiation",
       "GO:0048666"="neuron development",
       "GO:0007268"="chemical synaptic transmission",
       "GO:0050808"="synapse organization",
       "GO:0045202"="synapse",
       "GO:0001764"="neuron migration",
       "GO:0007411"="axon guidance",
       "GO:0021885"="forebrain cell migration",
       "GO:0097485"="neuron projection guidance")
cat("\n===== neuronal gene sets: how many are in the DOWN list? =====\n")
cat(sprintf("%-34s %6s %8s %8s %8s\n","GO term","inGO","inDown","expected","ratio"))
res<-list()
for(id in names(gos)){
  g<-tryCatch(unique(AnnotationDbi::select(org.Hs.eg.db,keys=id,keytype="GOALL",columns="ENTREZID")$ENTREZID),
              error=function(e)character(0))
  g<-intersect(g,uni)
  hit<-intersect(g,dn$entrez)
  exp<-length(g)*nrow(dn)/length(uni)
  cat(sprintf("%-34s %6d %8d %8.1f %8.2f\n",gos[id],length(g),length(hit),exp,length(hit)/max(exp,0.01)))
  res[[gos[id]]]<-hit
}
# which specific neuro genes are in the down list, ranked by p
allneuro<-unique(unlist(res))
sub<-dn[dn$entrez %in% allneuro,]
sub<-sub[order(sub$p_FE),]
logf("\ntotal distinct neuronal genes in down list: %d of %d (%.1f%%)",
     nrow(sub),nrow(dn),100*nrow(sub)/nrow(dn))
cat("\n===== top 30 neuronal genes in the down list =====\n")
print(head(sub[,c("sym","beta_1","beta_2","beta_FE","p_FE","padj_FE")],30),row.names=FALSE)
write.csv(sub[,c("sym","ens","beta_1","beta_2","beta_FE","p_FE","padj_FE","I2")],
          file.path(OUT,"META_neuronal_genes_in_downset.csv"),row.names=FALSE)
logf("\nwrote META_neuronal_genes_in_downset.csv")
