#!/usr/bin/env Rscript
# Volcano for the CONFOUND-FREE pool (GSE102741 + PRJNA434002, n=45).
# Reads the saved DE table - no DESeq2 recomputation.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(ggplot2); library(ggrepel)})
logf<-function(...) cat(sprintf(...),"\n")
BASE<-"/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
IN  <-file.path(BASE,"POOLED_cortex","DE_ASD_child_confoundfree_full_QC.csv")
PLT <-file.path(BASE,"DISSERTATION_final","01_figures"); if(!dir.exists(PLT)) PLT<-file.path(BASE,"DISSERTATION_final","figures")
dir.create(PLT,showWarnings=FALSE,recursive=TRUE)
FDR<-0.05; LFC<-1
stopifnot(file.exists(IN))
r<-read.csv(IN,stringsAsFactors=FALSE)
v<-r[!is.na(r$padj),]
v$diffexpressed<-"NO"
v$diffexpressed[v$padj<FDR & v$log2FoldChange>= LFC]<-"UP"
v$diffexpressed[v$padj<FDR & v$log2FoldChange<=-LFC]<-"DOWN"
v$diffexpressed<-factor(v$diffexpressed,levels=c("DOWN","UP","NO"))
nUp<-sum(v$diffexpressed=="UP"); nDn<-sum(v$diffexpressed=="DOWN"); nF<-sum(v$padj<FDR)
logf("confound-free: FDR<0.05 = %d | DOWN=%d UP=%d",nF,nDn,nUp)
key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2","PCDH19","PCDHGC3")
lab<-v[v$sym %in% key,]
p<-ggplot(v,aes(log2FoldChange,-log10(padj),color=diffexpressed))+
 geom_point(size=1.1,alpha=.55)+
 scale_color_manual(values=c(DOWN="#2c7bb6",UP="#d7191c",NO="grey72"),name="diffexpressed",drop=FALSE)+
 geom_vline(xintercept=c(-LFC,LFC),linetype=2,color="grey45")+
 geom_hline(yintercept=-log10(FDR),linetype=2,color="grey45")+
 geom_text_repel(data=lab,aes(label=sym),size=3,max.overlaps=30,show.legend=FALSE,seed=1,
                 min.segment.length=0,box.padding=0.5,color="grey20")+
 labs(title="ASD vs Control - confound-free pool (GSE102741 + PRJNA434002)",
      subtitle=sprintf("childhood cortex, n=45 (16 ASD / 29 Control); age balanced without matching (p=0.678)\n%d genes at FDR<0.05  |  candidate genes labelled, all non-significant",nF),
      x="log2 fold change (ASD / Control)",y="-log10 adjusted p")+
 theme_bw(base_size=12)+theme(legend.position="right",plot.title=element_text(size=12.5),
                              plot.subtitle=element_text(size=10,lineheight=1.05))
ggsave(file.path(PLT,"volcano_ASD_child_confoundfree.png"),p,width=8.3,height=6.5,dpi=200,type="cairo")
ggsave(file.path(PLT,"volcano_ASD_child_confoundfree.pdf"),p,width=8.3,height=6.5)
logf("wrote volcano_ASD_child_confoundfree.png/.pdf -> %s",PLT)

# effect-size comparison: confounded cohort vs confound-free, on the candidate genes
G<-file.path(BASE,"POOLED_cortex","DE_ASD_child_GSE59288only_full.csv")
if(file.exists(G)){
 g<-read.csv(G,stringsAsFactors=FALSE)
 m<-merge(g[g$sym%in%key,c("sym","log2FoldChange","padj")],
          v[v$sym%in%key,c("sym","log2FoldChange","padj")],by="sym",suffixes=c("_confounded","_clean"))
 m$fold_shrink<-round(abs(m$log2FoldChange_confounded/m$log2FoldChange_clean),1)
 m<-m[order(-abs(m$log2FoldChange_confounded)),]
 write.csv(m,file.path(BASE,"DISSERTATION_final","tables","effectsize_confounded_vs_confoundfree.csv"),row.names=FALSE)
 cat("\n===== EFFECT SIZE COLLAPSE (candidate genes) =====\n"); print(m,row.names=FALSE)
}
