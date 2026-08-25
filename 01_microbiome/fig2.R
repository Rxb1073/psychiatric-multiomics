library(ggplot2)
d <- data.frame(var=rep(c("Study (batch)","Diagnosis"),each=2),
  when=factor(rep(c("Before","After"),2),levels=c("Before","After")),
  R2=c(0.3927,0.1212,0.0074,0.0077))
ggplot(d, aes(when,R2,fill=var)) + geom_col(position="dodge",width=0.7) +
  geom_text(aes(label=sprintf("%.3f",R2)),position=position_dodge(0.7),vjust=-0.4,size=3.5) +
  scale_fill_manual(values=c("Study (batch)"="#e67e22","Diagnosis"="#2980b9")) +
  labs(y="PERMANOVA R2", x=NULL, fill=NULL) + theme_minimal(base_size=13)
ggsave("fig2_batch.png", width=5, height=4, dpi=300)
