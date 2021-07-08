# Date: 06.07.2021
# Author: Anastasiya Boersch
# E-mail: anastasiya.boersch@unibas.ch

# R version 3.6.0 (2019-04-26)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS Linux 7 (Core)

# Load libraries
library(matrixStats) # matrixStats_0.54.0
library(grid) # base package
library(gridExtra) # gridExtra_2.3
library(ggplot2) # ggplot2_3.1.1

# Load Table_S3
table_s3 <- read.delim("Figure_S3/Table_S3.txt",sep="\t",header=T)

# Cumulative distribution of logFCs for targets and non-targets for proteomics
# CFIm25 KD
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"Prot_log2FC_25_KD"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"Prot_log2FC_25_KD"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF125_KD=c(tar_lf,nontar_lf))

p1<-ggplot(df, aes(x=log2ratio_CF125_KD, color=type)) +
  stat_ecdf() + xlim(-0.5,0.5) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(25-KD/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,2)),size=3, family="Helvetica")
p1

# CFIm68 KD
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"Prot_log2FC_68_KD"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"Prot_log2FC_68_KD"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF168_KD=c(tar_lf,nontar_lf))

p2<-ggplot(df, aes(x=log2ratio_CF168_KD, color=type)) +
  stat_ecdf()+ xlim(-0.5,0.5) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(68-KD/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,5)),size=3, family="Helvetica")
p2

# CFIm25 OE
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"Prot_log2FC_25_OE"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"Prot_log2FC_25_OE"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF125_OE=c(tar_lf,nontar_lf))

p3<-ggplot(df, aes(x=log2ratio_CF125_OE, color=type)) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(25-OE/Control)")+
  ylab("Cumulative distribution")+
  stat_ecdf()+ xlim(-0.5,0.5) +
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,2)),size=3, family="Helvetica")
p3

# CFIm68 OE
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"Prot_log2FC_68_OE"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"Prot_log2FC_68_OE"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF168_OE=c(tar_lf,nontar_lf))

p4<-ggplot(df, aes(x=log2ratio_CF168_OE, color=type)) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(68-OE/Control)")+
  ylab("Cumulative distribution")+
  stat_ecdf()+ xlim(-0.5,0.5) +
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,2)),size=3, family="Helvetica")
p4

pdf("Figure_S3/Cumulative_distribution_proteomics_targets_nontargets.pdf",height=4,width=7.8)
grid.arrange(p1,p2,p3,p4,widths = c(0.5,0.5),layout_matrix = rbind(c(1,2),c(3,4)))
dev.off()

# Make cumulative distribution of the expression changes in targets/not targets
# CFIm25 KD
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"RNA_log2FC_25_KD"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"RNA_log2FC_25_KD"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF125_KD=c(tar_lf,nontar_lf))

p1<-ggplot(df, aes(x=log2ratio_CF125_KD, color=type)) +
  stat_ecdf() + xlim(-1,1) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(25-KD/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,8)),size=3, family="Helvetica")
p1

# CFIm68 KD
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"RNA_log2FC_68_KD"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"RNA_log2FC_68_KD"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF168_KD=c(tar_lf,nontar_lf))

p2<-ggplot(df, aes(x=log2ratio_CF168_KD, color=type)) +
  stat_ecdf() + xlim(-1,1) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(68-KD/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,12)),size=3, family="Helvetica")
p2

# CFIm25 OE
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"RNA_log2FC_25_OE"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"RNA_log2FC_25_OE"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF125_OE=c(tar_lf,nontar_lf))

p3<-ggplot(df, aes(x=log2ratio_CF125_OE, color=type)) +
  stat_ecdf() + xlim(-1,1) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(25-OE/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,2)),size=3, family="Helvetica")
p3

# CFIm68 OE
tar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==1),"RNA_log2FC_68_OE"])
nontar_lf <- na.omit(table_s3[which(table_s3$CFIm_Target==0),"RNA_log2FC_68_OE"])
t <- wilcox.test(tar_lf,nontar_lf)

df <- data.frame(type=c(rep("target",length(tar_lf)),rep("nontarget",length(nontar_lf))),log2ratio_CF168_OE=c(tar_lf,nontar_lf))

p4<-ggplot(df, aes(x=log2ratio_CF168_OE, color=type)) +
  stat_ecdf() + xlim(-1,1) +
  theme_bw()+
  scale_color_manual(values=c("red","darkblue"))+
  theme(legend.title = element_blank(),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black", size=8, family="Helvetica"),
        axis.text.y = element_text(color="black", size=8, family="Helvetica"),
        axis.title.x = element_text(color="black", size=8, family="Helvetica"),
        axis.title.y = element_text(color="black", size=8, family="Helvetica"))+
  xlab("log2(68-OE/Control)")+
  ylab("Cumulative distribution")+
  annotate("text", x=-0.3, y=0.8, label= paste0("p=",round(t$p.value,2)),size=3, family="Helvetica")
p4

pdf("Figure_S3/Cumulative_distribution_transcriptomics_targets_nontargets.pdf",height=4,width=7.8)
grid.arrange(p1,p2,p3,p4,widths = c(0.5,0.5),layout_matrix = rbind(c(1,2),c(3,4)))
dev.off()
