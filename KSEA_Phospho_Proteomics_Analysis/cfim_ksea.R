library(devtools)
#devtools::install_github("evocellnet/ksea")
library(ksea)
library(stringr)
library(ComplexHeatmap)
library(factoextra)
library(NbClust)


############OVEREXPRESSION RESULTS############

###M25 vs Control###
data_clean <- read.table("input_for_ksea_cfim_m25_control_OE.txt", header = TRUE)
peptides=data_clean$Peptide_ID
peptides=as.character(peptides)
#sites_logFC=data_clean$Log2FC
sites=data_clean$Log2FC
kinases=data_clean$Predicted_Kinase
kinases=as.character(kinases)
unique_kinases <-unique(kinases)
ind <- which (kinases %in% unique_kinases[1])
regulons <- list(PKCE=peptides[ind])

for (i in 1:length(unique_kinases)) {
  ind <- which (kinases %in% unique_kinases[i])
  regulons[[unique_kinases[i]]] <-peptides[ind]
}
names(sites) <- peptides
#names(sites_logFC) <- peptides
sites <- sites[order(sites, decreasing=TRUE)]
options("mc.cores" = 16)
kinases_ksea_kd_cfim <- ksea_batchKinases(names(sites), sites, regulons, trial=10000)
logFC_m25_cont=c()
for (i in 1:length(unique_kinases)) {
  logFC_m25_cont=c(logFC_m25_cont, as.numeric(mean(sites[regulons[[unique_kinases[i]]]])))
}
result_OE_m25_cont =data.frame("pval"=kinases_ksea_kd_cfim, "ES"=logFC_m25_cont, row.names = unique_kinases)
write.csv(result_OE_m25_cont, file = 'logFC_m25_cont_OE')
###M68 vs Control###
data_clean <- read.table("input_for_ksea_cfim_m68_control_OE.txt", header = TRUE)
peptides=data_clean$Peptide_ID
peptides=as.character(peptides)
#sites_logFC=data_clean$Log2FC
sites=data_clean$Log2FC
kinases=data_clean$Predicted_Kinase
kinases=as.character(kinases)
unique_kinases <-unique(kinases)
ind <- which (kinases %in% unique_kinases[1])
regulons <- list(LRRK2=peptides[ind])

for (i in 1:length(unique_kinases)) {
  ind <- which (kinases %in% unique_kinases[i])
  regulons[[unique_kinases[i]]] <-peptides[ind]
}
names(sites) <- peptides
#names(sites_logFC) <- peptides
sites <- sites[order(sites, decreasing=TRUE)]
options("mc.cores" = 16)
kinases_ksea_kd_cfim <- ksea_batchKinases(names(sites), sites, regulons, trial=10000)
logFC_m68_cont=c()
for (i in 1:length(unique_kinases)) {
  logFC_m68_cont=c(logFC_m68_cont, as.numeric(mean(sites[regulons[[unique_kinases[i]]]])))
}
result_OE_m68_cont =data.frame("pval"=kinases_ksea_kd_cfim, "ES"=logFC_m68_cont, row.names = unique_kinases)
write.csv(result_OE_m68_cont, file = 'logFC_m68_cont_OE')

############KNOCK-DOWN RESULTS############

###M25 vs Control###

data_clean <- read.table("input_for_ksea_cfim_m25_control_KD.txt", header = TRUE)
peptides=data_clean$Peptide_ID
peptides=as.character(peptides)
#sites_logFC=data_clean$Log2FC
sites=data_clean$Log2FC
kinases=data_clean$Predicted_Kinase
kinases=as.character(kinases)
unique_kinases <-unique(kinases)
ind <- which (kinases %in% unique_kinases[1])
regulons <- list(LRRK2=peptides[ind])

for (i in 1:length(unique_kinases)) {
  ind <- which (kinases %in% unique_kinases[i])
  regulons[[unique_kinases[i]]] <-peptides[ind]
}
names(sites) <- peptides
#names(sites_logFC) <- peptides
sites <- sites[order(sites, decreasing=TRUE)]
options("mc.cores" = 16)
kinases_ksea_kd_cfim <- ksea_batchKinases(names(sites), sites, regulons, trial=10000)
logFC_m25_cont=c()
for (i in 1:length(unique_kinases)) {
  logFC_m25_cont=c(logFC_m25_cont, as.numeric(mean(sites[regulons[[unique_kinases[i]]]])))
}
result_KD_m25_cont =data.frame("pval"=kinases_ksea_kd_cfim, "ES"=logFC_m25_cont, row.names = unique_kinases)
write.csv(result_KD_m25_cont, file = 'logFC_m25_cont_KD')

###M68 vs Control###
data_clean <- read.table("input_for_ksea_cfim_m68_control_KD.txt", header = TRUE)
peptides=data_clean$Peptide_ID
peptides=as.character(peptides)
#sites_logFC=data_clean$Log2FC
sites=data_clean$Log2FC
kinases=data_clean$Predicted_Kinase
kinases=as.character(kinases)
unique_kinases <-unique(kinases)
ind <- which (kinases %in% unique_kinases[1])
regulons <- list(LRRK2=peptides[ind])

for (i in 1:length(unique_kinases)) {
  ind <- which (kinases %in% unique_kinases[i])
  regulons[[unique_kinases[i]]] <-peptides[ind]
}
names(sites) <- peptides
#names(sites_logFC) <- peptides
sites <- sites[order(sites, decreasing=TRUE)]
options("mc.cores" = 16)
kinases_ksea_kd_cfim <- ksea_batchKinases(names(sites), sites, regulons, trial=10000)
logFC_m68_cont=c()
for (i in 1:length(unique_kinases)) {
  logFC_m68_cont=c(logFC_m68_cont, as.numeric(mean(sites[regulons[[unique_kinases[i]]]])))
}
result_KD_m68_cont =data.frame("pval"=kinases_ksea_kd_cfim, "ES"=logFC_m68_cont, row.names = unique_kinases)
write.csv(result_KD_m68_cont, file = 'logFC_m68_cont_KD')

kkk=1
to_keep=list()
for (i in 1:length(unique_kinases)) {
  ind <- which (kinases %in% unique_kinases[i])
  if (length(ind)>10){
    to_keep[kkk]=i
    kkk=kkk+1
  }
}

all_kinases = unique_kinases
to_keep_vector <- unlist(to_keep)
all_kinases=all_kinases[to_keep_vector]

col_names <- c("result_OE_m25_cont", "result_KD_m25_cont", "result_OE_m68_cont", "result_KD_m68_cont")
res_for_heatmap <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(res_for_heatmap) <- col_names
for (kinase in all_kinases) {
  for (col_name in col_names) {
    df=eval(parse(text = col_name))
    if (!kinase %in% rownames(df)) {
      nes_score=0
      fdr=1
    }
    else {
      nes_score=df[kinase, "ES"]
      fdr=df[kinase, "pval"]
    }
    col_name_for_plot=str_replace(col_name, "result_", "")
    if (nes_score < 0) {
      if (fdr==0) {
        res_for_heatmap[kinase, col_name] = log10(0.000101)
      }
      else {
        res_for_heatmap[kinase, col_name] = log10(fdr)
      }
    }
    else {
      if (fdr==0) {
        res_for_heatmap[kinase, col_name] = -log10(0.000101)
      }
      else {
        res_for_heatmap[kinase, col_name] = -log10(fdr)
      }
    }
  }
  sig_check=abs(res_for_heatmap[kinase, ])>2
  if (sum(sig_check)==0) {
    res_for_heatmap=res_for_heatmap[!(row.names(res_for_heatmap) %in% c(kinase)), ]
  }
}

#res_for_heatmap[res_for_heatmap==0] <- NA
res_for_heatmap<-res_for_heatmap[complete.cases(res_for_heatmap),]
labels_fdr <- c("OE_m25_cont", "KD_m25_cont", "OE_m68_cont", "KD_m68_cont")

colnames(res_for_heatmap) <- labels_fdr
#set.seed(123)
#fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
#  labs(subtitle = "Gap statistic method")

set.seed(123)
res_for <- res_for_heatmap
pdf(paste0("heatmap_CF1_m25_m65_OE_KD_filtered.pdf"), width = 10, height = 0.1969*90 + 1.3150)
draw(Heatmap(as.matrix(res_for), #column_split = rep(c("m25", "m68"), 2),
             row_names_gp = gpar(fontsize = 22), column_names_gp = gpar(fontsize = 22),
             km = 7, row_dend_width = unit(3, "cm"), column_dend_height = unit(3, "cm"), row_title_gp = gpar(fontsize = 24), column_title_gp = gpar(fontsize = 24), heatmap_legend_param = list(at=c(-4,-3,-2,-1,0,1,2,3,4), color = c("blue", "white", "red"), color_bar="continuous", direction = "horizontal", title = "log10(pval)", 
                                                                                                                                                                                                 legend_width = unit(10, "cm"), labels_gp = gpar(fontsize = 20),
                                                                                                                                                                                                 title_gp = gpar(fontsize = 20)), 
             height = unit(5, "mm")*60, width = unit(5, "mm")*15, 
             column_dend_reorder = FALSE), heatmap_legend_side = "top")
dev.off()

data_cmgc <- read.table("DE_GSEA_KSEA_txt_FILES/CMGC_kinase_family.txt", header = TRUE)
cmgc_kin=data_cmgc$Kinases

kkk=1
to_keep=list()
for (i in 1:length(cmgc_kin)) {
  ind <- which (kinases %in% cmgc_kin[i])
  if (length(ind)>5){
    to_keep[kkk]=i
    kkk=kkk+1
  }
}

all_cmgc_kin = cmgc_kin
to_keep_vector <- unlist(to_keep)
cmgc_kin=cmgc_kin[to_keep_vector]

col_names <- c("result_OE_m25_cont", "result_OE_m68_cont", "result_KD_m25_cont", "result_KD_m68_cont")
res_for_heatmap_cmgc <- data.frame(matrix(ncol = length(col_names), nrow = 0))
colnames(res_for_heatmap_cmgc) <- col_names
for (kinase in cmgc_kin) {
  for (col_name in col_names) {
    df=eval(parse(text = col_name))
    if (!kinase %in% rownames(df)) {
      nes_score=0
      fdr=1
    }
    else {
      nes_score=df[kinase, "ES"]
      fdr=df[kinase, "pval"]
    }
    col_name_for_plot=str_replace(col_name, "result_", "")
    if (nes_score < 0) {
      if (fdr==0) {
        res_for_heatmap_cmgc[kinase, col_name] = log10(0.000101)
      }
      else {
        res_for_heatmap_cmgc[kinase, col_name] = log10(fdr)
      }
      
    }
    else {
      if (fdr==0) {
        res_for_heatmap_cmgc[kinase, col_name] = -log10(0.000101)
      }
      else {
        res_for_heatmap_cmgc[kinase, col_name] = -log10(fdr)
      }
    }
  }
  sig_check=abs(res_for_heatmap_cmgc[kinase, ])>2
  if (sum(sig_check)==0) {
    res_for_heatmap_cmgc=res_for_heatmap_cmgc[!(row.names(res_for_heatmap_cmgc) %in% c(kinase)), ]
  }
}

res_for_heatmap_cmgc[res_for_heatmap_cmgc==0] <- NA
res_for_heatmap_cmgc<-res_for_heatmap_cmgc[complete.cases(res_for_heatmap_cmgc),]
labels_fdr <- c("m25_cont_OE", "m68_cont_OE", "m25_cont_KD", "m68_cont_KD")
colnames(res_for_heatmap_cmgc) <- labels_fdr
set.seed(123)
#fviz_nbclust(res_for_heatmap_cmgc, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
#  labs(subtitle = "Gap statistic method")

res_for_cmgc <- res_for_heatmap_cmgc
pdf(paste0("cmgc_heatmap_CF1_m25_m65_OE_KD.pdf"), width = 20, height = 0.1969*100 + 1.3150)
draw(Heatmap(as.matrix(res_for_cmgc), row_names_gp = gpar(fontsize = 22), column_names_gp = gpar(fontsize = 22),
             row_dend_width = unit(3, "cm"), column_dend_height = unit(3, "cm"), row_title_gp = gpar(fontsize = 24), column_title_gp = gpar(fontsize = 24), heatmap_legend_param = list(at=c(-4,-3,-2,-1,0,1,2,3,4), color = c("blue", "white", "red"), color_bar="continuous",direction = "horizontal", title = "log10(pval)", 
                                                                                                                                                                                        labels_gp = gpar(fontsize = 20),
                                                                                                                                                                                        title_gp = gpar(fontsize = 20)), 
             height = unit(10, "mm")*25, width = unit(12, "mm")*12, 
             column_dend_reorder = FALSE), heatmap_legend_side = "top")
dev.off()