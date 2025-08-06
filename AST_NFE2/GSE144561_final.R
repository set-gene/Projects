
gse144561= read.table("GSE144561_rawCountsAllsamples.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE)

library(tibble)

CTCs <- gse144561[,c(22:81)]

########################

library(dplyr)

library(EnsDb.Hsapiens.v86)


library(tibble)

edb <- EnsDb.Hsapiens.v86
genes_ensemble <- genes(edb)
gene_length<- as.data.frame(genes_ensemble) %>% dplyr::select(gene_id, gene_name, width)



library(dplyr)

library(EnsDb.Hsapiens.v86)


library(tibble)

edb <- EnsDb.Hsapiens.v86
genes_ensemble <- genes(edb)
gene_length<- as.data.frame(genes_ensemble) %>% dplyr::select(gene_id, gene_name, width)

all_clean<- CTCs[rowSums(CTCs) > 5, ]


gene_length_in_mat<- left_join(data.frame(gene_name = rownames(all_clean)), gene_length) %>% dplyr::filter(!is.na(width))

all_Matrix_sub<- all_clean[rownames(all_clean) %in% gene_length_in_mat$gene_name, ]

all.equal(rownames(all_Matrix_sub), gene_length_in_mat$gene_name)


sub_gene <- data.frame(gene_name = rownames(all_Matrix_sub))


gene_length_in_mat2 <- left_join(sub_gene ,gene_length_in_mat,by = "gene_name")

gene_length_in_mat2 <-gene_length_in_mat2[-which(duplicated(gene_length_in_mat2$gene_name)),]


###########################################################################


countToTpm <- function(counts, effLen)
{
  rate <- log(counts + 1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

all_TPM<- apply(all_Matrix_sub,2, countToTpm, effLen = gene_length_in_mat2$width)


tpm_log = log2(all_TPM+1)



#correlation-plot
ctc_all <- tpm_log

selects <- c("NFE2","ITGA2B")

cor_select <- ctc_all[selects,]

nn_cor <- data.frame(t(cor_select))


library(ggpubr)


library(envalysis)



scatter_cor <- ggscatter(nn_cor, x = "NFE2", y = "ITGA2B", 
                         xlab = "NFE2 expression level",
                         ylab = "ITGA2B expression level",
                         conf.int = TRUE,
                         ggtheme = theme_bw()# Add confidence interval
)+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4) # Add confidence interval



ggsave(filename = "final_result/CTC_NFE2_ITGA2B_COR.png", plot = scatter_cor , width = 12, height = 10, dpi = 300, units = "cm")



selects2 <- c("IKZF1","ITGA2B")

cor_select2<- ctc_all[selects2,]

nn_cor2 <- data.frame(t(cor_select2))


library(ggpubr)
scatter_cor2 <- ggscatter(nn_cor2, x = "IKZF1", y = "ITGA2B", 
                          xlab = "IKZF1 expression level",
                          ylab = "ITGA2B expression level",
                          conf.int = TRUE
                          ,ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval


ggsave(filename = "final_result/CTC_IKZF1_ITGA2B_COR.png", plot = scatter_cor2 , width = 12, height = 10, dpi = 300, units = "cm")


selects3 <- c("IRF8","ITGA2B")

cor_select3<- ctc_all[selects3,]

nn_cor3 <- data.frame(t(cor_select3))


scatter_cor3 <- ggscatter(nn_cor3, x = "IRF8", y = "ITGA2B", 
                          xlab = "IRF8 expression level",
                          ylab = "ITGA2B expression level",
                          conf.int = TRUE,
                          ggtheme = theme_bw() )+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidenceinterval


ggsave(filename = "final_result/CTC_IRF8_ITGA2B_COR.png", plot = scatter_cor3, width = 12, height = 10, dpi = 300, units = "cm")




selects4 <- c("BTG2","ITGA2B")

cor_select4<- ctc_all[selects4,]

nn_cor4 <- data.frame(t(cor_select4))


library(ggpubr)
scatter_cor4 <- ggscatter(nn_cor4, x = "BTG2", y = "ITGA2B", 
                          xlab = "BTG2 expression level",
                          ylab = "ITGA2B expression level",
                          conf.int = TRUE,
                          ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

ggsave(filename = "final_result/CTC_BTG2_ITGA2B_COR.png", plot = scatter_cor4, width = 12, height = 10, dpi = 300, units = "cm")


#############################################################################dendrogram_reclustering

library(dplyr)

set.seed(1234)

ctc_all

tpm_ctcs_t <- t(ctc_all )

tpm_ctcs_t <- as.data.frame(tpm_ctcs_t)

nom <- scale(tpm_ctcs_t)

d <- dist(nom, method = "maximum")

hc <- hclust(d, method = "ward.D2")

plot(hc,hang = -1,cex = 0.8)


library(factoextra)

ctc_group_dend <- fviz_dend(hc,k=2,lwd = 0.5,cex = 0.2, k_colors = c("#F8766D", '#619CFF'))+ theme_publish(base_size = 11, base_linewidth = 0.7)


ggsave(filename = "final_result/CTC_groups_dendrogram.png", plot = ctc_group_dend , width = 12, height = 10, dpi = 300, units = "cm")


library(factoextra)

clusters <- cutree(hc, k = 2)

factor(clusters)

library(tibble)

cluster_complete <- as.data.frame(nom)



clusters <- as.data.frame(clusters)

clusters <- rownames_to_column(clusters,"sample")

clusters2 <- clusters

clusters3 <- subset(x = clusters2,clusters2$clusters == 1)

clusters4 <- subset(x = clusters2,clusters2$clusters == 2)

cluster_complete$clusters <- clusters$clusters


#############################################re_groupping

tpm_ctcs <- ctc_all

tpm_1 <- tpm_ctcs[,clusters3$sample]

tpm_2 <- tpm_ctcs[,clusters4$sample]

tpm_re <- cbind(tpm_1 ,tpm_2)

tpm_re_t <- t(tpm_re)



tpm_re_df <- data.frame(tpm_re)

tpm_re_df <- rownames_to_column(tpm_re_df,"Gene")


ast = c("NFE2","HBA1","HBA2","HBB","GATA1","ITGA2B","GP1BA","SELP","SELPLG")


AST_exp_tpm2 = tpm_re[ast,]


AST_exp_tpm2 <- as.data.frame(t(AST_exp_tpm2))

AST_exp_tpm2$Group <- c(rep("CTC_1",38), rep("CTC_2",22))


AST_exp_tpm_GROUP <- data.frame(Group = AST_exp_tpm2$Group)

rownames(AST_exp_tpm_GROUP) <- rownames(AST_exp_tpm2)


library(tidyr)

library(tibble)

dat_tpm2 <- AST_exp_tpm2%>%
  rownames_to_column("Sample") %>% 
  gather(key = Genes,value = Proportion,-Sample,-Group)


library(ggplot2)

library(ggpubr)

library(rstatix)


dat_tpm2$Group <- factor(dat_tpm2$Group,levels = c("CTC_1","CTC_2"))


library(ggplot2)

library(ggpubr)

library(rstatix)
ctc_groups_violin <- ggplot(dat_tpm2 , aes(x=Genes, y=Proportion,fill = Group)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(aes(group = Group), method = "t_test", label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)  + ylab(" log2(TPM+1) expression level") + scale_fill_manual(values = c('CTC_1'='#F8766D','CTC_2'='#619CFF')) + theme_publish(base_size = 11, base_linewidth = 0.7)

ggsave(filename = "final_result/CTC_groups_violine.png", plot = ctc_groups_violin , width = 12, height = 10, dpi = 300, units = "cm")





tpm_re_t <- t(tpm_re)


Group <- factor(c( rep("CTC_1",38), rep("CTC_2",22)), levels = c("CTC_1","CTC_2"))


color <- c( rep("CTC_1",38), rep("CTC_2",22))

df_pca <- prcomp(tpm_re_t)

df_pcs2 <-data.frame(df_pca$x, Group = Group) 

percentage<-round(df_pca$sdev / sum(df_pca$sdev) * 100,2)

percentage<-paste(colnames(df_pcs2),"(", paste(as.character(percentage), "%", ")", sep=""))

library(ggplot2)

PCA_groups <- ggplot(df_pcs2,aes(x=PC1,y=PC2,color=Group ))+ geom_point(size = 1)+ geom_point() +theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line= element_line(colour = "black"))+ geom_point()+stat_ellipse()+ geom_hline(yintercept = c(0,0),linetype = "dashed",size = 1)+ geom_vline(xintercept = c(0,0),linetype = "dashed",size = 1) + scale_color_discrete(breaks=c("CTC_1","CTC_2"))+ scale_color_manual(values = c( "#F8766D", "#619CFF"))+ theme_publish(base_size = 11, base_linewidth = 0.7)


ggsave(filename = "final_result/CTC_groups_PCA.png", plot = PCA_groups , width = 12, height = 10, dpi = 300, units = "cm")


####################################################################################EMtype

all_TPM_ctc1 <- all_TPM[,clusters3$sample]

all_TPM_ctc2 <- all_TPM[,clusters4$sample]

all_tpm_re <- cbind(all_TPM_ctc1,all_TPM_ctc2)

all_tpm_re_df <- data.frame(all_tpm_re)
################
##using TPM counts
Epithelial_type <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                     "MUC1","MUC2")

Mesenchymal_type <- c("CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM")


Gene_sets = list(Epithelial_type,Mesenchymal_type)

names(Gene_sets) = c("Epithelial_type","Mesenchymal_type")



library(GSVA)

ssgsea = gsva(all_tpm_re , Gene_sets,method = "ssgsea")


Group <- factor(c(rep("CTC_1",38),rep("CTC_2",22)), levels = c("CTC_1","CTC_2"))

colData <- data.frame(row.names = colnames(all_tpm_re ),Group )


library(pheatmap)



ctc_em_heat <- pheatmap(ssgsea,
                        annotation_col = colData,
                        cluster_rows = F,
                        cluster_cols = T,
                        show_rownames = T,
                        show_colnames = F,
                        border_color = NA,
                        #fontsize = 11,
                        scale = "row",
                        #fontsize_row = 11,
                        cellwidth = 8, 
                        cellheight=80,
                        annotation_colors = list(Group=c(CTC_1 = "#F8766D",CTC_2="#619CFF")),
                        color=colorRampPalette(c("navy", "white", "red"))(50))

ggsave(filename = "final_result/CTC_groups_EM_heatmap.png", plot = ctc_em_heat, width = 30, height = 10, dpi = 300, units = "cm", background= "transparent")

#########################################################################################################################################

library(ggplot2)

library(xCell)

xcell_scores = xCellAnalysis(all_tpm_re)

xcell_select <- xcell_scores[c("Fibroblasts","Adipocytes","Platelets","Neutrophils","Macrophages"),]

ctc_celltype_heat <- pheatmap(xcell_select ,
                              annotation_col = colData,
                              annotation_row = rowData ,
                              cluster_rows = F,
                              cluster_cols = T,
                              show_rownames = T,
                              show_colnames = F,
                              border_color = NA,
                              cellwidth = 8, 
                              cellheight=80,
                              scale = "row",
                              annotation_colors = list(Group=c(CTC_1 = "#F8766D",CTC_2="#619CFF"),
                                                       Types = c(Stromal = "gold",Immune = "magenta")),
                              color=colorRampPalette(c("navy", "white", "red"))(50))



Types <- factor(c(rep("Stromal",2),rep("Immune",3)), levels = c("Stromal","Immune"))

rowData <- data.frame(row.names = rownames(xcell_select),Types  )

ggsave(filename = "final_result/CTC_groups_celltype_heatmap.png", plot = ctc_celltype_heat  , width = 30, height = 20, dpi = 300, units = "cm", background= "transparent")
##############################################################################################################

library(fgsea)

library(GSVA)


Platelet_activation <- gmtPathways("GOBP_REGULATION_OF_PLATELET_ACTIVATION.v2023.2.Hs.gmt")

Platelet_activation <- Platelet_activation$GOBP_REGULATION_OF_PLATELET_ACTIVATION

Platelet_aggregation <- gmtPathways("GOBP_PLATELET_AGGREGATION.v2023.2.Hs.gmt")

Platelet_aggregation <- Platelet_aggregation$GOBP_PLATELET_AGGREGATION


ROS_negative_response <-gmtPathways("gene sets for reactive oxygen.gmt")

ROS_negative_response <- ROS_negative_response $`Negative Regulation Of Response To Reactive Oxygen Species (GO:1901032)`

Anoikis_resistance <- gmtPathways("gene sets for anoikis.gmt")

Anoikis_resistance  <- Anoikis_resistance  $`Negative Regulation Of Anoikis (GO:2000811)`



Gene_sets3 = list(Platelet_activation,Platelet_aggregation,
                  ROS_negative_response,Anoikis_resistance )


names(Gene_sets3) = c("Platelet_activation","Platelet_aggregation",
                      "ROS_negative_response","Anoikis_resistance")


ssgsea3 = GSVA::gsva(all_tpm_re, Gene_sets3,method = "ssgsea")

library(pheatmap)

ctc_ssgsea_heat <- pheatmap(ssgsea3,
         annotation_col = colData,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         #fontsize = 11,
         scale = "row",
         #fontsize_row = 11,
         cellwidth = 8, 
         cellheight=80,
         annotation_colors = list(Group=c(CTC_1 = "#F8766D",CTC_2="#619CFF")),
         color=colorRampPalette(c("navy", "white", "red"))(50))

library(ggpubr)

ggsave(filename = "final_result/CTC_groups_ssgsea_heatmap.png", plot = ctc_ssgsea_heat  , width = 30, height = 20, dpi = 300, units = "cm",background= "transparent")

#####################################################################################################candidate genes for ROC

library(tibble)

sample_pheno = rownames_to_column(colData)

library(pROC)

tpm_re

ctc_tpm_log_re <- tpm_re

rld_RE = data.frame(t(ctc_tpm_log_re))

rld_RE2 = rownames_to_column(rld_RE)

rld_RE2= merge(sample_pheno,rld_RE2,by = "rowname")

rld_RE2 = rld_RE2[,c(-1)]



rld_RE2$Group = relevel(rld_RE2$Group, ref = "CTC_2")

library(pROC)
library(ggplot2)
library(tidyverse)

rld_RE2$Group
# example data
roc.list <- roc(Group ~ NFE2 + ITGA2B + GP1BA, data = rld_RE2)

roc.list


##########################################

rld_RE2_DF <-  data.frame(Attribute=c(colnames(rld_RE2)[2:39794]), AUC=NA)


for(i in 1:nrow(rld_RE2_DF)){
  roc_result <- roc(rld_RE2$Group, rld_RE2[,as.character(rld_RE2_DF$Attribute[i])])   # 확진 결과에 대한 데이터(type)와 진단 방법에 대한 후보 변수를 입력하여 AUC를 계산합니다. 
  rld_RE2_DF[i,'AUC'] <- roc_result$auc   # AUC 값을 입력합니다.
}

rld_RE2_DF <- rld_RE2_DF[order(-rld_RE2_DF$AUC),]

rownames(rld_RE2_DF) <- NULL


rld_RE2_DF2 <- column_to_rownames(rld_RE2_DF,"Attribute")

rld_RE2_SELECT <- data.frame(AUC = rld_RE2_DF2[c("IKZF1","IRF8","BTG2","NFE2","ITGA2B","GP1BA","HBA1","HBA2","HBB","GATA1","SELPLG"),])

rownames(rld_RE2_SELECT) <- c("IKZF1","IRF8","BTG2","NFE2","ITGA2B","GP1BA","HBA1","HBA2","HBB","GATA1","SELPLG")

rld_RE2_SELECT <- rownames_to_column(rld_RE2_SELECT,"Attribute")


rld_RE2_SELECT  <- rld_RE2_SELECT[order(-rld_RE2_SELECT$AUC),]

rownames(rld_RE2_SELECT) <- NULL

#####################################################

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(xtable)


library(envalysis)

roc_plot1 <- ggroc(roc.list,   
                   legacy.axes = TRUE) +  scale_colour_manual(values=c("#F8766D", "#619CFF", "#00BA38"))  + xlab("1-specificity") + ylab("Sensitivity") +theme_bw() + guides(color = guide_legend(title = "Genes"),size = 2) + annotate("text", x=0.22, y=0.95, label="AUC = 0.95",size=3) + annotate("text", x=0.26, y=0.84, label="AUC = 0.90",size=3) +  annotate("text", x=0.3, y=0.63, label="AUC = 0.75",size=3) + theme_publish(base_size = 11, base_linewidth = 0.7)

library(ggpubr)

roc_table_re <- rld_RE2_SELECT

colnames(roc_table_re)[2] <- "AUC 
  (Area Under the ROC Curve)"


roc_table <- ggtexttable(roc_table_re, rows = NULL, theme = ttheme("light"))

roc_table 

roc_group <- ggarrange(roc_plot1, roc_table,
                       ncol = 2, nrow = 1,
                       heights = c(1, 1.5),
                       widths = c(1.5,1),vjust = 10)



ggsave(filename = "final_result/CTC_groups_roc_plot.png", plot =roc_group  , width = 12, height = 10, dpi = 300, units = "cm")

############################################################






library(dplyr)

set.seed(7777)

platelet <- c("ITGA2B","GP1BA")


ctc_plt_counts<- ctc_all[platelet,]

tpm_ctc_plt_t <- t(ctc_plt_counts)

tpm_ctc_plt_t  <- as.data.frame(tpm_ctc_plt_t )

nom2 <- scale(tpm_ctc_plt_t )

d2 <- dist(nom2, method = "maximum")

hc2 <- hclust(d2, method = "ward.D2")

plot(hc2,hang = -1,cex = 0.8)

clusterss <- cutree(hc2, k = 2)

library(tibble)

cluster_complete2 <- as.data.frame(nom2)



clusterss2 <- as.data.frame(clusterss)

clusterss2 <- rownames_to_column(clusterss2,"sample")

clusterss3 <- clusterss2

plt <- subset(x = clusterss3 ,clusterss3$clusterss == 1)

ctc <- subset(x = clusterss3,clusterss3$clusterss == 2)


############################################re_groupping for platelet marker

tpm_plt <- tpm_ctcs[,plt$sample]

tpm_ctc <- tpm_ctcs[,ctc$sample]


tpm_re_plt <- cbind(tpm_plt,tpm_ctc)

Group <- factor(c(rep("CTC_plt",32),rep("CTC",28)), levels = c("CTC_plt","CTC"))

colData_plt <- data.frame(row.names = colnames(tpm_re_plt ),Group )



group_plt_heatmap <- pheatmap(ctc_plt_counts,
                              cluster_cols = hc2
                              ,scale = "row", 
                              annotation_col = colData_plt,
                              cellwidth = 8, 
                              cellheight=80,
                              cutree_cols = 2,
                              annotation_colors = list(Group=c(CTC_plt = "#F8766D",CTC="#619CFF")),
                              color=colorRampPalette(c("navy", "white", "red"))(50))


ggsave(filename = "final_result/plt/CTC_plt_heatmap.png", plot = group_plt_heatmap, width = 30, height = 20, dpi = 300, units = "cm",background= "transparent")


AST_exp_plt_tpm = tpm_re_plt[ast,]

#AST_exp_tpm2 = AST_exp_tpm2[-9,]

AST_exp_plt_tpm  <- as.data.frame(t(AST_exp_plt_tpm ))

AST_exp_plt_tpm$Group <- c( rep("CTC_plt",32), rep("CTC",28))

AST_exp_plt_tpm_GROUP <- data.frame(Group = AST_exp_plt_tpm$Group)

rownames(AST_exp_plt_tpm_GROUP) <- rownames(AST_exp_plt_tpm)

library(tidyr)

library(tibble)

dat_plt_tpm <- AST_exp_plt_tpm%>%
  rownames_to_column("Sample") %>% 
  gather(key = Genes,value = Proportion,-Sample,-Group)


library(ggplot2)

library(ggpubr)

library(rstatix)


dat_plt_tpm$Group<- factor(dat_plt_tpm$Group,levels = c("CTC_plt","CTC"))
###############################################

ctc_plt_violin <- ggplot(dat_plt_tpm  , aes(x=Genes, y=Proportion,fill = Group)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(aes(group = Group), method = "t_test", label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)  + ylab(" log2(TPM+1) expression level") + scale_fill_manual(values = c('CTC_plt'='#F8766D','CTC'="#619CFF")) + theme_publish(base_size = 11, base_linewidth = 0.7)



ggsave(filename = "final_result/plt/CTC_PDAC_plt_violine.png", plot = ctc_plt_violin , width = 12, height = 10, dpi = 300, units = "cm")


########################################################ctc_plt


all_TPM_plt <- all_TPM[,plt$sample]

all_TPM_ctc <- all_TPM[,ctc$sample]

all_tpm_plt_re <- cbind(all_TPM_plt,all_TPM_ctc)

#############################################################

library(GSVA)

ssgsea4 = gsva(all_tpm_plt_re, Gene_sets,method = "ssgsea")

colData_plt


library(pheatmap)

em_plt_ssgsea <- pheatmap(ssgsea4,
                          annotation_col = colData_plt,
                          cluster_rows = F,
                          cluster_cols = T,
                          show_rownames = T,
                          show_colnames = F,
                          border_color = NA,
                          scale = "row",
                          cellwidth = 8, 
                          cellheight=80,
                          annotation_colors = list(Group=c(CTC_plt = "#F8766D",CTC="#619CFF")),
                          color=colorRampPalette(c("navy", "white", "red"))(50))


ggsave(filename = "final_result/plt/CTC_plt_em_heatmap.png", plot = em_plt_ssgsea  , width = 30, height = 10, dpi = 300, units = "cm",background= "transparent")

###########################################################3

ssgsea5 = GSVA::gsva(all_tpm_plt_re , Gene_sets3,method = "ssgsea")


library(pheatmap)

ssgsea_pathway_plt <- pheatmap(ssgsea5,
                               annotation_col = colData_plt,
                               cluster_rows = F,
                               cluster_cols = T,
                               show_rownames = T,
                               show_colnames = F,
                               border_color = NA,
                               scale = "row",
                               cellwidth = 8, 
                               cellheight=80,
                               annotation_colors = list(Group=c(CTC_plt = "#F8766D",CTC="#619CFF")),
                               color=colorRampPalette(c("navy", "white", "red"))(50))


ggsave(filename = "final_result/plt/CTC_plt_ssgsea_heatmap.png", plot = ssgsea_pathway_plt  , width = 30, height = 20, dpi = 300, units = "cm",background= "transparent")


#################candidate genes for ROC

library(tibble)

sample_pheno_plt = rownames_to_column(colData_plt)


ctc_log_tpm_re_plt <-tpm_re_plt
library(pROC)

rld_RE_plt = data.frame(t(ctc_log_tpm_re_plt))

rld_RE2_plt = rownames_to_column(rld_RE_plt)

rld_RE2_plt= merge(sample_pheno_plt,rld_RE2_plt,by = "rowname")

rld_RE2_plt = rld_RE2_plt[,c(-1)]


rld_RE2_plt$Group = relevel(rld_RE2_plt$Group, ref = "CTC")


roc.plt.list <- roc(Group ~ NFE2 + ITGA2B + GP1BA, data = rld_RE2_plt)

roc.plt.list 
