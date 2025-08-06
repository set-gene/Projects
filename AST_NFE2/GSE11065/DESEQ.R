setwd("C:/Users/jeeh9/Documents/sm/AST/GSE111065/hs_single")

filenames <- list.files(path = getwd())

numfiles <- length(filenames)

ctc_single <- do.call(cbind,lapply(filenames,read.table))

names(ctc_single) <- ctc_single[1,]

ctc_single <- ctc_single[-1,]

library(dplyr)
library(tibble)

rownames(ctc_single) <- NULL

ctc_single <- column_to_rownames(ctc_single,"Geneid")

ctc_single <- ctc_single[ , !colnames(ctc_single)%in%c("Geneid")]

setwd("C:/Users/jeeh9/Documents/sm/AST/GSE111065/hs_cluster")

filenames <- list.files(path = getwd())

numfiles <- length(filenames)

ctc_cluster <- do.call(cbind,lapply(filenames,read.table))

names(ctc_cluster) <- ctc_cluster[1,]

ctc_cluster <- ctc_cluster[-1,]

library(dplyr)
library(tibble)

rownames(ctc_cluster) <- NULL

ctc_cluster <- column_to_rownames(ctc_cluster,"Geneid")

ctc_cluster <- ctc_cluster[ , !colnames(ctc_cluster)%in%c("Geneid")]

ctc_all <- cbind(ctc_single,ctc_cluster)


chars <- sapply(ctc_all, is.character)

#convert all character columns to numeric
ctc_all[ , chars] <- as.data.frame(apply(ctc_all[ , chars], 2, as.numeric))



Group <- factor(c(rep("SC",48),rep("CC",21)), levels = c("SC","CC"))


colData <- data.frame(row.names = colnames(ctc_all),Group)


colData$Group = relevel(colData$Group, ref = "SC")

library(DESeq2)

ctc_all <- as.matrix(ctc_all)

dds <-DESeqDataSetFromMatrix(ctc_all,colData,design = ~ Group)


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq2::DESeq(dds)


resultsNames(dds)

res = DESeq2::results(dds ,c("Group", "CC","SC"))

res2 = res[order(res$padj),]

res2


cut_pvalue <- 0.05

cut_lfc <- 0.5

significant_results <- res2[which(res2$padj < cut_pvalue),]

significant_results

significant_results_EXP <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_EXP = data.frame(significant_results_EXP)

significant_results_High <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange>cut_lfc)),]

significant_results_High <- rownames_to_column(as.data.frame(significant_results_High), "Gene")

library(openxlsx)

write.xlsx(significant_results_High,"sig_H.xlsx")


significant_results_Low<- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) )),]

significant_results_Low <- rownames_to_column(as.data.frame(significant_results_Low), "Gene")

res2_df <- data.frame(res2)

res2_df <- rownames_to_column(res2_df,"Gene")


rld <- DESeq2::vst(dds, blind = FALSE)

rlogMat <- assay(rld)
#################################################################################
pheatmap(rlogMat[c("IKZF1","IRF8","NFE2","BTG2",
                   "HBA1","HBA2"),],
         show_colnames = T,
         show_rownames = T,
         cluster_cols = F,
         annotation_col = colData,
         fontsize = 8,
         scale = "row", 
         color=colorRampPalette(c("navy", "white", "red"))(50))



#################################################################################
rlog_df <- as.data.frame(rlogMat)

rlog_df <- rownames_to_column(rlog_df,"Gene")

Epithelial_type <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                     "MUC1","MUC2")

Mesenchymal_type <- c("CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM")


Gene_sets = list(Epithelial_type,Mesenchymal_type)

names(Gene_sets) = c("Epithelial_type","Mesenchymal_type")

rlog_cc <- rlogMat[,49:69]

library(GSVA)

ssgsea = gsva(rlogMat, Gene_sets,method = "ssgsea")

ssgsea_cc <- ssgsea[,49:69]

library(pheatmap)

pheatmap(ssgsea,
         show_colnames = F,
         annotation_col = colData,
         cluster_rows = T,
         cluster_cols = T,
         fontsize = 8,
         scale = "row", 
         color=colorRampPalette(c("navy", "white", "red"))(50))

###########################################################################

ssgsea_df = as.data.frame(ssgsea)

ssgsea_df = as.data.frame(t(ssgsea_df))

library(tibble)

ssgsea_df = rownames_to_column(ssgsea_df,"Sample")

colData_df <- rownames_to_column(colData,"Sample")

ssgsea_df = merge(ssgsea_df,colData_df,by = "Sample")

ssgsea_df = column_to_rownames(ssgsea_df,"Sample")

library(ggplot2)

library(ggpubr)

ggboxplot(ssgsea_df,x = "Group",
          y = "Mesenchymal_type", color = "Group", add = "jitter",ylab = "Ferroptosis") + 
  stat_compare_means(label = "p.signif", method = "t.test") 

###########################################################################3


ALL_rlog_EM <- rlogMat[c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                        "MUC1","MUC2","CDH2","SNAI1","SPARC","TWIST1","VIM"),]



pheatmap(ALL_rlog_EM ,
         annotation_col = colData ,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         color=colorRampPalette(c("navy", "white", "red"))(50))




##############################################################

cc_mesenchymal <- data.frame(Mesenchymal_type = ssgsea_cc["Mesenchymal_type",])


cc_nfe <- data.frame(NFE2 = rlog_cc["NFE2",])


cc_cor <- cbind(cc_nfe,cc_mesenchymal)

library(ggpubr)


scatter_cor_ctc <- ggscatter(cc_cor, x = "NFE2", y = "Mesenchymal_type", xlab = "NFE2 expression level in CTC_cluster",
                             ylab = "Mesenchymal_type in CtC_cluster",
                             conf.int = TRUE,size = 5,
                             ggtheme = theme_bw())# Add confidence interval


scatter_cor_ctc + stat_cor(method = "pearson",label.sep = "\n", size = 8,label.y = 0.5)+  theme(text = element_text(face = "bold",size=20))

###############################################################

cc_rlog_EM <- rlog_cc[c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                        "MUC1","MUC2","CDH2","SNAI1","SPARC","TWIST1","VIM"),]



cols.cor2 <- cor(cc_rlog_EM, use = "pairwise.complete.obs", method = "pearson")

rows.cor <- cor(t(cc_rlog_EM), use = "pairwise.complete.obs", method = "pearson")

library(pheatmap)

pheatmap(cc_rlog_EM ,
         annotation_row = rowData ,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         clustering_distance_cols = as.dist(1 - cols.cor2),
         clustering_distance_rows = as.dist(1 - rows.cor),
         labels_row = as.expression(newnames3),
         annotation_colors = list(Types = c(Epithelial_type = "gold",Mesenchymal_type = "magenta")),
         color=colorRampPalette(c("navy", "white", "red"))(50))



cc_rlog_M <- rlog_cc[c("CDH2","SNAI1","SPARC","TWIST1","VIM"),]


cols.cor3 <- cor(cc_rlog_M, use = "pairwise.complete.obs", method = "pearson")

rows.cor2 <- cor(t(cc_rlog_M), use = "pairwise.complete.obs", method = "pearson")



pheatmap(cc_rlog_EM ,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         color=colorRampPalette(c("navy", "white", "red"))(50))





newnames3 <- lapply(
  rownames(cc_rlog_EM),
  function(x) bquote(bold(.(x))))


Types <- factor(c(rep("Epithelial_type",8),rep("Mesenchymal_type",5)), levels = c("Epithelial_type","Mesenchymal_type"))

rowData <- data.frame(row.names = rownames(cc_rlog_EM),Types  )


#################################################################

################################################################


rlog_EM <- rlogMat[c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                        "MUC1","MUC2","CDH2","SNAI1","SPARC","TWIST1","VIM"),]





cols.cor3 <- cor(rlog_EM, use = "pairwise.complete.obs", method = "pearson")

rows.cor2 <- cor(t(rlog_EM), use = "pairwise.complete.obs", method = "pearson")

library(pheatmap)

pheatmap(rlog_EM ,
         annotation_col = colData,
         annotation_row = rowData ,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         annotation_colors = list(Types = c(Epithelial_type = "gold",Mesenchymal_type = "magenta")),
         color=colorRampPalette(c("navy", "white", "red"))(50))




newnames4 <- lapply(
  rownames(rlog_EM),
  function(x) bquote(bold(.(x))))




Types <- factor(c(rep("Epithelial_type",8),rep("Mesenchymal_type",5)), levels = c("Epithelial_type","Mesenchymal_type"))

rowData2 <- data.frame(row.names = rownames(rlog_EM),Types  )

############################################################





#################################################################
library(ggplot2)

library(dplyr)

library(tibble)

CC_SC <- res2_df %>%
  mutate(Group = case_when(log2FoldChange >  1 & padj < 0.05 ~ "CC",
                           log2FoldChange< -1 & padj < 0.05 ~ "SC",
                           TRUE ~ "Stable"))   

CC_SC 

CC_SC$Group <- factor(CC_SC$Group,levels = c("SC","Stable","CC"))

significnat_exp_genes <- rownames(significant_results)

CC_SC <- column_to_rownames(CC_SC,"Gene")

options(ggrepel.max.overlaps=Inf)

library(ggplot2)

g

volcano_all <- ggplot(data = CC_SC,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = Group), 
             alpha = 0.2, 
             shape = 16,
             size =3) +
  scale_color_manual(values = c("SC" = "#F8766D","Stable" = "grey","CC" = "#00BA38")) +geom_point(data = CC_SC["NFE2", ],
             shape = 21,
             size = 3, 
             fill = "purple", 
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             size = 1) +
  geom_vline(xintercept = c(-1,0, 1),
             linetype = "dashed",
             size = 1) +
  geom_label_repel(data = CC_SC[significnat_exp_genes,],     
                   aes(label = significnat_exp_genes),
                   force = 1,
                   nudge_y = 1,
                   size = 3) +  theme(text = element_text(size=11)) + ylab("-log10(adj.P_value)") + xlab("Log2FoldChange")+ guides(colour = guide_legend(override.aes = list(size=3))) + theme_bw() + theme_publish(base_size = 11, base_linewidth = 0.7)

library(envalysis)

ggsave(filename = "CTC_CC_single_volcano.png", plot = volcano_all , width = 15, height = 14, dpi = 300, units = "cm")

ggplot(data = CC_SC,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  geom_point(aes(colour = Group), 
             alpha = 0.2, 
             shape = 16,
             size =5) +
  scale_color_manual(values = c("SC" = "#F8766D","Stable" = "grey","CC" = "#00BA38")) +
  geom_point(data = CC_SC["NFE2", ],
             shape = 21,
             size = 5, 
             fill = "purple", 
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             size = 1) +
  geom_vline(xintercept = c(-1,0, 1),
             linetype = "dashed",
             size = 1) +
  geom_label_repel(data = CC_SC[significnat_exp_genes,],     
                   aes(label = significnat_exp_genes),
                   force = 1,
                   nudge_y = 1,
                   size = 3) +  theme(text = element_text(size=11)) + ylab("-log10(adj.P_value)") + xlab("Log2FoldChange")+ guides(colour = guide_legend(override.aes = list(size=5))) + theme_bw() + theme_publish(base_size = 11, base_linewidth = 0.7)

significnat_exp_genes

library(ggrepel)
##########################


library(tibble)

library(dplyr)

library(fgsea)

res22 = rownames_to_column(data.frame(res2), "SYMBOL")


res3 <- res22 %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))


ranks <- deframe(res3)



canonical = gmtPathways("c2.cp.kegg.v2023.1.Hs.symbols (2).gmt")

canonical$regulation

canonical %>% 
  head() %>% 
  lapply(head)

canonical_Res <- fgsea(pathways=canonical, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)



GO = gmtPathways("c5.go.bp.v2023.1.Hs.symbols (1).gmt")

GO_Res <- fgsea(pathways=GO, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)


go_ResT <- GO_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

go_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

go_ResT_up = subset(x = go_ResT, go_ResT$NES > 1)

go_ResT_up <- go_ResT_up[order(go_ResT_up$padj),]

go_ResT_up_sig <- subset(go_ResT_up,go_ResT_up$padj < 0.05)

library(openxlsx)

write.xlsx(go_ResT_up_sig,"GO_BP_up_significant_GSEA.xlsx")

go_ResT_dn = subset(x = go_ResT, go_ResT$NES < -1)


go_ResT_dn <- go_ResT_dn[order(go_ResT_dn$padj),]


go_ResT_dn_sig <- subset(go_ResT_dn,go_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(go_ResT_dn_sig,"GO_BP_dn_significant_GSEA.xlsx")

#################################################################

library(xCell)
xCell.data$genes


xcell_scores = xCellAnalysis(rlogMat)

xcell_select <- xcell_scores[c("Fibroblasts","Adipocytes","Platelets","Neutrophils","Macrophages"),]

pheatmap(xcell_select ,
         annotation_col = colData,
         annotation_row = rowData ,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         fontsize = 10,
         cellwidth = 8, 
         cellheight=80,
         scale = "row",
         fontsize_row = 10,
         annotation_colors = list(Group=c(SC = "#F8766D",CC="#00BA38"),
                                  Types = c(Stromal = "gold",Immune = "magenta")),
         color=colorRampPalette(c("navy", "white", "red"))(50))


newnames2 <- lapply(
  rownames(xcell_select),
  function(x) bquote(bold(.(x))))


Types <- factor(c(rep("Stromal",2),rep("Immune",3)), levels = c("Stromal","Immune"))

rowData <- data.frame(row.names = rownames(xcell_select),Types  )
#######################################################################################

###################################################################################FOR_GSEA_Prelanked
for_gsea_pre = data.frame(res2)

for_gsea_pre = rownames_to_column(for_gsea_pre,"Gene")

for_gsea_pre$fcsign = sign(for_gsea_pre$log2FoldChange)

for_gsea_pre$logP = -log10(for_gsea_pre$pvalue)

for_gsea_pre$metric = for_gsea_pre$logP/for_gsea_pre$fcsign

for_gsea = for_gsea_pre[,c("Gene","metric")]

for_gsea = na.omit(for_gsea)

#write.table(for_gsea,file="DE_sig_genes.rnk",quote=F,sep="\t",row.names=F)


write.table(for_gsea,file="DE_sig_genes_cc_sc.rnk",quote=F,sep="\t",row.names=F,col.names = F)



######################################################################################


library(clusterProfiler)
library(enrichplot)
library(ggplot2)

library(tibble)

res2


significant_results_up<- res2[which(res2$padj<cut_pvalue & res2$log2FoldChange>cut_lfc),]

significant_results_up

significant_results_sorted_up<- significant_results_up[order(significant_results_up$padj),]

significant_results_down<- res2[which(res2$padj<cut_pvalue & res2$log2FoldChange<(-cut_lfc)),]

significant_results_sorted_down<- significant_results_down[order(significant_results_down$padj),]


#Input genes; convert to ENTREZID 
eg_u= bitr(rownames(significant_results_sorted_up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

eg_d = bitr(rownames(significant_results_sorted_down), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

write.xlsx(eg_u,"up_sig_genes.xlsx")

go_u_bp <- enrichGO(eg_u$SYMBOL, OrgDb = "org.Hs.eg.db", ont="BP",keyType = "SYMBOL")

go_d_bp <- enrichGO(eg_d$SYMBOL, OrgDb = "org.Hs.eg.db", ont="BP",keyType = "SYMBOL")


go_u_bp_df <- as.data.frame(go_u_bp@result)

go_u_bp_df <- subset(go_u_bp_df,go_u_bp_df$p.adjust < 0.05)

################################################################

Epithelial_type <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                     "MUC1","MUC2")

Mesenchymal_type <- c("CDH2","SNAI1","SPARC","TWIST1","VIM")

Desmosome <- c("JUP","DSP")

Platelet <- c("ITGA2B","GP1BA","SELP")

Macrophage <- c("CD14","CD68")

Neutrophil <- c("ITGAM","PTPRC","CSF3R")

all_types <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
               "MUC1","MUC2","CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM",
               "JUP","DSP","ITGA2B","GP1BA","SELP", "CD14","CD68",
               "ITGAM","PTPRC","CSF3R")


rlog_all_types <- rlogMat[c("JUP","DSP","ITGA2B","GP1BA","SELP", "CD14","CD68",
                            "ITGAM","PTPRC","CSF3R"),]



rlog_cc_types <- rlog_cc[c("JUP","DSP","ITGA2B","GP1BA","SELP", "CD14","CD68",
                            "ITGAM","PTPRC","CSF3R"),]

library(pheatmap)

pheatmap(rlog_all_types ,
         annotation_col = colData,
         annotation_row = rowData ,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         annotation_colors = list(Types = c(Epithelial_type = "gold",Mesenchymal_type = "magenta")),
         color=colorRampPalette(c("navy", "white", "red"))(50))



pheatmap(rlog_all_types ,
         annotation_col = colData,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         color=colorRampPalette(c("navy", "white", "red"))(50))


newnames5 <- lapply(
  rownames(rlog_EM),
  function(x) bquote(bold(.(x))))




pheatmap(rlog_cc_types ,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 10, 
         cellheight=10,
         color=colorRampPalette(c("navy", "white", "red"))(50))


Types <- factor(c(rep("Epithelial_type",8),rep("Mesenchymal_type",5),rep("Desmosome",2)), levels = c("Epithelial_type","Mesenchymal_type"))

rowData2 <- data.frame(row.names = rownames(rlog_EM),Types  )
