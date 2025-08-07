library(openxlsx)

mycounts <- read.xlsx("READ.xlsx")

genes <- read.xlsx("GSE60407_platform.xlsx")

genes <- genes[,c(1,4)]

mycounts2 <- merge(genes,mycounts,"ID")

mycounts2 <- mycounts2[,c(-1)]


mycounts3= aggregate(mycounts2[, -c(1)],
                      by = list(Gene = mycounts2$symbol),
                      FUN = mean,
                      na.rm = TRUE)

mycounts3 <- column_to_rownames(mycounts3,"Gene")

mycounts3 <- round(mycounts3)

#############################################################################################################TPM

library(dplyr)

library(EnsDb.Hsapiens.v75)

library(tibble)

edb <- EnsDb.Hsapiens.v75
genes_ensemble <- genes(edb)
gene_length<- as.data.frame(genes_ensemble) %>% dplyr::select(gene_id, gene_name, width)


all_clean<- mycounts3

gene_length_in_mat<- left_join(data.frame(gene_name = rownames(all_clean)), gene_length) %>% dplyr::filter(!is.na(width))

all_Matrix_sub<- all_clean[rownames(all_clean) %in% gene_length_in_mat$gene_name, ]

all.equal(rownames(all_Matrix_sub), gene_length_in_mat$gene_name)


sub_gene <- data.frame(gene_name = rownames(all_Matrix_sub))


gene_length_in_mat2 <- left_join(sub_gene ,gene_length_in_mat,by = "gene_name")

gene_length_in_mat2 <-gene_length_in_mat2[-which(duplicated(gene_length_in_mat2$gene_name)),]


countToTpm <- function(counts, effLen)
{
  rate <- log(counts + 1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

all_TPM<- apply(all_Matrix_sub,2, countToTpm, effLen = gene_length_in_mat2$width)


tpm_log = log2(all_TPM+1)


rlog_all_types <- tpm_log[c("IKZF1","IRF8","NFE2","BTG2","HBA1","HBA2","GATA1","SELPLG","JUP","DSP","ITGA2B","GP1BA","SELP", "CD14","CD68", "ITGAM","PTPRC","CSF3R"),]

tpm_df <- data.frame(tpm_log)

tpm_df <- rownames_to_column(tpm_df,"Gene")

library(pheatmap)

rlog_all_types2 <- rlog_all_types[,-c(6,7)]

pheatmap(rlog_all_types ,
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

####################################################################################################
Epithelial_type <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                     "MUC1","MUC2")

Mesenchymal_type <- c("CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM")


Gene_sets = list(Epithelial_type,Mesenchymal_type)

names(Gene_sets) = c("Epithelial_type","Mesenchymal_type")


library(GSVA)

ssgsea = gsva(tpm_log , Gene_sets,method = "ssgsea")


library(pheatmap)

pheatmap(ssgsea,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         cellwidth = 8, 
         cellheight=80,
         labels_row = as.expression(newnames),
         color=colorRampPalette(c("navy", "white", "red"))(50))


newnames <- lapply(
  rownames(ssgsea),
  function(x) bquote(bold(.(x))))


#######################################################################


tpm_ctc_plt <- tpm_log[c("ITGA2B","ITGB3","GP1BA","SELP"),]


set.seed(1234)

tpm_ctc_plt_t <- t(tpm_ctc_plt)


nom <- scale(tpm_ctc_plt_t)


d <- dist(nom, method = "maximum")

hc <- hclust(d, method = "ward.D")


plot(hc,hang = -1,cex = 0.8)


clusters <- cutree(hc, k = 2)


library(gplots)
library(dplyr)

heatmap.2(tpm_ctc_plt,col = colorRampPalette(c("navy", "white", "red")), Colv=as.dendrogram(hc), scale="row", density.info="none", trace="none",cexRow = 0.75,cexCol = 0.75) 

library(tibble)

clusters <- as.data.frame(clusters)

clusters <- rownames_to_column(clusters,"sample")


clusters2 <- clusters

clusters2 

PLT_low <- subset(x = clusters2,clusters2$clusters == 1)

PLT_high <- subset(x = clusters2,clusters2$clusters == 2)



PLT_low_1 <- tpm_log[,PLT_low$sample]

PLT_high_2 <- tpm_log[,PLT_high$sample]


re_ctc_plt <- cbind(PLT_low_1,PLT_high_2)


Group <- factor(c(rep("PLT_low",7), rep("PLT_high",8)), levels = c("PLT_low","PLT_high"))


colData2 <- data.frame(row.names = colnames(re_ctc_plt),Group )


Epithelial_type <- c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                     "MUC1","MUC2")

Mesenchymal_type <- c("CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM")



tpm_ctc_plt_EM <- re_ctc_plt[c("CD276","CEACAM5","CDH1","KRT14","KRT19","KRT7",
                               "MUC1","MUC2",
                               "CDH2","SNAI2","SNAI1","SPARC","TWIST1","VIM" ),]


tpm_ctc_plt_M <- re_ctc_plt[c("VIM","NFE2"),]


tpm_ctc_plt_NFE2 <- re_ctc_plt[c("NFE2","GATA1","VIM","CDH1"),]


library(pheatmap)

pheatmap(tpm_ctc_plt_NFE2 ,
         annotation_col = colData2,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         color=colorRampPalette(c("navy", "white", "red"))(50))
#########################################################################

library(xCell)

xcell_scores = xCellAnalysis(re_ctc_plt)

xcell_select <- xcell_scores[c("Fibroblasts","Adipocytes","Platelets","Neutrophils","Macrophages","Monocytes"),]

library(pheatmap)

pheatmap(xcell_select ,
         annotation_col = colData2,
         annotation_row = rowData ,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         annotation_colors = list(Group=c(PLT_low = "#F8766D",PLT_high="#00BA38"),
                                  Types = c(Stromal = "gold",Immune = "magenta")),
         color=colorRampPalette(c("navy", "white", "red"))(50))

newnames2 <- lapply(
  rownames(xcell_select),
  function(x) bquote(bold(.(x))))


Types <- factor(c(rep("Stromal",2),rep("Immune",4)), levels = c("Stromal","Immune"))

rowData <- data.frame(row.names = rownames(xcell_select),Types  )
