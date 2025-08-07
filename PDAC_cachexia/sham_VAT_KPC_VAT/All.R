sham_VAT_1 <- read.table("SRR14136106.txt", sep = "\t",header = T,quote = "")

sham_VAT_2 <- read.table("SRR14136107.txt", sep = "\t",header = T,quote = "")

sham_VAT_3 <-  read.table("SRR14136108.txt", sep = "\t",header = T,quote = "")

sham_VAT_4 <-  read.table("SRR14136109.txt", sep = "\t",header = T,quote = "")

sham_VAT_5 <-  read.table("SRR14136110.txt", sep = "\t",header = T,quote = "")

KPC_VAT_1 <- read.table("SRR14136111.txt", sep = "\t",header = T,quote = "")

KPC_VAT_2 <- read.table("SRR14136112.txt", sep = "\t",header = T,quote = "")

KPC_VAT_3 <-  read.table("SRR14136113.txt", sep = "\t",header = T,quote = "")

KPC_VAT_4 <-  read.table("SRR14136114.txt", sep = "\t",header = T,quote = "")

KPC_VAT_5 <-  read.table("SRR14136115.txt", sep = "\t",header = T,quote = "")


colnames(sham_VAT_1)[7]<-("sham_VAT_1")
colnames(sham_VAT_2)[7]<-("sham_VAT_2")
colnames(sham_VAT_3)[7]<-("sham_VAT_3")
colnames(sham_VAT_4)[7]<-("sham_VAT_4")
colnames(sham_VAT_5)[7]<-("sham_VAT_5")
colnames(KPC_VAT_1)[7]<-("KPC_VAT_1")
colnames(KPC_VAT_2)[7]<-("KPC_VAT_2")
colnames(KPC_VAT_3)[7]<-("KPC_VAT_3")
colnames(KPC_VAT_4)[7]<-("KPC_VAT_4")
colnames(KPC_VAT_5)[7]<-("KPC_VAT_5")

sham_VAT_1 <- sham_VAT_1[,c(1,7)]
sham_VAT_2 <- sham_VAT_2[,c(1,7)]
sham_VAT_3 <- sham_VAT_3[,c(1,7)]
sham_VAT_4 <- sham_VAT_4[,c(1,7)]
sham_VAT_5 <- sham_VAT_5[,c(1,7)]
KPC_VAT_1 <- KPC_VAT_1[,c(1,7)]
KPC_VAT_2 <- KPC_VAT_2[,c(1,7)]
KPC_VAT_3 <- KPC_VAT_3[,c(1,7)]
KPC_VAT_4 <- KPC_VAT_4[,c(1,7)]
KPC_VAT_5 <- KPC_VAT_5[,c(1,7)]

library(dplyr)

n_raw_count_1 <- merge(sham_VAT_1,sham_VAT_2,by="Geneid") %>%
  merge(.,sham_VAT_3,by = "Geneid") %>%
  merge(.,sham_VAT_4,by = "Geneid") %>%
  merge(.,sham_VAT_5,by = "Geneid") %>%
  merge(.,KPC_VAT_1,by = "Geneid") %>%
  merge(.,KPC_VAT_2,by = "Geneid") %>%
  merge(.,KPC_VAT_3,by = "Geneid") %>%
  merge(.,KPC_VAT_4,by = "Geneid") %>%
  merge(.,KPC_VAT_5,by = "Geneid") 
 
library(tibble)

n_raw_count_1 <- column_to_rownames(n_raw_count_1,"Geneid")


sham_VAT_6 <- read.table("SRR8278867.txt", sep = "\t",header = T,quote = "")

sham_VAT_7 <- read.table("SRR8278868.txt", sep = "\t",header = T,quote = "")

sham_VAT_8 <-  read.table("SRR8278869.txt", sep = "\t",header = T,quote = "")

sham_VAT_9 <-  read.table("SRR8278870.txt", sep = "\t",header = T,quote = "")

KPC_VAT_6 <- read.table("SRR8278871.txt", sep = "\t",header = T,quote = "")

KPC_VAT_7 <- read.table("SRR8278872.txt", sep = "\t",header = T,quote = "")

KPC_VAT_8 <-  read.table("SRR8278873.txt", sep = "\t",header = T,quote = "")

KPC_VAT_9 <-  read.table("SRR8278874.txt", sep = "\t",header = T,quote = "")




colnames(sham_VAT_6)[7]<-("sham_VAT_6")
colnames(sham_VAT_7)[7]<-("sham_VAT_7")
colnames(sham_VAT_8)[7]<-("sham_VAT_8")
colnames(sham_VAT_9)[7]<-("sham_VAT_9")
colnames(KPC_VAT_6)[7]<-("KPC_VAT_6")
colnames(KPC_VAT_7)[7]<-("KPC_VAT_7")
colnames(KPC_VAT_8)[7]<-("KPC_VAT_8")
colnames(KPC_VAT_9)[7]<-("KPC_VAT_9")

sham_VAT_6 <- sham_VAT_6[,c(1,7)]
sham_VAT_7 <- sham_VAT_7[,c(1,7)]
sham_VAT_8 <- sham_VAT_8[,c(1,7)]
sham_VAT_9 <- sham_VAT_9[,c(1,7)]

KPC_VAT_6 <- KPC_VAT_6[,c(1,7)]
KPC_VAT_7 <- KPC_VAT_7[,c(1,7)]
KPC_VAT_8 <- KPC_VAT_8[,c(1,7)]
KPC_VAT_9 <- KPC_VAT_9[,c(1,7)]


library(dplyr)

n_raw_count <- merge(sham_VAT_6,sham_VAT_7,by="Geneid") %>%
  merge(.,sham_VAT_8,by = "Geneid") %>%
  merge(.,sham_VAT_9,by = "Geneid") %>%

  merge(.,KPC_VAT_6,by = "Geneid") %>%
  merge(.,KPC_VAT_7,by = "Geneid") %>%
  merge(.,KPC_VAT_8,by = "Geneid") %>%
  merge(.,KPC_VAT_9,by = "Geneid")

library(tibble)

n_raw_count <- column_to_rownames(n_raw_count,"Geneid")


save("n_raw_count",file = "vat2_raw_count.RData")


raw_count_all <- merge(n_raw_count_1,n_raw_count,"Geneid")


library(tibble)

raw_count_all <- column_to_rownames(raw_count_all,"Geneid")

raw_count_all <- raw_count_all[,c(1:5,11:14, 6:10, 15:18)]

Condition <- factor(c(rep("sham_VAT",9),rep("KPC_VAT",9)), levels = c("sham_VAT","KPC_VAT"))


colData <- data.frame(row.names = colnames(raw_count_all),Condition)

colData$Batch <- factor(c(rep("Exp_1",5),rep("Exp_2",4),rep("Exp_1",5),rep("Exp_2",4)),levels = c("Exp_1","Exp_2"))

colData2 <- data.frame(row.names = colnames(raw_count_all),Condition)


expr_count_combat <- ComBat_seq(counts = as.matrix(raw_count_all), 
                                batch = colData$Batch  ) 
library(tinyarray)
library(tidyverse)

library(rtracklayer)


Gle <- read.table("SRR8278867.txt",skip = 1,sep="\t",header = T)

Gle  <- Gle[,c(1,6)]

le = Gle[match(rownames(expr_count_combat),Gle$Geneid),"Length"]

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

tpms <- apply(expr_count_combat,2,countToTpm,le)

write.table(exp_tpm,file = "exp_cac_VAT_tpm.txt",row.names = F,quote = F,sep = "\t")


exp_tpm = as.data.frame(tpms)

exp_tpm = rownames_to_column(exp_tpm)

##############################################
library(DESeq2)

dds <-DESeqDataSetFromMatrix(expr_count_combat,colData,design = ~ Condition)


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq2::DESeq(dds)

dat  <- counts(dds, normalized = TRUE)

idx  <- rowMeans(dat) > 1



resultsNames(dds)

res = DESeq2::results(dds ,c("Condition", "KPC_VAT","sham_VAT"))

res2 = res[order(res$padj),]


cut_pvalue <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < cut_pvalue),]

significant_results

significant_results_EXP <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_EXP = data.frame(significant_results_EXP)

significant_results_High <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange>cut_lfc)),]

significant_results_High <- rownames_to_column(as.data.frame(significant_results_High), "Gene")

significant_results_Low<- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) )),]

significant_results_Low <- rownames_to_column(as.data.frame(significant_results_Low), "Gene")



library(openxlsx)

write.xlsx(significant_results_High,"KPC_significant_High.xlsx")

write.xlsx(significant_results_Low,"KPC_significant_Low.xlsx")
###################################################################################################

########################################################
cut_lfc <- 1

cut_pvalue <- 0.05



topT <- as.data.frame(res)

topT <- rownames_to_column(topT,"Gene")

# Adjusted P values

rld <- DESeq2::rlog(dds, blind = FALSE)

rlogMat <- assay(rld)

png(filename="final_result/batch_box_plot.png",width=800,height=600,unit="px",bg="transparent",res = 100)

par(mar=c(7,5,5,2))

boxplot(rlogMat, las = 2, cex.axis = 1,main = "Normalized expression barplot", ylab= "Normalized expression value") +
  theme(axis.title.y = element_text(face="bold"),axis.text.y = element_text(face="bold"),legend.title = element_text(face="bold"),axis.text.x = element_text(face="bold"),axis.title.x =element_text(face="bold")) + 
  theme(text = element_text(face = "bold"))

dev.off()
####################################################################################################


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



canonical = gmtPathways("m2.cp.reactome.v2023.1.Mm.symbols.gmt")

canonical %>% 
  head() %>% 
  lapply(head)


canonical_Res <- fgsea(pathways=canonical, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)

GO = gmtPathways("m5.go.bp.v2023.1.Mm.symbols.gmt")

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

topGOPathwaysUp <- GO_Res[ES > 0][head(order(pval), n=20), pathway]

topGOPathwaysDown <- GO_Res[ES < 0][head(order(pval), n=20), pathway]

topPathways <- c(topGOPathwaysUp , rev(topGOPathwaysDown))

go_ResT_up_FOR = subset(go_ResT_up,go_ResT_up$pathway %in% topGOPathwaysUp)


go_ResT_dn_FOR = subset(go_ResT_dn,go_ResT_dn$pathway %in% topGOPathwaysDown)

go_ResT_all <- rbind(go_ResT_up_FOR,go_ResT_dn_FOR)

##go_ResT_all_df <- rbind(go_ResT_up,go_ResT_dn)


png(filename="final_result/VAT_ALLGO_BP_GSEA.png",width=1000,height=800,unit="px",bg="transparent",res = 100)


ggplot(go_ResT_all  , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES > 1.0)) +
  coord_flip() +
  labs(x="Pathway Terms", y="Normalized Enrichment Score",
       title="GO Biological Process from GSEA in KPC VAT") + 
  theme_bw()+ scale_fill_manual(name = "NES",labels = c("Negative","Positive")
                                ,values = c("TRUE" = "#F8766D","FALSE"= "#619CFF"))+
  theme(axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"), axis.text.x = element_text(face="bold"),axis.text.y = element_text(face="bold")) + 
  theme(text = element_text(face = "bold"))


dev.off()

#####################################################################################################
library(ggplot2)

library(envalysis)

library(fgsea)


GOBP_ERK <- plotEnrichment(GO[["GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE"]],stats = ranks)+ 
  labs(title="GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE") + theme_publish(base_size = 11, base_linewidth = 0.7) + theme(plot.title = element_text(size=6))

ggsave(filename = "final_result/ERK1_ERK2_GSEA.png", plot = GOBP_ERK, width = 12, height = 10, dpi = 300, units = "cm")


GOBP_PLC <-plotEnrichment(GO[["GOBP_PHOSPHOLIPASE_C_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"]],stats = ranks)+ 
  labs(title="GOBP_PHOSPHOLIPASE_C_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY") + theme_publish(base_size = 11, base_linewidth = 0.7) + theme(plot.title = element_text(size=6))


ggsave(filename = "final_result/PLC_GSEA.png", plot = GOBP_PLC, width = 12, height = 10, dpi = 300, units = "cm")



GOBP_ion <-plotEnrichment(GO[["GOBP_POSITIVE_REGULATION_OF_CYTOSOLIC_CALCIUM_ION_CONCENTRATION"]],stats = ranks)+ 
  labs(title="GOBP_POSITIVE_REGULATION_OF_CYTOSOLIC_CALCIUM_ION_CONCENTRATION") + theme_publish(base_size = 11, base_linewidth = 0.7) + theme(plot.title = element_text(size=6))


ggsave(filename = "final_result/ION_GSEA.png", plot = GOBP_ion, width = 12, height = 10, dpi = 300, units = "cm")

######################################################################################################

canonical_ResT <- canonical_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

canonical_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

canonical_ResT_up = subset(x = canonical_ResT, canonical_ResT$NES > 1)

canonical_ResT_up <- canonical_ResT_up[order(canonical_ResT_up$padj),]

canonical_ResT_up_sig <- subset(canonical_ResT_up,canonical_ResT_up$padj < 0.05)

library(openxlsx)

write.xlsx(canonical_ResT_up_sig,"canonical_up_significant_GSEA.xlsx")

canonical_ResT_dn = subset(x = canonical_ResT, canonical_ResT$NES < -1)

canonical_ResT_dn <- canonical_ResT_dn[order(canonical_ResT_dn$padj),]


canonical_ResT_dn <- canonical_ResT_dn[order(canonical_ResT_dn$padj),]

canonical_ResT_dn_sig <- subset(canonical_ResT_dn,canonical_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(canonical_ResT_dn_sig,"canonical_dn_significant_GSEA.xlsx")


topcanonicalPathwaysUp <- canonical_Res[ES > 0][head(order(pval), n=20), pathway]

topcanonicalPathwaysDown <- canonical_Res[ES < 0][head(order(pval), n=20), pathway]

topPathways <- c(topcanonicalPathwaysUp , rev(topcanonicalPathwaysDown))

canonical_ResT_up_FOR = subset(canonical_ResT_up,canonical_ResT_up$pathway %in% topcanonicalPathwaysUp)


canonical_ResT_dn_FOR = subset(canonical_ResT_dn,canonical_ResT_dn$pathway %in% topcanonicalPathwaysDown)

canonical_ResT_all <- rbind(canonical_ResT_up_FOR,canonical_ResT_dn_FOR)
########################################################################################################
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

library(tibble)



significant_results_up<- res2[which(res2$padj<cut_pvalue & res2$log2FoldChange>cut_lfc),]

significant_results_up

significant_results_sorted_up<- significant_results_up[order(significant_results_up$padj),]

significant_results_down<- res2[which(res2$padj<cut_pvalue & res2$log2FoldChange<(-cut_lfc)),]

significant_results_sorted_down<- significant_results_down[order(significant_results_down$padj),]


#Input genes; convert to ENTREZID 
eg_u= bitr(rownames(significant_results_sorted_up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

eg_d = bitr(rownames(significant_results_sorted_down), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


go_u_bp <- enrichGO(eg_u$SYMBOL, OrgDb = "org.Mm.eg.db", ont="BP",keyType = "SYMBOL")

go_d_bp <- enrichGO(eg_d$SYMBOL, OrgDb = "org.Mm.eg.db", ont="BP",keyType = "SYMBOL")


png(filename="final_result/GO_BP_KPC_UP_dotplot.png",width=1000,height=800,unit="px",bg="transparent")

dotplot(go_u_bp,showCategory = 20,font.size = 15,)

dev.off()

go_u_bp_df <- as.data.frame(go_u_bp@result)

go_u_bp_df <- subset(go_u_bp_df,go_u_bp_df$p.adjust < 0.05)

library(openxlsx)

write.xlsx(go_u_bp_df,"significant_KPC_VAT_UP_GO_BP.xlsx",rowNames = FALSE)



go_d_bp_df <- as.data.frame(go_d_bp@result)

go_d_bp_df <- subset(go_d_bp_df,go_d_bp_df$p.adjust < 0.05)


write.xlsx(go_d_bp_df,"significant_KPC_VAT_DN_GO_BP.xlsx",rowNames = FALSE)


png(filename="final_result/GO_BP_KPC_DN_dotplot.png",width=1000,height=800,unit="px",bg="transparent")

dotplot(go_d_bp,showCategory = 20,font.size = 15)

dev.off()

go_d_bp_df <- as.data.frame(go_d_bp@result)

go_d_bp_df <- subset(go_d_bp_df,go_d_bp_df$p.adjust < 0.05)

write.xlsx(go_d_bp_df,"significant_down_GO_BP.xlsx",rowNames = FALSE)


####CAA_VAT_UP <- read.xlsx("VAT_CAA_select_path_up_go.xlsx")


select_path <- CAA_VAT_UP$Pathway

go_u_select_path <- go_u_bp

go_u_select_path@result <-  go_u_bp@result[ go_u_bp@result$Description %in% select_path, ]
library(ggplot2)
library(enrichplot)

# dotplot 생성


png(filename="final_result/GO_BP_select_VAT.png",width=1000,height=900,unit="px",bg="transparent",res = 100)


dotplot(go_u_select_path, showCategory = length(select_path)) +
  ggtitle("Selected GO Biological Processes in KPC VAT") +
  theme_minimal(base_size = 14) +  # 전체 폰트 크기
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold")
  ) +
  scale_color_gradient(low = "#56B1F7", high = "#132B43") +  # 색 조절
  labs(
    x = "Gene Ratio", y = NULL, color = "Adjusted p-value", size = "Gene Count"
  )

dev.off()

#################################################################################################################

#remotes::install_github("omnideconv/immunedeconv")

library(immunedeconv)

h_tpms <- immunedeconv::convert_human_mouse_genes(tpms,convert_to = "human")

write.table(exp_tpm,file = "exp_cac_VAT_tpm.txt",row.names = F,quote = F,sep = "\t")


exp_tpm = as.data.frame(h_tpms)

exp_tpm = rownames_to_column(exp_tpm)

select_tpm <- h_tpms["P2RY6",c(10:18)]

source("CIBERSORT.R")

TME = CIBERSORT("LM22.txt", 
                "exp_cac_VAT_tpm.txt", 
                perm = 100, 
                QN = T)


cibersort_res <- TME[,-(23:25)]

library(pheatmap)

k <- apply(cibersort_res,2,function(x) {sum(x == 0) < nrow(TME)/2})

table(k)

cibersort_res2 <- as.data.frame(t(cibersort_res[,k]))

cibersort_res2_df <- rownames_to_column(cibersort_res2 ,"Celltypes")

nb.cols <- 18

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

mycolors <- c("#66C2A5", "#A3AC89","#E1917A" ,"#DCC199", "#929ECA", "#B795C7", "#C99592","#C7BAA6", "#EBC87C", "#B0D84F", "#D5D840", "#B3B3B3","#F5D152")
  

nb.cols <- 25

mycolors3 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

 
 library(tidyr)

png(filename="final_result/fractions_barplot_immunecells.png",width=1000,height=800,unit="px",bg="transparent",res = 150)


cibersort_res2_df%>%
  gather(sample, fraction, -Celltypes) %>%
  # plot as stacked bar chart
  ggplot(aes(x = sample, y = fraction, fill = Celltypes)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = mycolors) + 
  theme(axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"),
                                              axis.text.x = element_text(face="bold")) + 
  theme(text = element_text(face = "bold"))+ 
  scale_x_discrete(limits = rev(levels(cibersort_res2_df)))+ xlab("") + ylab("Proportion")

dev.off()



library(RColorBrewer)


cibersort_res_file <- t(cibersort_res2 )%>% as.data.frame()


cibersort_res_file <-   rownames_to_column(cibersort_res_file,"Sample") 

library(tibble)

cibersort_res2_df2 <- t(cibersort_res2 )%>% as.data.frame()

cibersort_res2_df2$Group <-  factor(c(rep("sham_VAT",9),rep("KPC_VAT",9)))

cibersort_res2_df2 <-   rownames_to_column(cibersort_res2_df2,"Sample") 

cell_types <- colnames(cibersort_res_file)[colnames(cibersort_res_file) != "Sample"]

t_test_results <- list()

for (cell in cell_types) {
  group1 <- cibersort_res2_df2[cibersort_res2_df2$Group == "sham_VAT", cell]
  group2 <- cibersort_res2_df2[cibersort_res2_df2$Group == "KPC_VAT", cell]
  
  test_result <- t.test(group1, group2)
  t_test_results[[cell]] <- test_result
}

for (cell in names(t_test_results)) {
  cat(sprintf("Cell type: %s\n", cell))
  print(t_test_results[[cell]])
  cat("\n")
}

summary_table <- data.frame(
  CellType = character(),
  p.value = numeric(),
  Mean_sham = numeric(),
  Mean_KPC = numeric(),
  stringsAsFactors = FALSE
)

for (cell in names(t_test_results)) {
  test <- t_test_results[[cell]]
  summary_table <- rbind(summary_table, data.frame(
    CellType = cell,
    p.value = test$p.value,
    Mean_sham = mean(cibersort_res2_df2[cibersort_res2_df2$Group == "sham_VAT", cell], na.rm = TRUE),
    Mean_KPC = mean(cibersort_res2_df2[cibersort_res2_df2$Group == "KPC_VAT", cell], na.rm = TRUE)
  ))
}

significant_summary <- subset(summary_table,summary_table$p.value < 0.05)


library(tibble)

library(dplyr)

library(reshape2)

dat_new2 = melt(cibersort_res2_df2)

colnames(dat_new2)=c("Sample","Group","Celltypes","Composition")

plot_order2 = dat_new2[dat_new2$Group=="KPC_VAT",] %>% 
  group_by(Celltypes) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltypes)

dat_new2$Celltypes = factor(dat_new2$Celltypes,levels = plot_order2)

dat_new2$Group<- relevel(dat_new2$Group, ref = "sham_VAT")

colnames(dat_new2)[2] <- "Condition"


png(filename="final_result/fractions_barplot2_immunecells.png",width=1000,height=800,unit="px",bg="transparent",res = 120)


ggplot(dat_new2, aes(x = Celltypes, y = Composition)) + 
  geom_boxplot(aes(fill = Condition)) + 
  theme_bw() + 
  labs(x = "Cell types", y = "Proportion of Immune infiltration") +
  theme(legend.position = "right") + 
  theme(axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"),
        axis.text.x = element_text(face="bold")) +
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  theme(text = element_text(face = "bold",size=10)) +
  scale_fill_manual(values = c("#619CFF", '#F8766D'))+ stat_compare_means(aes(group =  Condition),
                                                                          label = "p.signif",
                                                                          method = "t.test",
                                                                          size = 5,
                                                                          hide.ns = T)

dev.off()

significant_dat_new2 <- subset(dat_new2,dat_new2$Celltypes %in%significant_summary$CellType )

colnames(significant_dat_new2)[2] <- "Condition"

library(ggplot2)

library(ggstats)

library(ggpubr)

png(filename="final_result/fractions_barplot_significant_immunecells.png",width=1000,height=800,unit="px",bg="transparent",res = 130)

ggplot(significant_dat_new2, aes(x = Condition, y = Composition,fill = Condition))+ 
  labs(y="Cell Composition",x= NULL)+  
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  geom_boxplot(lwd = 1.5,position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)+ 
  theme_bw() + scale_fill_manual(values = c('sham_VAT'="#619CFF", 'KPC_VAT'='#F8766D')) +
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~ Celltypes,scales = "free",ncol = 4) + 
  theme(axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"),
        axis.text.x = element_text(face="bold")) +
  theme(text = element_text(face = "bold",size=10)) +
  stat_compare_means(aes(group =  Condition),
                     label = "p.format",
                     method = "t.test",
                     size = 3.5,
                     hide.ns = T,label.y.npc = 1.0,
                     label.x.npc = 0.4) 

dev.off()


select_tpm_log <- log2(select_tpm +1)

select_tpm_df <- data.frame(t(select_tpm_log))

cibersort_res2_df2_cac <- cibersort_res2_df2[10:18,]

select_tpm_df_cor <- cbind(select_tpm_df,cibersort_res2_df2_cac[,c(11,15)])

library(ggplot2)
library(ggpubr)
library(envalysis)


scatter_cor_M2 <- ggscatter( select_tpm_df_cor , x = "P2RY6", y = "Macrophages M2", xlab = "P2RY6 expression level in KPC VAT",
                              ylab = "Macrophages M2 in KPC VAT",
                              conf.int = TRUE,
                             add = "reg.line",  # Add regressin line
                             add.params = list(color = "red", fill = "lightgray"), # Customize reg. line
                              ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval




ggsave(filename = "final_result/P2RY6_M2_cor.png", plot = scatter_cor_M2, width = 12, height = 10, dpi = 300, units = "cm")
####################################################

library(openxlsx)

selected_GO_BP_UP <- read.xlsx("GO_selected_GSEA_UP_res.xlsx")

selected_GO_BP_DN <- read.xlsx("GO_selected_GSEA_DN_res.xlsx")

selected_GO_GSEA_all <- rbind(selected_GO_BP_UP,selected_GO_BP_DN)

###########################################################################

library(ggplot2)


png(filename="final_result/selected_GO_BP_GSEA.png",width=1200,height=1000,unit="px",bg="transparent",res = 100)


ggplot(selected_GO_GSEA_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES > 1.0)) +
  coord_flip() +
  labs(x="Pathway Terms", y="Normalized Enrichment Score",
       title="Selected GO Biological Process from GSEA in KPC VAT") + 
  theme_bw()+ scale_fill_manual(name = "NES",labels = c( "Negative","Positive")
                                     ,values = c("TRUE" = "#F8766D","FALSE"= "#619CFF"))+
  theme(axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"), axis.text.x = element_text(face="bold"),axis.text.y = element_text(face="bold")) + 
  theme(text = element_text(face = "bold"))


dev.off()


go_gsea_ggplot <-  ggplot(selected_GO_GSEA_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES > 1.0)) +
  coord_flip() +
  labs(x="Pathway Terms", y="Normalized Enrichment Score",
       title="GO Biological Process from GSEA") + 
  theme_bw()+ scale_fill_manual(name = "NES",labels = c("Negative", "Positive")
                                ,values = c("TRUE" = "#F8766D","FALSE"= "#619CFF"))+
  theme(axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"),legend.title = element_text(face="bold"),axis.text.x =element_text(face="bold") )



ggsave(filename = "final_result/GO_GSEA2.png", plot = go_gsea_ggplot, width = 12, height = 10, dpi = 300, units = "cm")
####################################################

overlap_high_sig_exp_genes <- read.table("overlap_high_sigExp.txt", sep = "\t",header = TRUE,stringsAsFactors = FALSE)



library(tibble)

toptt <- column_to_rownames(topT,"Gene")


VAT_res <- toptt %>%
  mutate(Condition= case_when(log2FoldChange >  1 & padj < 0.05 ~ "KPC_VAT",
                               log2FoldChange< -1 & padj < 0.05 ~ "sham_VAT",
                               TRUE ~ "Stable"))   



VAT_res$Condition <- factor(VAT_res$Condition,levels = c("KPC_VAT","Stable","sham_VAT"))

library(ggplot2)

library(ggpubr)

library(ggrepel)

options(ggrepel.max.overlaps = Inf)


png(filename="final_result/KPC_volcano_plot2.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

ggplot(data =VAT_res,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  theme_classic() +
  geom_point(aes(colour = Condition), 
             alpha = 0.2, 
             shape = 16,
             size = 3) +
  scale_color_manual(values = c("sham_VAT" = "#619CFF","Stable" = "grey","KPC_VAT" = "#F8766D")) +
  geom_point(data = VAT_res[overlap_high_sig_exp_genes$Gene, ],
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
  geom_label_repel(data = VAT_res[overlap_high_sig_exp_genes$Gene, ],     
                   aes(label = overlap_high_sig_exp_genes$Gene),
                   force = 2,
                   nudge_y = 1,
                   size = 3)  + ylab("-log10(adj p value)") + xlab("Log2FoldChange")+ guides(colour = guide_legend(override.aes = list(size=3)))+ ggtitle('KPC VAT vs sham VAT volcano plot')+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))  +
  theme(axis.title.y = element_text(face="bold",size = 15),axis.text.y = element_text(face="bold",size = 15),legend.title = element_text(face="bold"),axis.text.x = element_text(face="bold",size = 15),axis.title.x =element_text(face="bold",size = 15))



dev.off()
##################################################2507IMMUNE_RESPONSEcor


ACTIVATION_OF_IMMUNE_RESPONSE <- GO$GOBP_ACTIVATION_OF_IMMUNE_RESPONSE
ACUTE_INFLAMMATORY_RESPONSE<- GO$GOBP_ACUTE_INFLAMMATORY_RESPONSE

Gene_set_ssgsea = list(ACTIVATION_OF_IMMUNE_RESPONSE,ACUTE_INFLAMMATORY_RESPONSE)


names(Gene_set_ssgsea) = c("ACTIVATION_OF_IMMUNE_RESPONSE","ACUTE_INFLAMMATORY_RESPONSE")
exp_tpm 


exp_tpm_log <- log2(tpms +1)


select_gene_cac_tpm_df <- data.frame( P2ry6 = (exp_tpm_log["P2ry6",10:18]))

select_gene_cac_tpm_df <- rownames_to_column(select_gene_cac_tpm_df,"Sample")

library(GSVA)

ssgsea_immune = gsva(exp_tpm_log, Gene_set_ssgsea,method = "ssgsea")



ssgsea_immune = as.data.frame(t(ssgsea_immune))

ssgsea_immune = rownames_to_column(ssgsea_immune,"Sample")


ssgsea_immune = merge(ssgsea_immune,select_gene_cac_tpm_df  ,by = "Sample")

ssgsea_immune  = column_to_rownames(ssgsea_immune ,"Sample")



library(tibble)


library(ggpubr)
library(envalysis)

scatter_cor_inflamm <- ggscatter( ssgsea_immune , x = "P2ry6", y = "ACUTE_INFLAMMATORY_RESPONSE", xlab = "P2RY6 expression level in KPC VAT",
                                 ylab = "INFLAMMATORY_RESPONSE in KPC VAT",
                                 conf.int = TRUE,
                                 add = "reg.line",  # Add regressin line
                                 add.params = list(color = "red", fill = "lightgray"), # Customize reg. line
                                 ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval


ggsave(filename = "final_result/P2RY6_inflamm_cor.png", plot = scatter_cor_inflamm, width = 12, height = 10, dpi = 300, units = "cm")


scatter_cor_immune <- ggscatter( ssgsea_immune , x = "P2ry6", y = "ACTIVATION_OF_IMMUNE_RESPONSE", xlab = "P2RY6 expression level in KPC VAT",
                             ylab = "Immune response in KPC VAT",
                             conf.int = TRUE,
                             add = "reg.line",  # Add regressin line
                             add.params = list(color = "red", fill = "lightgray"), # Customize reg. line
                             ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval




ggsave(filename = "final_result/P2RY6_immune_cor.png", plot = scatter_cor_immune, width = 12, height = 10, dpi = 300, units = "cm")

####################################################
select_tpm <- h_tpms["P2RY6",c(10:18)]




inflammatory = c("Tnf","Il6","Il1b")

inflammatory_exp  = rlogMat[inflammatory,]

inflammatory_exp  <- as.data.frame(t(inflammatory_exp))

inflammatory_exp$Condition <- c(rep("sham_VAT",9),rep("KPC_VAT",9))


library(tidyr)

library(tibble)


inflammatory_exp<- inflammatory_exp%>%
  rownames_to_column("Sample") %>% 
  gather(key = Gene,value = Proportion,-Sample,-Condition)


inflammatory_exp$Condition <- factor(inflammatory_exp$Condition,levels = c("sham_VAT","KPC_VAT"))



library(ggplot2)

library(ggpubr)

library(rstatix)


png(filename="final_result/inflammatory_violin.png",width=800,height=600,unit="px",bg="transparent",res = 100)

ggplot(inflammatory_exp, aes(x=Gene, y=Proportion,fill = Condition)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(
    aes(group = Condition),
    method = "t_test", label = "p.signif",size = 1,label.size = 10) +
  geom_boxplot(lwd = 1.5,position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE) +  theme(text = element_text(face = "bold",size=20)) + ylab("Normalized expression level") + scale_fill_manual(values = c('sham_VAT'="#619CFF", 'KPC_VAT'='#F8766D'))

dev.off()
