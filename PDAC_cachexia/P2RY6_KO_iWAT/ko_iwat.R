P2RY6_KO_1 <- read.table("SRR11103992.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_2 <- read.table("SRR11103991.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_3 <-  read.table("SRR11104003.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_4  <- read.table("SRR11104002.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_5  <- read.table("SRR11104001.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_6 <-  read.table("SRR11104000.txt", sep = "\t",header = T,quote = "")

P2RY6_KO_7 <-  read.table("SRR11103999.txt", sep = "\t",header = T,quote = "")

WAT_1 <- read.table("SRR11104005.txt", sep = "\t",header = T,quote = "")

WAT_2 <- read.table("SRR11104004.txt", sep = "\t",header = T,quote = "")

WAT_3 <-  read.table("SRR11103998.txt", sep = "\t",header = T,quote = "")

WAT_4  <- read.table("SRR11103997.txt", sep = "\t",header = T,quote = "")

WAT_5  <- read.table("SRR11103996.txt", sep = "\t",header = T,quote = "")

WAT_6 <-  read.table("SRR11103995.txt", sep = "\t",header = T,quote = "")

WAT_7 <-  read.table("SRR11103994.txt", sep = "\t",header = T,quote = "")

WAT_8 <-  read.table("SRR11103993.txt", sep = "\t",header = T,quote = "")




colnames(P2RY6_KO_1)[7]<-("P2RY6_KO_1")
colnames(P2RY6_KO_2)[7]<-("P2RY6_KO_2")
colnames(P2RY6_KO_3)[7]<-("P2RY6_KO_3")
colnames(P2RY6_KO_4)[7]<-("P2RY6_KO_4")
colnames(P2RY6_KO_5)[7]<-("P2RY6_KO_5")
colnames(P2RY6_KO_6)[7]<-("P2RY6_KO_6")
colnames(P2RY6_KO_7)[7]<-("P2RY6_KO_7")



colnames(WAT_1)[7]<-("WAT_1")
colnames(WAT_2)[7]<-("WAT_2")
colnames(WAT_3)[7]<-("WAT_3")
colnames(WAT_4)[7]<-("WAT_4")
colnames(WAT_5)[7]<-("WAT_5")
colnames(WAT_6)[7]<-("WAT_6")
colnames(WAT_7)[7]<-("WAT_7")
colnames(WAT_8)[7]<-("WAT_8")


P2RY6_KO_1 <- P2RY6_KO_1[,c(1,7)]
P2RY6_KO_2 <- P2RY6_KO_2[,c(1,7)]
P2RY6_KO_3<- P2RY6_KO_3[,c(1,7)]
P2RY6_KO_4<- P2RY6_KO_4[,c(1,7)]
P2RY6_KO_5<- P2RY6_KO_5[,c(1,7)]
P2RY6_KO_6<- P2RY6_KO_6[,c(1,7)]
P2RY6_KO_7 <- P2RY6_KO_7[,c(1,7)]


WAT_1 <- WAT_1[,c(1,7)]
WAT_2 <- WAT_2[,c(1,7)]
WAT_3<- WAT_3[,c(1,7)]
WAT_4<- WAT_4[,c(1,7)]
WAT_5<- WAT_5[,c(1,7)]
WAT_6<- WAT_6[,c(1,7)]
WAT_7 <- WAT_7[,c(1,7)]
WAT_8 <- WAT_8[,c(1,7)]


library(dplyr)


n_raw_count <- merge(WAT_1,WAT_2,by="Geneid") %>%
  merge(.,WAT_3,by = "Geneid") %>%
  merge(.,WAT_4,by = "Geneid") %>%
  merge(.,WAT_5,by = "Geneid") %>%
  merge(.,WAT_6,by = "Geneid") %>%
  merge(.,WAT_7,by = "Geneid") %>%
  merge(.,WAT_8,by = "Geneid") %>%
  merge(.,P2RY6_KO_1,by = "Geneid") %>%
  merge(.,P2RY6_KO_2,by = "Geneid") %>%
  merge(.,P2RY6_KO_3,by = "Geneid") %>%
  merge(.,P2RY6_KO_4,by = "Geneid") %>%
  merge(.,P2RY6_KO_5,by = "Geneid") %>%
  merge(.,P2RY6_KO_6,by = "Geneid") %>%
  merge(.,P2RY6_KO_7,by = "Geneid")


library(tibble)

n_raw_count <- column_to_rownames(n_raw_count,"Geneid")


condition <- factor(c(rep("WAT",8),rep("P2RY6_KO",7)), levels = c("WAT","P2RY6_KO"))

colData <- data.frame(row.names = colnames(n_raw_count),condition)


colData$condition = relevel(colData$condition, ref = "WAT")


library(DESeq2)


dds <- DESeq2::DESeqDataSetFromMatrix(n_raw_count,colData,design = ~ condition )
  

dds <- DESeq2::DESeq(dds)



res = DESeq2::results(dds ,c("condition", "P2RY6_KO","WAT"))

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
eg_u= bitr(rownames(significant_results_sorted_up), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

eg_d = bitr(rownames(significant_results_sorted_down), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")


go_u_bp <- enrichGO(eg_u$SYMBOL, OrgDb = "org.Mm.eg.db", ont="BP",keyType = "SYMBOL")


#go_u_all <- enrichGO(eg_u$SYMBOL, OrgDb = "org.Mm.eg.db", ont="all",keyType = "SYMBOL")


go_d_bp <- enrichGO(eg_d$SYMBOL, OrgDb = "org.Mm.eg.db", ont="BP",keyType = "SYMBOL")



dotplot(go_u_all , split="ONTOLOGY", showCategory=10) + facet_grid(ONTOLOGY~., scale="free")

dd <- go_u_all@result

dotplot(go_u_bp,showCategory = 30,font.size = 15)

go_u_bp_df <- as.data.frame(go_u_bp@result)

go_u_bp_df <- subset(go_u_bp_df,go_u_bp_df$p.adjust < 0.05)


library(openxlsx)

write.xlsx(go_u_bp_df,"significant_UP_GO_BP.xlsx",rowNames = FALSE)

dotplot(go_u_bp,showCategory = 30,font.size = 10)


dotplot(go_d_bp,showCategory = 30,font.size = 10)


go_d_bp_df <- as.data.frame(go_d_bp@result)

go_d_bp_df <- subset(go_d_bp_df,go_d_bp_df$p.adjust < 0.05)


write.xlsx(go_d_bp_df,"significant_down_GO_BP.xlsx",rowNames = FALSE)

significant_results_EXP_df = rownames_to_column(significant_results_EXP,"Gene")


significant_results_EXP_High <- subset(significant_results_EXP_df,significant_results_EXP_df$log2FoldChange > 1)

library(openxlsx)

write.xlsx(significant_results_EXP_High,"significant_High_P2RY6_KO.xlsx")


significant_results_EXP_Low <- subset(significant_results_EXP_df,significant_results_EXP_df$log2FoldChange < -1)

write.xlsx(significant_results_EXP_Low,"significant_Low_P2RY6_KO.xlsx")

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


go_ResT_dn = subset(x = go_ResT, go_ResT$NES < -1)

go_ResT_dn <- go_ResT_dn[order(go_ResT_dn$padj),]

go_ResT_dn_sig <- subset(go_ResT_dn,go_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(go_ResT_up_sig,"GO_BP_up_significant_GSEA.xlsx")

write.xlsx(go_ResT_dn_sig,"GO_BP_dn_significant_GSEA.xlsx")

