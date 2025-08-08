sham_muscle_1 <- read.table("SRR8278856.txt", sep = "\t",header = T,quote = "")

sham_muscle_2 <- read.table("SRR8278857.txt", sep = "\t",header = T,quote = "")

sham_muscle_3 <-  read.table("SRR8278858.txt", sep = "\t",header = T,quote = "")

KPC_muscle_1 <-  read.table("SRR8278859.txt", sep = "\t",header = T,quote = "")

KPC_muscle_2<- read.table("SRR8278860.txt", sep = "\t",header = T,quote = "")

KPC_muscle_3<- read.table("SRR8278861.txt", sep = "\t",header = T,quote = "")

KPC_muscle_4 <-  read.table("SRR8278862.txt", sep = "\t",header = T,quote = "")



colnames(sham_muscle_1)[7]<-("sham_muscle_1")
colnames(sham_muscle_2)[7]<-("sham_muscle_2")
colnames(sham_muscle_3)[7]<-("sham_muscle_3")
colnames(KPC_muscle_1)[7]<-("KPC_muscle_1")
colnames(KPC_muscle_2)[7]<-("KPC_muscle_2")
colnames(KPC_muscle_3)[7]<-("KPC_muscle_3")
colnames(KPC_muscle_4)[7]<-("KPC_muscle_4")


sham_muscle_1 <- sham_muscle_1[,c(1,7)]
sham_muscle_2 <- sham_muscle_2[,c(1,7)]
sham_muscle_3 <- sham_muscle_3[,c(1,7)]

KPC_muscle_1 <- KPC_muscle_1[,c(1,7)]
KPC_muscle_2 <- KPC_muscle_2[,c(1,7)]
KPC_muscle_3 <- KPC_muscle_3[,c(1,7)]
KPC_muscle_4 <- KPC_muscle_4[,c(1,7)]


library(dplyr)

n_raw_count <- merge(sham_muscle_1,sham_muscle_2,by="Geneid") %>%
  merge(.,sham_muscle_3,by = "Geneid") %>%
  
  merge(.,KPC_muscle_1,by = "Geneid") %>%
  merge(.,KPC_muscle_2,by = "Geneid") %>%
  merge(.,KPC_muscle_3,by = "Geneid") %>%
  merge(.,KPC_muscle_4,by = "Geneid")


save("n_raw_count",file = "muscle_raw_count.RData")

library(tibble)

n_raw_count <- column_to_rownames(n_raw_count,"Geneid")

##############################################################


Groups <- factor(c(rep("sham_muscle",3),rep("KPC_muscle",4)), levels = c("sham_muscle","KPC_muscle"))


colData <- data.frame(row.names = colnames(n_raw_count),Groups)


colData$Groups = relevel(colData$Groups, ref = "sham_muscle")


library(DESeq2)

dds <-DESeqDataSetFromMatrix(n_raw_count,colData,design = ~ Groups )


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq2::DESeq(dds)

dat  <- counts(dds, normalized = TRUE)


res = DESeq2::results(dds ,c("Groups", "KPC_muscle","sham_muscle"))

res2 = res[order(res$padj),]

res2_df <- data.frame(res2)

res2_df <- rownames_to_column(res2_df,"Gene")

cut_pvalue <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < cut_pvalue),]

significant_results

significant_results_EXP <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_EXP = data.frame(significant_results_EXP)

significant_results_EXP = rownames_to_column(significant_results_EXP,"Gene")

significant_results_High <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange>cut_lfc)),]

significant_results_High <- rownames_to_column(as.data.frame(significant_results_High), "Gene")


significant_results_Low<- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) )),]

significant_results_Low <- rownames_to_column(as.data.frame(significant_results_Low), "Gene")
############################################################################################
library(openxlsx)

cocul_KPC_high <- read.xlsx("Adipocyte_cocultured_significant_kpc_High2.xlsx") 


gene <- cocul_KPC_high$Gene


toptt <- column_to_rownames(topT,"Gene")


gene_df <- data.frame(Gene = gene)

cocul_KPC_high_muscle_high <- merge(significant_results_High,gene_df)

write.xlsx(cocul_KPC_high_muscle_high,"cocul_KPC_high_muscle_high.xlsx")

cocul_KPC_high_muscle_low <- merge(significant_results_Low,gene_df)


write.xlsx(cocul_KPC_high_muscle_low,"cocul_KPC_high_muscle_low.xlsx")
################################################################################################

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

#################################################################################

library(DESeq2)

rld <- DESeq2::rlog(dds, blind = FALSE)

rlogMat <- assay(rld)

muscle_atrophy = c("Fbxo32","Trim63")


atrophy_exp  = rlogMat[muscle_atrophy,]

atrophy_exp <- as.data.frame(t(atrophy_exp))

atrophy_exp$Condition <- c(rep("sham_muscle",3),rep("KPC_muscle",4))


library(tidyr)

library(tibble)


atrophy_exp<- atrophy_exp%>%
  rownames_to_column("Sample") %>% 
  gather(key = Gene,value = Proportion,-Sample,-Condition)

atrophy_exp$Condition <- factor(atrophy_exp$Condition,levels = c("sham_muscle","KPC_muscle"))

library(ggplot2)

library(ggpubr)

library(rstatix)


png(filename="final_result/muscle_related_violin.png",width=800,height=600,unit="px",bg="transparent")

atrophys <- ggplot(atrophy_exp, aes(x=Gene, y=Proportion,fill = Condition)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(
    aes(group = Condition),
    method = "t_test", label = "p.signif",size = 1,label.size = 10) +
  geom_boxplot(lwd = 1.5,position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE) +  theme(text = element_text(face = "bold",size=20)) + ylab("Normalized expression level") + scale_fill_manual(values = c('sham_muscle'="#619CFF", 'KPC_muscle'='#F8766D'))


muscle_receptor = c("Acvr1b","Il6ra","Tnfrsf12a")


receptor_exp  = rlogMat[muscle_receptor,]

receptor_exp<- as.data.frame(t(receptor_exp))

receptor_exp$Condition <- c(rep("sham_muscle",3),rep("KPC_muscle",4))


library(tidyr)

library(tibble)

library(dplyr)

receptor_exp<- receptor_exp%>%
  rownames_to_column("Sample") %>% 
  gather(key = Gene,value = Proportion,-Sample,-Condition)


receptor_exp$Condition <- factor(receptor_exp$Condition,levels = c("sham_muscle","KPC_muscle"))

library(ggplot2)

library(ggpubr)

library(rstatix)
receptors <- ggplot(receptor_exp, aes(x=Gene, y=Proportion,fill = Condition)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(
    aes(group = Condition),
    method = "t_test", label = "p.signif",size = 1,label.size = 10) +
  geom_boxplot(lwd = 1.5,position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE) +  theme(text = element_text(face = "bold",size=20)) + ylab("Normalized expression level") + scale_fill_manual(values = c('sham_muscle'="#619CFF", 'KPC_muscle'='#F8766D'))


ggarrange(atrophys , receptors,
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 20))

dev.off()

Fibrosis <- c("Col6a1","Col3a1","Col1a1","Fn1","Cd63","Timp1","P2ry6")

pheatmap(rlogMat[Fibrosis,]  ,
         annotation_col = annot_colData,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         fontsize = 8,
         scale = "row",
         fontsize_row = 8,
         annotation_colors = list(Groups=c(sham_muscle = "#619CFF",KPC_muscle='#F8766D')),
         color=colorRampPalette(c("navy", "white", "red"))(50))

         inflammatory = c("Tnf","Il6","Il1b")

inflammatory_exp  = rlogMat[inflammatory,]

inflammatory_exp  <- as.data.frame(t(inflammatory_exp))

inflammatory_exp$Condition <- c(rep("sham_muscle",3),rep("KPC_muscle",4))


library(tidyr)

library(tibble)


inflammatory_exp<- inflammatory_exp%>%
  rownames_to_column("Sample") %>% 
  gather(key = Gene,value = Proportion,-Sample,-Condition)


inflammatory_exp$Condition <- factor(inflammatory_exp$Condition,levels = c("sham_muscle","KPC_muscle"))



library(ggplot2)

library(ggpubr)

library(rstatix)


png(filename="final_result/inflammatory_violin.png",width=800,height=600,unit="px",bg="transparent",res = 100)

ggplot(inflammatory_exp, aes(x=Gene, y=Proportion,fill = Condition)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(
    aes(group =Condition),
    method = "t_test", label = "p.signif",size = 1,label.size = 10) +
  geom_boxplot(lwd = 1.5,position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE) +  theme(text = element_text(face = "bold",size=20)) + ylab("Normalized expression level") + scale_fill_manual(values = c('sham_muscle'="#619CFF", 'KPC_muscle'='#F8766D'))

dev.off()

library(openxlsx)

high_overlapped <- read.xlsx("overlap_high_sigExp.xlsx")

library(dplyr)
library(tibble)

toptt <- column_to_rownames(topT,"Gene")


muscle_res <- toptt %>%
  mutate(Condition= case_when(log2FoldChange >  1 & padj < 0.05 ~ "KPC_muscle",
                              log2FoldChange< -1 & padj < 0.05 ~ "sham_muscle",
                              TRUE ~ "Stable"))   




muscle_res$Condition <- factor(muscle_res$Condition,levels = c("KPC_muscle","Stable","sham_muscle"))

library(ggplot2)

library(ggpubr)

library(ggrepel)


options(ggrepel.max.overlaps = Inf)


png(filename="final_result/KPC_muslce_plot2.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

ggplot(data =muscle_res,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
  theme_classic() +
  geom_point(aes(colour = Condition), 
             alpha = 0.2, 
             shape = 16,
             size = 3) +
  scale_color_manual(values = c("sham_muscle" = "#619CFF","Stable" = "grey","KPC_muscle" = "#F8766D")) +
  geom_point(data = muscle_res[high_overlapped$Gene, ],
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
  geom_label_repel(data = muscle_res[high_overlapped$Gene, ],     
                   aes(label = high_overlapped$Gene),
                   force = 2,
                   nudge_y = 1,
                   size = 3) +  theme(text = element_text(size=3)) +  theme(text = element_text(size=10)) + ylab("-log10(adj p value)") + xlab("Log2FoldChange")+ guides(colour = guide_legend(override.aes = list(size=3)))+ ggtitle('KPC muscle vs sham muscle volcano plot')+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))  +
  theme(axis.title.y = element_text(face="bold",size = 15),axis.text.y = element_text(face="bold",size = 15),legend.title = element_text(face="bold"),axis.text.x = element_text(face="bold",size = 15),axis.title.x =element_text(face="bold",size = 15))

dev.off()



