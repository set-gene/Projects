adipo_alone_1 <- read.table("SRR8329153.txt", sep = "\t",header = T,quote = "")

adipo_alone_2 <- read.table("SRR8329154.txt", sep = "\t",header = T,quote = "")

adipo_alone_3 <-  read.table("SRR8329155.txt", sep = "\t",header = T,quote = "")

adipo_co_1 <- read.table("SRR8329156.txt", sep = "\t",header = T,quote = "")

adipo_co_2 <- read.table("SRR8329157.txt", sep = "\t",header = T,quote = "")

adipo_co_3 <-  read.table("SRR8329158.txt", sep = "\t",header = T,quote = "")

colnames(adipo_alone_1)[7]<-("adipo_alone_1")
colnames(adipo_alone_2)[7]<-("adipo_alone_2")
colnames(adipo_alone_3)[7]<-("adipo_alone_3")
colnames(adipo_co_1)[7]<-("adipo_co_1")
colnames(adipo_co_2)[7]<-("adipo_co_2")
colnames(adipo_co_3)[7]<-("adipo_co_3")

adipo_alone_1 <- adipo_alone_1[,c(1,7)]
adipo_alone_2 <- adipo_alone_2[,c(1,7)]
adipo_alone_3 <- adipo_alone_3[,c(1,7)]
adipo_co_1 <- adipo_co_1[,c(1,7)]
adipo_co_2 <- adipo_co_2[,c(1,7)]
adipo_co_3 <- adipo_co_3[,c(1,7)]

library(dplyr)

n_raw_count <- merge(adipo_alone_1,adipo_alone_2,by="Geneid") %>%
  merge(.,adipo_alone_3,by = "Geneid") %>%
  merge(.,adipo_co_1,by = "Geneid") %>%
  merge(.,adipo_co_2,by = "Geneid") %>%
  merge(.,adipo_co_3,by = "Geneid")

library(tibble)

n_raw_count <- column_to_rownames(n_raw_count,"Geneid")


condition <- factor(c(rep("Alone",3),rep("Cocultured",3)), levels = c("Alone","Cocultured"))

colData <- data.frame(row.names = colnames(n_raw_count),condition)


colData$condition = relevel(colData$condition, ref = "Alone")


library(DESeq2)


dds <- DESeq2::DESeqDataSetFromMatrix(n_raw_count,colData,design = ~ condition )


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq2::DESeq(dds)



res = DESeq2::results(dds ,c("condition", "Cocultured","Alone"))

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


go_u_bp_df <- as.data.frame(go_u_bp@result)

go_u_bp_df <- subset(go_u_bp_df,go_u_bp_df$p.adjust < 0.05)


library(openxlsx)

write.xlsx(go_u_bp_df,"significant_UP_GO_BP.xlsx",rowNames = FALSE)

go_d_bp_df <- as.data.frame(go_d_bp@result)

go_d_bp_df <- subset(go_d_bp_df,go_d_bp_df$p.adjust < 0.05)


write.xlsx(go_d_bp_df,"significant_down_GO_BP.xlsx",rowNames = FALSE)

library(openxlsx)

library(tibble)
##########################################################################
KPC_high <- read.xlsx("KPC_significant_High.xlsx")

KPC_low <- read.xlsx("KPC_significant_Low.xlsx")


KPC_high_genes = data.frame(Gene = KPC_high$Gene)

KPC_low_genes = data.frame(Gene = KPC_low$Gene)

significant_results_KPC_High_select <- merge(significant_results_EXP_High,KPC_high_genes,"Gene")

significant_results_KPC_low_select <- merge(significant_results_EXP_Low,KPC_low_genes,"Gene")

significant_results_KPC <- rbind(significant_results_KPC_High_select, significant_results_KPC_low_select)

significant_results_KPC_genes <- data.frame(Gene = significant_results_KPC$Gene)

library(openxlsx)

write.xlsx(significant_results_KPC_High_select,"final_result/Adipocyte_cocultured_significant_kpc_High.xlsx")

write.xlsx(significant_results_KPC_low_select,"Adipocyte_cocultured_significant_kpc_Low.xlsx")

############################################################################

topT <- as.data.frame(res)

topT <- rownames_to_column(topT,"Gene")

topT_kpc <- merge(topT , significant_results_KPC_genes,by = "Gene")

topT_kpc

$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

library(openxlsx)

library(dplyr)



CAA_VAT_UP <- read.xlsx("VAT_CAA_select_path_up_go.xlsx")


select_path <- CAA_VAT_UP$Pathway

go_u_select_path <- go_u_bp

go_u_select_path@result <-  go_u_bp@result[ go_u_bp@result$Description %in% select_path, ]
library(ggplot2)
library(enrichplot)

png(filename="final_result/GO_BP_select_CAA.png",width=1000,height=900,unit="px",bg="transparent",res = 100)


dotplot(go_u_select_path, showCategory = length(select_path)) +
  ggtitle("Selected GO Biological Processes in CAA") +
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

#################################################################################################

CAA_high_genes <- data.frame(Gene = significant_results_EXP_High$Gene)

KPC_high_genes


dat_high = list(
  A = KPC_high_genes$Gene,
  B = CAA_high_genes$Gene)

names(dat_high) <- c("KPC_VAT vs Sham_VAT","CAA vs Alone ")

library(ggVennDiagram)

vat_caa_venn <- ggVennDiagram(dat_high, label_alpha = 0, set_size = 3,
              category.names = c("KPC_VAT vs VAT","CAA vs Alone"))+
  scale_fill_gradient(low = "skyblue", high = "red")+
  scale_x_continuous(expand = expansion(mult = .6)) 

library(officer)
library(rvg)
#####################################
read_pptx() %>% 
  print(target = "./cacVAT_CAA_venn.pptx")

read_pptx( "./cacVAT_CAA_venn.pptx")


editable_graph <- rvg::dml(ggobj = vat_caa_venn )

read_pptx( "./cacVAT_CAA_venn.pptx") %>%
  add_slide("Title Slide","Office Theme") %>%
  ph_with(editable_graph,
          location = ph_location(left=0, top=0,width = 5, height = 5, bg="transparent")) %>%
  print(target = "./cacVAT_CAA_venn.pptx")

###################################
library(ggvenn)

library(openxlsx)

significant_results_EXP_KO <- read.xlsx("significant_Low_P2RY6_KO.xlsx")


names(dat_overlap) <- c("KPC_VAT vs Sham_VAT","CAA vs Alone ","P2RY6_KO_iWAT vs iWAT")


dat_overlap = list(
  A = KPC_high_genes$Gene,
  B = CAA_high_genes$Gene,
  C = significant_results_EXP_KO$Gene)

read_pptx() %>% 
  print(target = "./cacVAT_CAA_venn_KOiwat.pptx")

read_pptx( "./cacVAT_CAA_venn_KOiwat.pptx")


editable_graph2 <- rvg::dml(ggobj = all_ovlap_venn )

read_pptx( "./cacVAT_CAA_venn_KOiwat.pptx") %>%
  add_slide("Title Slide","Office Theme") %>%
  ph_with(editable_graph2,
          location = ph_location(left=0, top=0,width = 5, height = 5, bg="transparent")) %>%
  print(target = "./cacVAT_CAA_venn_KOiwat.pptx")
###############################################
significant_results_KPC_High_select
significant_results_EXP_KO
overlap_genes <- merge(significant_results_KPC_High_select,significant_results_EXP_KO,"Gene")
###############################################


#convertMouseGeneList(Gene)

library(biomaRt)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = Gene , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

genesV2

select_gene <- unique(genesV2$MGI.symbol)

select_gene

overlap_high_sig_exp_df_only <- subset(overlap_high_sig_exp_df,overlap_high_sig_exp_df$Gene %in% select_gene) 

overlap_high_sig_exp_df_only <- overlap_high_sig_exp_df_only [c(order(overlap_high_sig_exp_df_only$CAA_L2FC,decreasing = TRUE)),]


overlap_high_sig_exp_df_only$Gene[1:10]


write.table(overlap_high_sig_exp_df,"overlap_high_sigExp.txt",row.names = FALSE,quote = FALSE,sep = "\t")

write.xlsx(overlap_high_sig_exp_df,"overlap_high_sigExp.xlsx")


library(ggrepel)


library(ggplot2)

library(dplyr)

CAA_res <- toptt %>%
  mutate(Condition = case_when(log2FoldChange >  1 & padj < 0.05 ~ "Cocultured",
                               log2FoldChange< -1 & padj < 0.05 ~ "Alone",
                               TRUE ~ "Stable"))   



CAA_res$Condition <- factor(CAA_res$Condition,levels = c("Cocultured","Stable","Alone"))


options(ggrepel.max.overlaps = Inf)


png(filename="final_result/caa2_volcano_plot2.png",width=1000,height=800,unit="px",bg="transparent",res = 100)

ggplot(data =CAA_res,
       aes(x = log2FoldChange,
           y = -log10(padj))) +
   theme_classic() +
  geom_point(aes(colour = Condition), 
             alpha = 0.2, 
             shape = 16,
             size = 3) +
  scale_color_manual(values = c("Alone" = "#619CFF","Stable" = "grey","Cocultured" = "#F8766D")) +
  geom_point(data = CAA_res[overlap_high_sig_exp_df$Gene, ],
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
  geom_label_repel(data = CAA_res[overlap_high_sig_exp_df$Gene, ],     
                   aes(label = overlap_high_sig_exp_df$Gene),
                   force = 2,
                   nudge_y = 1,
                   size = 3)  + ylab("-log10(adj p value)") + xlab("Log2FoldChange")+ guides(colour = guide_legend(override.aes = list(size=3)))+ ggtitle('Cocultured vs Alone volcano plot')+
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold'))  +
  theme(axis.title.y = element_text(face="bold",size = 15),axis.text.y = element_text(face="bold",size = 15),legend.title = element_text(face="bold"),axis.text.x = element_text(face="bold",size = 15),axis.title.x =element_text(face="bold",size = 15))



dev.off()

library(openxlsx)

write.xlsx(genesV2,"overlap_VAT_CAA_high_mtoh_genes.xlsx")

write.xlsx(overlap_high_sig_exp_df_only,"overlap_high_sig_exp_df_only.xlsx")



