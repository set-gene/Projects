library(TCGAbiolinks)


library(reshape)

library(DT)
library(SummarizedExperiment)

library(dplyr)

query <- GDCquery(project = "TCGA-KIRC", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts")

GDCdownload(query, method = "api")


expdat <- GDCprepare(query = query)


mrna = assay(expdat)

tpm <- assay(expdat,4) %>% as.data.frame()

###############################################################change to gene_symbol

library(stringr)

tpm$Ensembl_ID <- rownames(tpm)

tpm$Ensembl_ID = str_sub(tpm$Ensembl_ID,1,15)


library(clusterProfiler)

library(org.Hs.eg.db)

library(limma)


keytypes(org.Hs.eg.db)

Ensembl_ID <- tpm$Ensembl_ID

gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")

head(gene_symbol)

gene_symbol= gene_symbol[match(tpm$Ensembl_ID,gene_symbol$ENSEMBL),]

tpm$SYMBOL <- gene_symbol$SYMBOL

tpm <- tpm[,-(ncol(tpm)-1)]

tpm2 <- avereps(tpm, ID=tpm$SYMBOL)

tpm2 <- as.data.frame(tpm2)

tpm2 <- tpm2[!is.na(tpm2$SYMBOL),]

rownames(tpm2) <- tpm2$SYMBOL

count_matrix <- tpm2[,-ncol(tpm2)]


#################################################################selection samples
ls <- unlist(substr(colnames(count_matrix),14,15))

table(ls) 

group_list = ifelse(as.numeric(substr(colnames(count_matrix),14,15)) < 2, 'Tumor', 'normal')

table(group_list) 

count_matrix = na.omit(count_matrix)

count_matrix_tu <-  count_matrix[,group_list=="Tumor"]

count_matrix_n <-  count_matrix[,group_list=="normal"]

count_matrix3 <- cbind(count_matrix_n,count_matrix_tu)

save(list = c("count_matrix_tu","count_matrix3","count_data"),file = "KIRC_TPM-tumor.RData")

#########################################################################

count_data <- as.data.frame(lapply(count_matrix_tu, as.numeric))

rownames(count_data) <- rownames(count_matrix_tu)

colnames(count_data) <- colnames(count_matrix_tu)
################################################################


samples_df <- data.frame(sample = colnames(count_data))

sam = str_sub(samples_df$sample,1,12)

count_data2 <- count_data

colnames(count_data2) <- sam



#######################################################################
clinical <- GDCquery_clinic(project = "TCGA-KIRC", type = "clinical")

clinical_trait <- clinical  %>%
  dplyr::select(bcr_patient_barcode,gender,vital_status,days_to_last_follow_up ,days_to_death,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_n,ajcc_pathologic_m,age_at_index) %>%distinct( bcr_patient_barcode, .keep_all = TRUE)  


dead_patient <- clinical_trait  %>%
  dplyr::filter(vital_status == 'Dead') %>%
  dplyr::select(-days_to_last_follow_up) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_death='OS.Time',
                    ajcc_pathologic_stage = 'Stage',
                    age_at_index = 'Age',
                    ajcc_pathologic_t = "Tstage",
                    ajcc_pathologic_n = "Nstage",
                    ajcc_pathologic_m = "Mstage")) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))


alive_patient <- clinical_trait  %>%
  dplyr::filter(vital_status == 'Alive') %>%
  dplyr::select(-days_to_death) %>%
  reshape::rename(c(bcr_patient_barcode = 'Barcode',
                    gender = 'Gender',
                    vital_status = 'OS',
                    days_to_last_follow_up='OS.Time',
                    ajcc_pathologic_stage = 'Stage',
                    age_at_index = 'Age',
                    ajcc_pathologic_t = "Tstage",
                    ajcc_pathologic_n = "Nstage",
                    ajcc_pathologic_m = "Mstage")) %>%
  mutate(OS=ifelse(OS=='Dead',1,0))

survival_data <- rbind(dead_patient,alive_patient)

intented_gene <- "SORBS1"

SORBS1_count <- data.frame(SORBS1 = t(count_data2[intented_gene,]))

library(tibble)

SORBS1_count <- rownames_to_column(SORBS1_count,var = 'Barcode')

SORBS1_count$Barcode <- gsub("[:.:]","-",SORBS1_count$Barcode)

SORBS1_count$SORBS1_group = ifelse(SORBS1_count$SORBS1 > median(SORBS1_count$SORBS1),"High","Low")

SORBS1_surv <- merge(SORBS1_count,survival_data,by = "Barcode")

SORBS1_surv2 <- SORBS1_surv[,-2]

colnames(SORBS1_surv2)[2] <- "SORBS1"

library(survival)


SORBS1_SURVIVAL = survfit(Surv(OS.Time, OS)~SORBS1,data=SORBS1_surv2)


library(survminer)
library(survival)

library(ggpubr)

library(rstatix)

ggsurv_SORBS1 <- ggsurvplot(SORBS1_SURVIVAL,
           legend = "top",
           legend.title = "Group",
           palette = c("purple","green"),
           xlab="Days",
           ylab="Survival probability",
           font.x = 11,
           font.y = 11,
           pval = TRUE,
           pval.method = TRUE,
           pval.size = 3,
           pval.coord = c(50,.24), 
           risk.table = FALSE,
           ggtheme = theme_classic2(base_size=11))



ggsurv_SORBS1$plot <- ggsurv_SORBS1$plot +
  ggplot2::annotate(
    "text",
    x = 600, y = 0.20,
    vjust = 1, hjust = 1,
    label = "HR = 1.4 \n p = 0.046",
    size = 3
  )

ggsave(filename = "SORBS1_survival_plot.png", plot = ggsurv_SORBS1$plot, width = 12, height = 10, dpi = 300, units = "cm")

dev.off()

table(survival_data$Stage)



by_cox = coxph(Surv(OS.Time, OS) ~ Age+Stage+SORBS1,data=SORBS1_surv2)


ggforest(by_cox,data = SORBS1_surv2, refLabel = 1, main = "Hazard ratio in KIRC")

####################################################################################

cor_count<- data.frame(SORBS1= t(log2(count_data2["SORBS1",]+1)),
                      GPX4 = t(log2(count_data2["GPX4",]+1)))

library(ggplot2)
library(envalysis)
library(ggpubr)
scatter_SORBS1<- ggscatter(cor_count , x = "SORBS1", y = "GPX4", xlab = "SORBS1 expression level",
                              ylab = "GPX4 expression level",
                              conf.int = TRUE,
                             ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

ggsave(filename = "KIRC_SORBS1_GPX4_COR.png", plot = scatter_SORBS1, width = 12, height = 10, dpi = 300, units = "cm")
#######################################################################

cor_count2<- data.frame(SORBS1= t(log2(count_data2["SORBS1",]+1)),
                       FTH1 = t(log2(count_data2["FTH1",]+1)))

library(ggplot2)
library(envalysis)
library(ggpubr)
scatter_FTH1<- ggscatter(cor_count2 , x = "SORBS1", y = "FTH1", xlab = "SORBS1 expression level",
                           ylab = "FTH1 expression level",
                           conf.int = TRUE,
                           ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

ggsave(filename = "KIRC_SORBS1_FTH1_COR.png", plot = scatter_SORBS1, width = 12, height = 10, dpi = 300, units = "cm")
######################################################################

cor_count3<- data.frame(SORBS1= t(log2(count_data2["SORBS1",]+1)),
                        TFRC = t(log2(count_data2["TFRC",]+1)))

scatter_TFRC<- ggscatter(cor_count3 , x = "SORBS1", y = "TFRC", xlab = "SORBS1 expression level",
                         ylab = "TFRC expression level",
                         conf.int = TRUE,
                         ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval
######################################################################
cor_count4<- data.frame(SORBS1= t(log2(count_data2["SORBS1",]+1)),
                        PTGS2 = t(log2(count_data2["PTGS2",]+1)))

scatter_PTGS2<- ggscatter(cor_count4 , x = "SORBS1", y = "PTGS2", xlab = "SORBS1 expression level",
                         ylab = "PTGS2 expression level",
                         conf.int = TRUE,
                         ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Addconfidence interval
########################################################################

cor_count<- data.frame(SORBS1= t(log2(count_data2["SORBS1",]+1)),
                       FGF2 = t(log2(count_data2["FGF2",]+1)))

library(ggplot2)
library(envalysis)
library(ggpubr)
scatter_FGF2<- ggscatter(cor_count , x = "SORBS1", y = "FGF2", xlab = "SORBS1 expression level",
                           ylab = "FGF2 expression level",
                           conf.int = TRUE,
                           ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

ggsave(filename = "KIRC_SORBS1_GPX4_COR.png", plot = scatter_SORBS1, width = 12, height = 10, dpi = 300, units = "cm")

#######################################################################
log_tpm <- log2(count_data2 + 1)

NN_cor  = as.matrix(t(log_tpm))

library("Hmisc")

ccor_select <- rcorr(NN_cor, type ="pearson")

ccor_select_g$P

cor_df = data.frame(ccor_EMT_MT$r)

cor_df = rownames_to_column(cor_df,"Gene")


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

filter_cor = flattenCorrMatrix(ccor_select_g$r, ccor_select_g$P)

write.table(filter_cor,"cor_CCN4_MYH11.txt", sep = "\t", quote = FALSE)



######################################################################

mrna <- as.data.frame(mrna)

###############################################################change to gene_symbol

library(stringr)

mrna$Ensembl_ID <- rownames(mrna)

mrna$Ensembl_ID = str_sub(mrna$Ensembl_ID,1,15)


library(clusterProfiler)

library(org.Hs.eg.db)

library(limma)


keytypes(org.Hs.eg.db)

Ensembl_ID <- mrna$Ensembl_ID

gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")

head(gene_symbol)

gene_symbol= gene_symbol[match(mrna$Ensembl_ID,gene_symbol$ENSEMBL),]

mrna$SYMBOL <- gene_symbol$SYMBOL

mrna <- mrna[,-(ncol(mrna)-1)]

mrna2 <- avereps(mrna, ID=mrna$SYMBOL)

mrna2 <- as.data.frame(mrna2)

mrna2 <- mrna2[!is.na(mrna2$SYMBOL),]

rownames(mrna2) <- mrna2$SYMBOL

count_mrna <- mrna2[,-ncol(mrna2)]

#################################################################selection samples
ls2 <- unlist(substr(colnames(count_mrna),14,15))

table(ls2) 

group_list2 = ifelse(as.numeric(substr(colnames(count_mrna),14,15)) < 2, 'Tumor', 'normal')

table(group_list2) 

count_mrna = na.omit(count_mrna)

count_mrna_tu <-  count_mrna[,group_list=="Tumor"]

count_mrna_n <-  count_mrna[,group_list=="normal"]

count_mrna3 <- cbind(count_mrna_n,count_mrna_tu)




#########################################################################

count_data_mrna <- as.data.frame(lapply(count_mrna_tu, as.numeric))

rownames(count_data_mrna) <- rownames(count_mrna_tu)

colnames(count_data_mrna) <- colnames(count_mrna_tu)
################################################################

samples_mrna_df <- data.frame(sample = colnames(count_data_mrna))

sam2 = str_sub(samples_mrna_df$sample,1,12)

count_data2_mrna <- count_data_mrna

colnames(count_data2_mrna) <- sam2

SORBS1_count

SORBS1_Group <-SORBS1_count[,c(1,3)]

library(tibble)

SORBS1_Group <- column_to_rownames(SORBS1_Group,"Barcode")

colnames(SORBS1_Group)[1] <- "Group"

SORBS1_Group$Group <-as.factor(SORBS1_Group$Group)


colnames(count_data2_mrna) <- rownames(SORBS1_Group)


library(DESeq2)

dds <-DESeqDataSetFromMatrix(count_data2_mrna,SORBS1_Group,design = ~ Group )


dds <- DESeq(dds)

res = DESeq2::results(dds  ,c("Group", "Low","High"))

res

res2 = res[order(res$padj),]

#########################################################################
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



#########################################################################

library(fgsea)

library(dplyr)

library(tibble)

res_df = data.frame(res)

res_df = rownames_to_column(res_df, "SYMBOL")  

res3 <- res_df %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))


res3$SYMBOL

ranks <- deframe(res3)


kegg = gmtPathways("c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")

kegg %>% 
  head() %>% 
  lapply(head)


kegg_Res <- fgsea(pathways=kegg, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)

GO = gmtPathways("c5.go.bp.v2023.2.Hs.symbols.gmt")

GO_Res <- fgsea(pathways=GO, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)


go_ResT <- GO_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

go_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

library(ggplot2)


##############################################################
go_ResT_up = subset(x = go_ResT, go_ResT$NES > 1)


go_ResT_dn = subset(x = go_ResT, go_ResT$NES < -1)


go_ResT_dn <- go_ResT_dn[order(go_ResT_dn$padj),]


topGOPathwaysUp <- GO_Res[ES > 0][head(order(pval), n=20), pathway]

topGOPathwaysDown <- GO_Res[ES < 0][head(order(pval), n=20), pathway]

topPathways <- c(topGOPathwaysUp , rev(topGOPathwaysDown))




go_ResT_up_FOR = subset(go_ResT_up,go_ResT_up$pathway %in% topGOPathwaysUp)


go_ResT_dn_FOR = subset(go_ResT_dn,go_ResT_dn$pathway %in% topGOPathwaysDown)

go_ResT_all <- rbind(go_ResT_up_FOR,go_ResT_dn_FOR)
##############################################################################




ggplot(go_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO BP pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))




plotGseaTable(GO[topPathways], ranks, GO_Res, 
              gseaParam=0.5)


GO[topPathways]

kegg_ResT <- kegg_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

kegg_ResT %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

kegg_ResT_up = subset(x = kegg_ResT, kegg_ResT$NES > 1)


kegg_ResT_dn = subset(x = kegg_ResT, kegg_ResT$NES < -1)

kegg_ResT_dn <- kegg_ResT_dn[order(kegg_ResT_dn$padj),]



topkeggPathwaysUp <- kegg_Res[ES > 0][head(order(pval), n=20), pathway]

topkeggPathwaysDown <- kegg_Res[ES < 0][head(order(pval), n=20), pathway]

topkeggPathways <- c(topkeggPathwaysUp , rev(topkeggPathwaysDown))


kegg_ResT_up_FOR = subset(kegg_ResT_up,kegg_ResT_up$pathway %in% topkeggPathwaysUp)


kegg_ResT_dn_FOR = subset(kegg_ResT_dn,kegg_ResT_dn$pathway %in% topkeggPathwaysDown)

kegg_ResT_all <- rbind(kegg_ResT_up_FOR,kegg_ResT_dn_FOR)
##############################################################################

ggplot(kegg_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))


#plotGseaTable(kegg[topkeggPathways], ranks, kegg_Res, 
              #gseaParam=0.5) + cowplot::draw_text(list(size=5))

library(dplyr)

hmrk = gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

hmrk %>% 
  head() %>% 
  lapply(head)


hmrk_Res <- fgsea(pathways=hmrk, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)


hmrk_ResT <- hmrk_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

hmrk_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

##############################################################
hmrk_ResT_up = subset(x = hmrk_ResT, hmrk_ResT$NES > 1)


hmrk_ResT_dn = subset(x = hmrk_ResT, hmrk_ResT$NES < -1)


hmrk_ResT_dn <- hmrk_ResT_dn[order(hmrk_ResT_dn$padj),]

tophmrkPathwaysUp <- hmrk_Res[ES > 0][head(order(pval), n=20), pathway]

tophmrkPathwaysDown <- hmrk_Res[ES < 0][head(order(pval), n=20), pathway]

topPathways <- c(tophmrkPathwaysUp , rev(tophmrkPathwaysDown))


hmrk_ResT_up_FOR = subset(hmrk_ResT_up,hmrk_ResT_up$pathway %in% tophmrkPathwaysUp)


hmrk_ResT_dn_FOR = subset(hmrk_ResT_dn,hmrk_ResT_dn$pathway %in% tophmrkPathwaysDown)

hmrk_ResT_all <- rbind(hmrk_ResT_up_FOR,hmrk_ResT_dn_FOR)


library(openxlsx)

write.xlsx(hmrk_ResT_all,"hm_dn.xlsx")

##############################################################################




ggplot(hmrk_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))

##############################################################################
log_tpm 

stage_log_tpm <- log_tpm[,clinical_stage_exp$submitter_id]


hmrk
############################################################################
log_tpm <- log2(count_data2 + 1)

select_gene = c("SORBS1","EPAS1","FGF2","SLC2A4")

rm(list = c("clinical_select_count"))


select_count <- as.data.frame(t(log_tpm[select_gene,]))

library(tibble)

select_count <- rownames_to_column(select_count,"submitter_id")

select_count$submitter_id <- gsub("[:.:]","-",select_count$submitter_id)

clinical_stage <- clinical[,c("submitter_id","ajcc_pathologic_stage")]


clinical_select_count <- merge(select_count,clinical_stage,"submitter_id")

colnames(clinical_select_count)[6] <- "Stage"

clinical_select_count<- clinical_select_count[!is.na(clinical_select_count$Stage),]

clinical_select_count$Stage <- as.factor(clinical_select_count$Stage)

rownames(clinical_select_count) <- NULL

clinical_select_count <- column_to_rownames(clinical_select_count,"submitter_id")


library(ggplot2)

library(ggpubr)

library(rstatix)

ggboxplot(clinical_select_count,x = "Stage",
          y = "SORBS1", color = "Stage", add = "jitter",ylab = "SORBS1 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 
##########################################################################################################

ggboxplot(clinical_select_count,x = "Stage",
          y = "EPAS1", color = "Stage", add = "jitter",ylab = "EPAS1 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 

##################################################################################

ggboxplot(clinical_select_count,x = "Stage",
          y = "KAT7", color = "Stage", add = "jitter",ylab = "KAT7 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 

####################################################################################

ggboxplot(clinical_select_count,x = "Stage",
          y = "SLC2A4", color = "Stage", add = "jitter",ylab = "SLC2A4(GLUT4) expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


##################################################################################

ggboxplot(clinical_select_count,x = "Stage",
          y = "FGF2", color = "Stage", add = "jitter",ylab = "FGF2 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


###################################################################################
stage_clinical_select_count <- subset(clinical_select_count,clinical_select_count$Stage == "Stage IV")

scatter_stage_cor<- ggscatter(stage_clinical_select_count  , x = "SORBS1", y = "FGF2", xlab = "SORBS1 expression level",
                         ylab = "FGF2 expression level",
                         conf.int = TRUE,
                         ggtheme = theme_bw()# Add confidence interval
) + theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

##########################################################################################################3

library(tibble)

SORBS1_Group_df <- rownames_to_column(SORBS1_Group,"submitter_id")

SORBS1_Group_df$submitter_id <- as.character(SORBS1_Group_df$submitter_id )


ssgsea_sorbs1_gr = gsva(as.matrix(count_data2), Gene_set2,method = "ssgsea")

ssgsea_df3 <- as.data.frame(t(ssgsea_sorbs1_gr))

ssgsea_df3 <- rownames_to_column(ssgsea_df3 ,"submitter_id")

ssgsea_df3$submitter_id <- gsub("[:.:]","-",ssgsea_df3$submitter_id)


ssgsea_df3 <- merge(ssgsea_df3,SORBS1_Group_df,"submitter_id")


ggboxplot(ssgsea_df3 ,x = "Group",
          y = "G6PMetabolicProcess", color = "Group", add = "jitter",ylab = "G6PMetabolicProcess ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Group),
            method = "t_test", label = "p.signif") 


ggboxplot(ssgsea_df3 ,x = "Group",
          y = "CellCycle", color = "Group", add = "jitter",ylab = "CellCycle ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Group),
            method = "t_test", label = "p.signif") 


###########################################################################################################

library(dplyr)

library(ggplot2)

library(GSVA)

clinical_stage <- clinical_stage[!is.na(clinical_stage$ajcc_pathologic_stage),]

clinical_stage_exp<- subset(clinical_stage,clinical_stage$submitter_id %in% colnames(count_data2))

colnames(clinical_stage_exp)[2] <- "Stage"

stage_count <- as.data.frame(t(count_data2))

stage_count <- rownames_to_column( stage_count , "submitter_id")


stage_count$submitter_id <- gsub("[:.:]","-",stage_count$submitter_id)


stage_count <- merge(clinical_stage_exp,stage_count,"submitter_id")

rownames(stage_count) <- NULL

stage_count <- stage_count[,-2]


stage_count <- column_to_rownames(stage_count,"submitter_id")

stage_count2 <- as.data.frame(t( stage_count))

#stage_count2 <- log2(stage_count2+1)

HYPOXIA <- hmrk$HALLMARK_HYPOXIA

#EMT <- hmrk$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

EpithelialMesenchymalTransition <- c("ZEB2","ZEB1","FN1","TWIST2","VIM","SNAI2","TWIST1","CDH2","SNAI1","MUC1",
                                     "TJP3","CLDN4","CLDN7","CDH1")

CellCycle <- kegg$KEGG_CELL_CYCLE

G2MPhase <- GO$GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION

G0G1Phase <-GO$GOBP_G0_TO_G1_TRANSITION

G2MChekpoint <- hmrk$HALLMARK_G2M_CHECKPOINT

G1SPhase <-GO$GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION

PhaseTransition <- GO$GOBP_CELL_CYCLE_PHASE_TRANSITION

Ferroptosis <- GO$GOBP_FERROPTOSIS

GlucoseMetabolism <-GO$GOBP_GLUCOSE_METABOLIC_PROCESS

G6PMetabolicProcess <- GO$GOBP_GLUCOSE_6_PHOSPHATE_METABOLIC_PROCESS

FattyAcidMetabolsim <- GO$GOBP_FATTY_ACID_METABOLIC_PROCESS

FattyAcidSynthesis <- GO$GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS

LipidStorage <- GO$GOBP_LIPID_STORAGE

Gene_set = list(HYPOXIA)

GlucoseMetabolism

G6PMetabolicProcess




Gene_set2 = list(CellCycle,Ferroptosis,GlucoseMetabolism,G6PMetabolicProcess,G1SPhase,G2MPhase,PhaseTransition,
                 G2MChekpoint,G0G1Phase,EMT,EpithelialMesenchymalTransition,FattyAcidMetabolsim,FattyAcidSynthesis,LipidStorage,FattyAcidTransport)



names(Gene_set) = c("HYPOXIA")


names(Gene_set2) = c("CellCycle","Ferroptosis","GlucoseMetabolism","G6PMetabolicProcess","G1SPhase","G2MPhase","PhaseTransition",
                     "G2MChekpoint","G0G1Phase","EMT","EpithelialMesenchymalTransition","FattyAcidMetabolsim","FattyAcidSynthesis",
                     "LipidStorage","FattyAcidTransport")

library(GSVA)

ssgsea2 = gsva(as.matrix(stage_count2), Gene_set2,method = "ssgsea")

ssgsea_df2 <- as.data.frame(t(ssgsea2))

ssgsea_df2 <- rownames_to_column(ssgsea_df2 ,"submitter_id")

ssgsea_df2 <- merge(ssgsea_df2,clinical_stage_exp,"submitter_id")

library(tibble)

ssgsea_df2$Stage <- as.factor(ssgsea_df2$Stage)

ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "G2MChekpoint", color = "Stage", add = "jitter",ylab = "G2MChekpoint ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "PhaseTransition", color = "Stage", add = "jitter",ylab = "PhaseTransition ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "G6PMetabolicProcess", color = "Stage", add = "jitter",ylab = "G6PMetabolicProcess ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "CellCycle", color = "Stage", add = "jitter",ylab = "CellCycle ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "GlucoseMetabolism", color = "Stage", add = "jitter",ylab = "GlucoseMetabolism ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "Ferroptosis", color = "Stage", add = "jitter",ylab = "Ferroptosis ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "G0G1Phase", color = "Stage", add = "jitter",ylab = "G0 and G1 Phase Transition ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "G1SPhase", color = "Stage", add = "jitter",ylab = "G1 and S Phase ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "G2MPhase", color = "Stage", add = "jitter",ylab = "G2 and M Phase ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "EMT", color = "Stage", add = "jitter",ylab = "EMT ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 

ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "EpithelialMesenchymalTransition", color = "Stage", add = "jitter",ylab = "EMT ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


library(ggpubr)

ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "FattyAcidMetabolsim", color = "Stage", add = "jitter",ylab = "FattyAcidMetabolsim ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 




ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "FattyAcidSynthesis", color = "Stage", add = "jitter",ylab = "FattyAcidSynthesis ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "LipidStorage", color = "Stage", add = "jitter",ylab = "LipidStorage ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 


ggboxplot(ssgsea_df2 ,x = "Stage",
          y = "FattyAcidTransport", color = "Stage", add = "jitter",ylab = "FattyAcidTransport ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 

FattyAcidMetabolsim <- GO$GOBP_FATTY_ACID_METABOLIC_PROCESS

FattyAcidSynthesis <- GO$GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS

LipidStorage <- GO$GOBP_LIPID_STORAGE

FattyAcidTransport <- GO$GOBP_FATTY_ACID_TRANSPORT


colData <- data.frame(row.names = NFE2_tpm2$Barcode,Group = NFE2_tpm2$NFE2_group)


library(pheatmap)

pheatmap(ssgsea2,
         annotation_col = clinical_stage_exp2,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 8,
         scale = "row")
library(ggplot2)

library(ggstats)

library(tibble)

library(GSVA)

library(dplyr)

#expr_s = bind_cols(exprSet_s1,exprSet_s2,exprSet_s3,exprSet_s4)

ssgsea = gsva(as.matrix(stage_count2), Gene_set,method = "ssgsea")

ssgsea_df <- as.data.frame(t(ssgsea))

ssgsea_df <- rownames_to_column(ssgsea_df ,"submitter_id")

ssgsea_df <- merge(ssgsea_df,clinical_stage_exp,"submitter_id")


ssgsea_df$Stage <- as.factor(ssgsea_df$Stage)

ssgsea_df <- column_to_rownames(ssgsea_df,"submitter_id")




library(ggplot2)

library(ggpubr)

library(rstatix)

ggboxplot(ssgsea_df ,x = "Stage",
          y = "HYPOXIA", color = "Stage", add = "jitter",ylab = "HYPOXIA ssgsea score")+ scale_fill_manual()+ geom_pwc(
            aes(Stage),
            method = "t_test", label = "p.signif") 



library(dplyr)
library(tibble)


clinical_select_count_df <- rownames_to_column(clinical_select_count,"submitter_id")

SORBS1_count_clinical <- clinical_select_count_df%>%
  select(submitter_id,SORBS1)


select_ssgsea <- merge(SORBS1_count_clinical,ssgsea_df2)

select_ssgsea_stage <- subset(select_ssgsea,select_ssgsea$Stage == "Stage IV")
  
select_ssgsea_cor <- ggscatter(select_ssgsea_stage, x = "SORBS1", y = "EMT", xlab = "SORBS1 expression level",
                              ylab = "EMT ssgsea score",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) 

select_ssgsea_cor + stat_cor(method = "pearson")

###############################################################
select_count_hypoxia <- stage_count2[HYPOXIA ,]

rownames(clinical_stage_exp) <- NULL

clinical_stage_exp2 <- column_to_rownames(clinical_stage_exp,"submitter_id")

clinical_stage_exp2$Stage <- as.factor(clinical_stage_exp2$Stage)

library(pheatmap)

pheatmap(select_count_hypoxia,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = clinical_stage_exp2,
         fontsize = 8)
