
library(TCGAbiolinks)


library(reshape)

library(DT)
library(SummarizedExperiment)

library(dplyr)

BiocManager::install("TCGAbiolinks")


query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query)


expdat <- GDCprepare(query = query)

count_matrix= as.data.frame(assay(expdat))


###############################################################change to gene_symbol

library(stringr)

count_matrix$Ensembl_ID <- rownames(count_matrix)

count_matrix$Ensembl_ID = str_sub(count_matrix$Ensembl_ID,1,15)


library(clusterProfiler)

library(org.Hs.eg.db)

library(limma)


keytypes(org.Hs.eg.db)

Ensembl_ID <- count_matrix$Ensembl_ID

gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL", "ENTREZID"), OrgDb="org.Hs.eg.db")

head(gene_symbol)

gene_symbol= gene_symbol[match(count_matrix$Ensembl_ID,gene_symbol$ENSEMBL),]

count_matrix$SYMBOL <- gene_symbol$SYMBOL

count_matrix <- count_matrix[,-(ncol(count_matrix)-1)]

count_matrix2 <- avereps(count_matrix, ID=count_matrix$SYMBOL)

count_matrix2 <- as.data.frame(count_matrix2)

count_matrix2 <- count_matrix2[!is.na(count_matrix2$SYMBOL),]

rownames(count_matrix2) <- count_matrix2$SYMBOL

count_matrix <- count_matrix2[,-ncol(count_matrix2)]


#################################################################selection samples
ls <- unlist(substr(colnames(count_matrix),14,15))

table(ls) 

group_list = ifelse(as.numeric(substr(colnames(count_matrix),14,15)) < 2, 'Tumor', 'normal')

table(group_list) 

count_matrix = na.omit(count_matrix)

count_matrix_tu <-  count_matrix[,group_list=="Tumor"]

count_matrix_n <-  count_matrix[,group_list=="normal"]

count_matrix3 <- cbind(count_matrix_n,count_matrix_tu)

save(list = c("count_matrix_tu","count_matrix3","count_data"),file = "PAAD-tumor.RData")

#########################################################################

count_data <- as.data.frame(lapply(count_matrix_tu, as.numeric))

rownames(count_data) <- rownames(count_matrix_tu)

colnames(count_data) <- colnames(count_matrix_tu)
################################################################
library(dplyr)

library(EnsDb.Hsapiens.v75)

edb <- EnsDb.Hsapiens.v75
genes_ensemble <- genes(edb)
gene_length<- as.data.frame(genes_ensemble) %>% dplyr::select(gene_id, gene_name, width)

############rid of 0 counts

PAADMatrix_clean<- count_data[rowSums(count_data) > 5, ]

gene_length_in_mat<- left_join(data.frame(gene_name = rownames(PAADMatrix_clean)), gene_length) %>% dplyr::filter(!is.na(width))

PAADMatrix_sub<- PAADMatrix_clean[rownames(PAADMatrix_clean) %in% gene_length_in_mat$gene_name, ]


all.equal(rownames(PAADMatrix_sub), gene_length_in_mat$gene_name)


sub_gene <- data.frame(gene_name = rownames(PAADMatrix_sub))

gene_length_in_mat2 <- left_join(sub_gene ,gene_length_in_mat,by = "gene_name")

gene_length_in_mat2 <-gene_length_in_mat2[-which(duplicated(gene_length_in_mat2$gene_name)),]

countToTpm <- function(counts, effLen)
{
  rate <- log(counts + 1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

PAAD_TPM<- apply(PAADMatrix_sub,2, countToTpm, effLen = gene_length_in_mat2$width)


#####################################################################



samples_df <- data.frame(sample = colnames(PAAD_TPM))

sam = str_sub(samples_df$sample,1,12)

PAAD_TPM2 <- PAAD_TPM

colnames(PAAD_TPM2) <- sam



#######################################################################
clinical <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")

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

intented_gene <- "NFE2"

PAAD_log_tpm <- log2(PAAD_TPM2+1)

NFE2_tpm <- data.frame(NFE2 = PAAD_log_tpm[intented_gene,])

library(tibble)

NFE2_tpm <- rownames_to_column(NFE2_tpm,var = 'Barcode')

NFE2_tpm$NFE2_group = ifelse(NFE2_tpm$NFE2 > median(NFE2_tpm$NFE2),"High","Low")

NFE2_surv <- merge(NFE2_tpm,survival_data,by = "Barcode")

NFE2_surv2 <- NFE2_surv[,-2]

colnames(NFE2_surv2)[2] <- "NFE2"

library(survival)


NFE2_SURVIVAL = survfit(Surv(OS.Time, OS)~NFE2,data=NFE2_surv2)


library(survminer)
library(survival)

library(ggpubr)

library(rstatix)

ggsurv_nfe2 <- ggsurvplot(NFE2_SURVIVAL,
           legend = "top",
           legend.title = "Group",
           palette = c("#F8766D","#619CFF"),
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

  

ggsurv_nfe2$plot <- ggsurv_nfe2$plot +
  ggplot2::annotate(
    "text",
    x = 450, y = 0.20,
    vjust = 1, hjust = 1,
    label = "HR = 1.5 \n p = 0.061",
    size = 3
  )

ggsurv_nfe2

library(ggsci)

ggsave(filename = "nfe2_survival_plot.png", plot = ggsurv_nfe2$plot, width = 12, height = 10, dpi = 300, units = "cm")

dev.off()

ggsurv_nfe2$plot <-ggsurv_nfe2$plot + 
  theme(legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 11, color = "black"))

ggsurv_nfe2$table <- ggsurv_nfe2$table +
  theme(axis.text.x = element_text(size = 11, color = "black" ),
        axis.text.y.left  = element_text(size = 11),
        axis.title = element_text(size = 11))



NFE2_surv2$Stages <- ifelse(NFE2_surv2$Stage == "Stage I"|NFE2_surv2$Stage == "Stage IA"| NFE2_surv2$Stage == "Stage IB","Stage I",ifelse(NFE2_surv2$Stage == "Stage IIA"|NFE2_surv2$Stage == "Stage IIB","Stage II",ifelse(NFE2_surv2$Stage == "Stage III","Stage III",ifelse(NFE2_surv2$Stage == "Stage IV","Stage IV",NA))) )



NFE2_surv$Stages <- ifelse(NFE2_surv$Stage == "Stage I"|NFE2_surv$Stage == "Stage IA"| NFE2_surv$Stage == "Stage IB","Stage I",ifelse(NFE2_surv$Stage == "Stage IIA"|NFE2_surv$Stage == "Stage IIB","Stage II",ifelse(NFE2_surv$Stage == "Stage III","Stage III",ifelse(NFE2_surv$Stage == "Stage IV","Stage IV",NA))) )

NFE2_surv<- NFE2_surv[,-7]

colnames(NFE2_surv)[8] <- "Stage"


NFE2_surv2$NFE2= as.factor(NFE2_surv2$NFE2)

NFE2_surv2$NFE2= relevel(NFE2_surv2$NFE2, ref = "Low")


by_cox = coxph(Surv(OS.Time, OS) ~ Age+Gender+NFE2,data=NFE2_surv2)

#by_cox2 = coxph(Surv(OS.Time, OS) ~ Age+Gender+NFE2+Stages,data=NFE2_surv)


ggforest(by_cox,data = NFE2_surv2, refLabel = 1, main = "Hazard ratio in PAAD")

#ggforest(by_cox2,data = NFE2_surv, refLabel = 1, main = "Hazard ratio in PAAD")


library(broom)


hr_table_nfe2 = data.frame(tidy(by_cox, exponentiate = T, conf.int = T)[, c("term", "estimate", "conf.low", "conf.high", "p.value")])




NFE2_surv3 = na.omit(NFE2_surv)



NFE2_surv3$Stages= as.factor(NFE2_surv3$Stages)



NFE2_surv3$Nstage= relevel(NFE2_surv3$Nstage, ref = "Stage I")


NFE2_surv3$Nstage <- gsub("N1b",replacement = "N1",NFE2_surv3$Nstage)


NFE2_surv3$Stages= as.factor(NFE2_surv3$Stages)



NFE2_surv3$Stages= relevel(NFE2_surv3$Stages, ref = "Stage I")



ggboxplot(NFE2_surv3  ,x = "Nstage",
          y = "NFE2", color = "Nstage", add = "jitter",ylab = " log2(TPM)+1 expression level") + geom_pwc(
            aes(group =Nstage),
            method = "wilcox_test", label = "p.signif") 



####################################################################################

cor_tpm <- data.frame(NFE2 = PAAD_log_tpm["NFE2",],
                      ITGA2B = PAAD_log_tpm["ITGA2B",])



scatter_nfe2_plt <- ggscatter(cor_tpm , x = "NFE2", y = "ITGA2B", xlab = "NFE2 expression level",
                              ylab = "ITGA2B expression level",
                              conf.int = TRUE,
                             ggtheme = theme_bw()# Add confidence interval
) 


scatter_nfe2_plt + stat_cor(method = "pearson",label.sep = "\n")

library(ggplot2)
library(envalysis)
library(ggpubr)
scatter_nfe2 <- ggscatter(cor_tpm, x = "NFE2", y = "ITGA2B", 
                          xlab = "NFE2 expression level",
                          ylab = "ITGA2B expression level",
                          conf.int = TRUE,
                          ggtheme = theme_bw())+ theme_publish(base_size = 11, base_linewidth = 0.7) + stat_cor(method = "pearson",label.sep = "\n", size = 4)# Add confidence interval

ggsave(filename = "PAAD_NFE2_ITGA2B_COR.png", plot = scatter_nfe2, width = 12, height = 10, dpi = 300, units = "cm")
#####################################################XCELL

library(xCell)

library(tibble)

scores = xCellAnalysis(PAAD_TPM2)


selects_df = as.data.frame(scores)

selects_df = rownames_to_column(selects_df,"Gene")

selects_i = subset(selects_df, selects_df$Gene == "Platelets")


rownames(selects_i) <- NULL

selects_i <- column_to_rownames(selects_i,"Gene")

selects_plt<- data.frame(t(selects_i))

selects_plt  <- rownames_to_column(selects_plt ,var = 'Barcode')

NFE2_plt <- merge(NFE2_tpm,selects_plt,by = "Barcode")

NFE2_plt_high <- subset(NFE2_plt,NFE2_plt$NFE2_group == "High")

library(ggpubr)


scatter_nfe2_plt <- ggscatter(NFE2_plt_high, x = "Platelets", y = "NFE2", xlab = "Platelets infiltration level",
                ylab = "NFE2 expression level",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) 

scatter_nfe2_plt + stat_cor(method = "pearson")

###################################################SSGSEA_SCORE_COR


Platelet = c("ITGA2B","GP1BA") 


Gene_set= list(Platelet)

names(Gene_set) = c("Platelet")


library(GSVA)

library(dplyr)

ssgsea = gsva(PAAD_TPM2, Gene_set,method = "ssgsea")


ssgsea_df = as.data.frame(ssgsea)

ssgsea_df = as.data.frame(t(ssgsea_df))

ssgsea_df <- rownames_to_column(ssgsea_df  ,var = 'Barcode')

NFE2_plt2 <- merge(NFE2_tpm,ssgsea_df,by = "Barcode")


NFE2_plt_high2 <- subset(NFE2_plt2,NFE2_plt2$NFE2_group == "High")

library(ggpubr)


scatter_nfe2_plt2 <- ggscatter(NFE2_plt2, x = "Platelet", y = "NFE2", xlab = "Platelets infiltration level",
                              ylab = "NFE2 expression level",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE # Add confidence interval
) 

scatter_nfe2_plt2 + stat_cor(method = "pearson")

#################################################################
AST = c("HBA1","HBA2","HBB","GATA1","ITGA2B","GP1BA")
AST_tpm <- as.data.frame(PAAD_log_tpm[AST,])

AST_tpm_t <- data.frame(t(AST_tpm))

AST_tpm_t  <- rownames_to_column(AST_tpm_t  ,var = 'Barcode')


NFE2_AST <- merge(NFE2_tpm,AST_tpm_t ,by = "Barcode")

colnames(NFE2_AST)[3] = "NFE2_Group"

library(tidyr)

dat_tpm <- NFE2_AST %>%
  gather(key = gene,value = Proportion,-Barcode,-NFE2_Group)



library(ggplot2)

library(ggpubr)

library(rstatix)

ggboxplot(dat_tpm ,x = "gene",
          y = "Proportion", color = "NFE2_Group", add = "jitter",ylab = " log2(TPM)+1 expression level") + geom_pwc(
            aes(group = NFE2_Group),
            method = "wilcox_test", label = "p.signif") 


#######################################################
NFE2_tpm2 <- NFE2_tpm

NFE2_tpm2 = NFE2_tpm2[order(NFE2_tpm2$NFE2_group,decreasing = TRUE),]

PAAD_TPM_LM <- PAAD_TPM2

PAAD_TPM_LM <- PAAD_TPM_LM[,NFE2_tpm2$Barcode]

group_list= NFE2_tpm2$NFE2_group

group_list <- factor(group_list,levels = c("Low","High"),ordered = F)

library(limma)

exprSet=normalizeBetweenArrays(PAAD_TPM_LM)

exprSet <- log2(exprSet+1)

design <- model.matrix(~0+group_list)

colnames(design) <- c("Low","High")

fit <- lmFit(exprSet  , design = design)


constrasts <- makeContrasts(High-Low,levels = design)

fitcontrasts <- contrasts.fit(fit,constrasts)

fitcontrasts2 <- eBayes(fitcontrasts)

dd_con <- topTable(fitcontrasts2 , n=Inf, sort="none")


library(tibble)

dd_con_df <- rownames_to_column(dd_con,"Gene")


NFE2_H_L <- dd_con

cut_pvalue <- 0.05

cut_lfc <- 1


significant_NFE2_H_L_EXP <- NFE2_H_L[which(NFE2_H_L$adj.P.Val< cut_pvalue & (NFE2_H_L$logFC<(-cut_lfc) | NFE2_H_L$logFC>cut_lfc)),]

##################################################
pathways.hallmark.ROS_meta <- read.table("GOBP_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS.v2022.1.Hs.tsv", sep = ",",header = T)

pathways.hallmark.ROS_res <- read.table("GOBP_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES.v2022.1.Hs.tsv", sep = ",",header = T)

pathways.hallmark.ROS_meta <- data.frame(gene.set = pathways.hallmark.ROS_meta[-c(1:256,262),])

pathways.hallmark.ROS_res <- data.frame(gene.set = pathways.hallmark.ROS_res[-c(1:229),])

ROS_metabolic_process <- pathways.hallmark.ROS_meta$gene.set

Response_to_ROS <- pathways.hallmark.ROS_res$gene.set


Gene_set_2 = list(ROS_metabolic_process,Response_to_ROS)

Gene_set_2

names(Gene_set_2 ) = c("ROS_metabolic_process","Response_to_ROS")

PAAD_TPM_LM

PAAD_TPM_LM_log <- log2(PAAD_TPM_LM+1)

library(GSVA)

ssgsea_ros = gsva(PAAD_TPM_LM_log, Gene_set_2,method = "ssgsea")


group_list2= NFE2_tpm2$NFE2_group


colData <- data.frame(row.names = NFE2_tpm2$Barcode,Group = NFE2_tpm2$NFE2_group)


library(pheatmap)

pheatmap(ssgsea_ros,
         annotation_col = colData,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         fontsize = 8,
         scale = "row")



ssgsea_ros = as.data.frame(ssgsea_ros)

ssgsea_ros = as.data.frame(t(ssgsea_ros))

ssgsea_ros = rownames_to_column(ssgsea_ros,"Sample")

colData_ros <- rownames_to_column(colData,"Sample")

ssgsea_ros = merge(ssgsea_ros,colData_ros ,by = "Sample")

ssgsea_ros = column_to_rownames(ssgsea_ros,"Sample")


library(ggplot2)

library(ggpubr)



ggboxplot(ssgsea_ros,x = "Group",
          y = "ROS_metabolic_process", color = "Group", add = "jitter",ylab = "ROS_metabolic_process",palette = "jco") +  stat_compare_means(method = "t.test")

