library(oligo)

library(pd.hugene.1.0.st.v1)

library(hugene10sttranscriptcluster.db)

library(dplyr)

library(data.table)

raw_data_dir <- getwd()

sdrf_location <- file.path(raw_data_dir, "E-MTAB-3267.sdrf.txt")

SDRF <- read.delim(sdrf_location)

T_SDRF <- subset(SDRF,SDRF$Characteristics.disease.=="Tumor")

T_SDRF=T_SDRF[ !(T_SDRF['Characteristics.sunitinib.response.']=='CLINICAL BENEFIT'),]

T_SDRF$Sunitinib_response <- ifelse(T_SDRF$Characteristics.sunitinib.response.=="PR","Response","non_Response")

table(T_SDRF$Characteristics.sunitinib.response.)


T_SDRF$Sunitinib_response <- as.factor(T_SDRF$Sunitinib_response)

rownames(T_SDRF) <- T_SDRF$Array.Data.File

T_SDRF  <- AnnotatedDataFrame(T_SDRF)

raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, T_SDRF$Array.Data.File), verbose = FALSE, phenoData = T_SDRF)

library(limma)

targets <- as.data.frame(strsplit2(rownames(raw_data@phenoData),"_")[,2:3])

colnames(targets) <- c("ID","Response")


targets$Response <-  T_SDRF@data$Sunitinib_response
  
rownames(targets) <- rownames(raw_data@phenoData)


Biobase::pData(raw_data) <- targets

eset <- oligo::rma(raw_data)

anno <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                              keys = (featureNames(eset)),
                              columns = c("SYMBOL", "GENENAME", "ENTREZID"), keytype = "PROBEID")

anno <- subset(anno, !is.na(SYMBOL))

anno_grouped <- group_by(anno, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(ENTREZID))

head(anno_summarized)

anno_filtered <- filter(anno_summarized, no_of_matches > 1)

head(anno_filtered)

ids_to_exlude <- (featureNames(eset) %in% anno_filtered$PROBEID)

table(ids_to_exlude)

eset_final <- subset(eset, !ids_to_exlude)

validObject(eset_final)

fData(eset_final)$PROBEID <- rownames(fData(eset_final))
fData(eset_final) <- dplyr::left_join(fData(eset_final), anno)

rownames(fData(eset_final)) <- fData(eset_final)$PROBEID 
validObject(eset_final)

keep <- !is.na(fData(eset_final)$ENTREZID)
eset_final <- eset_final[keep,]
validObject(eset_final)

norm_expr = log2(Biobase::exprs(eset_final))

rownames(norm_expr) <- pheno$SYMBOL

pheno = fData(eset_final)

select_gene = c("CD44","CDK1","CDK2","KDM5B","KLF4","MKI67","PCNA","SORBS1")


select_count <- as.data.frame(t(norm_expr[select_gene,]))

library(tibble)

select_count_t <- rownames_to_column(select_count,"sample_ID")

target_t <- rownames_to_column(targets,"sample_ID")

suni_select_count <- merge(select_count_t,target_t,"sample_ID")

suni_select_count <- column_to_rownames(suni_select_count ,"sample_ID")

library(ggplot2)

library(ggpubr)


ggboxplot(suni_select_count,x = "Response",
          y = "SORBS1", color = "Response", add = "jitter",ylab = "SORBS1 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Response),
            method = "t_test", label = "p.adj") 


ggboxplot(suni_select_count,x = "Response",
          y = "CD44", color = "Response", add = "jitter",ylab = "CD44 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Response),
            method = "t_test", label = "p.signif") 


ggboxplot(suni_select_count,x = "Response",
          y = "PCNA", color = "Response", add = "jitter",ylab = "PCNA expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Response),
            method = "t_test", label = "p.signif") 

ggboxplot(suni_select_count,x = "Response",
          y = "KDM5B", color = "Response", add = "jitter",ylab = "KDM5B expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Response),
            method = "t_test", label = "p.signif") 



ggboxplot(suni_select_count,x = "Response",
          y = "CDK2", color = "Response", add = "jitter",ylab = "CDK2 expression counts")+ scale_fill_manual()+ geom_pwc(
            aes(Response),
            method = "t_test", label = "p.signif") 



colData <- data.frame(row.names = rownames(targets),Response = targets$Response)

select_mat <- norm_expr[select_gene,]

library(pheatmap)


colData$Response = relevel(colData$Response, ref = "Response")


pheatmap(select_mat,
         annotation_col = colData,
         cluster_rows = F,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         scale = "row",
         color=colorRampPalette(c("navy", "white", "red"))(50))



library(tidyr)
dat_re <- suni_select_count %>%
  rownames_to_column("Sample") %>% 
  gather(key = Genes,value = Proportion,-Sample,-Response)


dat_re$Response <- factor(dat_re$Response ,levels = c("non_Response","Response"))


library(ggplot2)

library(ggpubr)

library(rstatix)
library(envalysis)




all_violin <- ggplot(dat_re , aes(x=Genes, y=Proportion,fill = Response)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(aes(group = Response), method = "t_test", label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)  + ylab(" Normalized expression level") + scale_fill_manual(values = c('Response'='#619CFF','non_Response' = '#F8766D'))+ theme_publish(base_size = 11, base_linewidth = 0.7)

all_violin

#########################################################################deg
Response <- factor(targets$Response)

Response

design <- model.matrix(~Response + 0)

table(Response)

colnames(design) <- levels(Response)


fit <- lmFit(eset_final, design)

contrast.matrix<-makeContrasts(paste0(c("non_Response","Response"),collapse = "-"),levels = design)


cont.matrix <- makeContrasts(nonvRes=non_Response-Response,
                             levels=design)


fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2,robust=TRUE, trend=TRUE) 

tempOutput = topTable(fit2, coef=1, n=Inf)

df <- topTable(fit2, coef=1, n=2000)

abs(tempOutput$logFC)

df <- tempOutput[tempOutput$P.Value < 0.05 & abs(tempOutput$logFC) > 1.2,]

tfit <- treat(fit2, lfc = 0.5)

tfit$coefficients

pval <- 0.05

fitSum <- summary(decideTests(tfit, p.value = pval))
fitSum
