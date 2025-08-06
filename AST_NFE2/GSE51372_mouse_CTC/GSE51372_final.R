gse51372 = read.table("GSE51372_readCounts.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, quote = "")

################################################3

CTC_plt_sample = c("MP4-1","MP4-3","MP4-4","MP4-6","MP4-7","MP4-8","MP4-13","MP4-14","MP4-17","MP4-20","MP4-22","MP4-28","MP4-29","MP4-31","MP4-32","MP6-2","MP6-3","MP6-4","MP6-5","MP6-7","MP6-9","MP6-10","MP6-17","MP6-20")

CTC_sample = c("MP2-1","MP2-2","MP2-4","MP2-11","MP2-17","MP2-18","MP2-20","MP2-21","MP2-24","MP2-26","MP2-30","MP2-32","MP2-36","MP4-24","MP6-6","MP6-11","MP6-15","MP6-16","MP6-18","MP7-1","MP7-3","MP7-4","MP7-7","MP7-8","MP7-9","MP7-12","MP7-13","MP7-16","MP7-18","MP7-20","MP7-25","MP7-29","MP7-30","MP7-31","MP7-33","MP7-34","MP7-37","MP7-40","MP7-41")


Tumor = gse51372[,c(4,128:147)]

#Tumor = gse51372[,c(4,148:181)]

Tumor_count <- Tumor
#20


#CTC_c_count = gse51372[,CTC_c_sample]
#39

CTC_count <- gse51372[,CTC_sample]
#39

CTC_plt_count = gse51372[,CTC_plt_sample]
#24



library(dplyr)

library(tibble)

multi_counts <- bind_cols(Tumor_count,CTC_count,CTC_plt_count)

multi_counts2 <- multi_counts[!(multi_counts$symbol == "" ), ]



multi_counts2 <- aggregate(multi_counts2[, -c(1)],
                           by = list(Gene = multi_counts2$symbol),
                           FUN = mean,
                           na.rm = TRUE)

rownames(multi_counts2) <- NULL

multi_counts2 <- column_to_rownames(multi_counts2,"Gene")

######################################################################

df<- data.frame(EntrezID = mm9_genes$gene_id, Symbol = mm9_genes$symbol, Gene_length = width(mm9_genes))

Gle  <- df[,c(2,3)]

all_Matrix_sub<-multi_counts2[rownames(multi_counts2 ) %in% df$Symbol, ]


le = Gle[match(rownames(all_Matrix_sub ),Gle$Symbol),"Gene_length"]

le

countToTpm <- function(counts, effLen)
{
  rate <- log(counts + 1) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

all_TPM<- apply(all_Matrix_sub,2, countToTpm, effLen =le)


tpm_log = log2(all_TPM+1)


ast = c("Nfe2","Itga2b","Gp1ba","Hba-a2","Hbb-b1")


AST_exp_tpm  = tpm_log[ast,]

AST_exp_tpm <- as.data.frame(t(AST_exp_tpm))

AST_exp_tpm$Group <- c(rep("PT",20), rep("CTC",39), rep("CTC_plt",24))


library(tidyr)

library(tibble)

dat_tpm <- AST_exp_tpm%>%
  rownames_to_column("Sample") %>% 
  gather(key = Genes,value = Proportion,-Sample,-Group)


dat_tpm$Group <- factor(dat_tpm$Group,levels = c("PT","CTC_plt","CTC"))


library(ggplot2)

library(ggpubr)

library(rstatix)
library(envalysis)

all_violin <- ggplot(dat_tpm , aes(x=Genes, y=Proportion,fill = Group)) + 
  geom_violin(trim = FALSE,position = position_dodge(width = 1),scale = 'width')+
  theme_bw() + geom_pwc(aes(group = Group), method = "t_test", label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 1),outlier.size = 0.7,width= 0.2,show.legend = FALSE)  + ylab(" log2(TPM+1) expression level") + scale_fill_manual(values = c('PT'='#619CFF','CTC_plt' = '#F8766D','CTC'=  "#00BA38"))+ theme_publish(base_size = 11, base_linewidth = 0.7)


ggsave(filename = "final_result/CTC_PT_PLT_all_violine.png", plot = all_violin , width = 12, height = 10, dpi = 300, units = "cm")
