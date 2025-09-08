RCC_786O_25mM_1 <- read.table("786-O_25mM-1.txt", sep = "\t",header = T,quote = "")

RCC_786O_25mM_2 <- read.table("786-O_25mM-2.txt", sep = "\t",header = T,quote = "")

RCC_786O_25mM_3 <-  read.table("786-O_25mM-3.txt", sep = "\t",header = T,quote = "")

RCC_786O_5mM_1 <- read.table("786-O_5mM-1.txt", sep = "\t",header = T,quote = "")

RCC_786O_5mM_2 <- read.table("786-O_5mM-2.txt", sep = "\t",header = T,quote = "")

RCC_786O_5mM_3 <-  read.table("786-O_5mM-3.txt", sep = "\t",header = T,quote = "")

colnames(RCC_786O_25mM_1)[7]<-("786O_25mM_1")
colnames(RCC_786O_25mM_2)[7]<-("786O_25mM_2")
colnames(RCC_786O_25mM_3)[7]<-("786O_25mM_3")
colnames(RCC_786O_5mM_1)[7]<-("786O_5mM_1")
colnames(RCC_786O_5mM_2)[7]<-("786O_5mM_2")
colnames(RCC_786O_5mM_3)[7]<-("786O_5mM_3")

RCC_786O_25mM_1 <- RCC_786O_25mM_1[,c(1,7)]
RCC_786O_25mM_2 <- RCC_786O_25mM_2[,c(1,7)]
RCC_786O_25mM_3 <- RCC_786O_25mM_3[,c(1,7)]
RCC_786O_5mM_1 <- RCC_786O_5mM_1[,c(1,7)]
RCC_786O_5mM_2 <- RCC_786O_5mM_2[,c(1,7)]
RCC_786O_5mM_3 <- RCC_786O_5mM_3[,c(1,7)]

library(dplyr)

n_raw_count <- merge(RCC_786O_25mM_1,RCC_786O_25mM_2,by="Geneid") %>%
  merge(.,RCC_786O_25mM_3,by = "Geneid") %>%
  merge(.,RCC_786O_5mM_1,by = "Geneid") %>%
  merge(.,RCC_786O_5mM_2,by = "Geneid") %>%
  merge(.,RCC_786O_5mM_3,by = "Geneid")
  
  
library(tibble)

n_raw_count <- column_to_rownames(n_raw_count,"Geneid")


Condition <- factor(c(rep("7860_25mM",3),rep("7860_5mM",3)), levels = c("7860_25mM","7860_5mM"))


colData <- data.frame(row.names = colnames(n_raw_count),Condition)


colData$Condition = relevel(colData$Condition, ref = "7860_25mM")

library(DESeq2)

dds <-DESeqDataSetFromMatrix(n_raw_count,colData,design = ~ Condition )


dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq2::DESeq(dds)

dat  <- counts(dds, normalized = TRUE)

idx  <- rowMeans(dat) > 1


resultsNames(dds)

res = DESeq2::results(dds ,c("Condition", "7860_5mM","7860_25mM"))

res2 = res[order(res$padj),]

value <- 0.05

cut_lfc <- 2

significant_results <- res2[which(res2$padj < cut_pvalue),]

significant_results


significant_results_EXP <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_EXP = data.frame(significant_results_EXP)

significant_results_High <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange>cut_lfc)),]

significant_results_High <- rownames_to_column(as.data.frame(significant_results_High), "Gene")

significant_results_Low<- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) )),]

significant_results_Low <- rownames_to_column(as.data.frame(significant_results_Low), "Gene")



library(openxlsx)

write.xlsx(significant_results_High,"786O_5mM_significant_LFC2_High.xlsx")

write.xlsx(significant_results_Low,"786O_5mM_significant_LFC2_Low.xlsx")


########################################################
cut_lfc <- 1

cut_pvalue <- 0.05



topT <- as.data.frame(res)

topT <- rownames_to_column(topT,"Gene")

# Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="786-O_5mM vs 786-O_25mM Volcano plot", col='grey', lwd = 5.0,cex = 4,cex.lab = 1,cex.axis = 1,cex.main=3, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value), xlim = c(-10,10)))

with(topT,legend("topright",c("786-O_25mM ","stable","786-O_5mM"),col = c("green","grey","purple"),pch = 20, title = "Condition",box.lwd = 2,cex = 1.0),font.lab = 2)


with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc),points(log2FoldChange, -log10(padj), pch=20, col='purple', cex=5))

with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=5))

## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=4, lwd=5)

abline(v=-cut_lfc, col='black', lty=4, lwd=5)

abline(v=cut_lfc, col='black', lty=4, lwd=5)

abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=5.0)



library(calibrate)

with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc), textxy(log2FoldChange, -log10(padj), labs=Gene, font = 2,cex = 1))

with(subset(topT, padj<cut_pvalue & log2FoldChange<(-cut_lfc)), textxy(log2FoldChange, -log10(padj), labs=Gene,  font = 2,cex = 1))


############################################################################################


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


go_u_bp <- enrichGO(eg_u$SYMBOL, OrgDb = "org.Hs.eg.db", ont="BP",keyType = "SYMBOL")

go_d_bp <- enrichGO(eg_d$SYMBOL, OrgDb = "org.Hs.eg.db", ont="BP",keyType = "SYMBOL")


dotplot(go_d_bp,showCategory = 30,font.size = 15)

go_u_bp_df <- as.data.frame(go_u_bp@result)

go_u_bp_df <- subset(go_u_bp_df,go_u_bp_df$p.adjust < 0.05)

library(openxlsx)

write.xlsx(go_u_bp_df,"significant_UP_LFC2_GO_BP.xlsx",rowNames = FALSE)

dotplot(go_u_bp,showCategory = 30,font.size = 15)


go_d_bp_df <- as.data.frame(go_d_bp@result)

go_d_bp_df <- subset(go_d_bp_df,go_d_bp_df$p.adjust < 0.05)

write.xlsx(go_d_bp_df,"significant_down_LFC2_GO_BP.xlsx",rowNames = FALSE)

###################################################################3
library(tibble)

library(dplyr)

library(fgsea)

library(DESeq2)

res22 = rownames_to_column(data.frame(res2), "SYMBOL")


res3 <- res22 %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))


ranks <- deframe(res3)



canonical = gmtPathways("c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt")

canonical %>% 
  head() %>% 
  lapply(head)


canonical_Res <- fgsea(pathways=canonical, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)

GO = gmtPathways("c5.go.bp.v2023.2.Hs.symbols (1).gmt")

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

go_ResT_all_df <- rbind(go_ResT_up,go_ResT_dn)

go_ResT_all_df <- subset(go_ResT_all_df,go)


library(ggplot2)

ggplot(go_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO BP pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))


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

write.xlsx(canonical_ResT_up,"canonical_up_significant_GSEA.xlsx")

canonical_ResT_dn = subset(x = canonical_ResT, canonical_ResT$NES < -1)

canonical_ResT_dn <- canonical_ResT_dn[order(canonical_ResT_dn$padj),]


canonical_ResT_dn <- canonical_ResT_dn[order(canonical_ResT_dn$padj),]

canonical_ResT_dn_sig <- subset(canonical_ResT_dn,canonical_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(canonical_ResT_dn,"canonical_dn_significant_GSEA.xlsx")


topcanonicalPathwaysUp <- canonical_Res[ES > 0][head(order(pval), n=20), pathway]

topcanonicalPathwaysDown <- canonical_Res[ES < 0][head(order(pval), n=20), pathway]

topPathways <- c(topcanonicalPathwaysUp , rev(topcanonicalPathwaysDown))


canonical_ResT_up_FOR = subset(canonical_ResT_up,canonical_ResT_up$pathway %in% topcanonicalPathwaysUp)


canonical_ResT_dn_FOR = subset(canonical_ResT_dn,canonical_ResT_dn$pathway %in% topcanonicalPathwaysDown)

canonical_ResT_all <- rbind(canonical_ResT_up_FOR,canonical_ResT_dn_FOR)
##############################################################################

library(ggplot2)

ggplot(canonical_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))

##############################################


GO = gmtPathways("c5.go.bp.v2023.2.Hs.symbols (1).gmt")

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

go_ResT_all_df <- rbind(go_ResT_up,go_ResT_dn)

go_ResT_all_df <- subset(go_ResT_all_df,go)


library(ggplot2)

ggplot(go_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO BP pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))
#############################################################

library(fgsea)

library(dplyr)

hm = gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

hm_Res <- fgsea(pathways=hm, stats=ranks, nperm=1000, minSize = 15, maxSize = 500)


hm_ResT <- hm_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

hm_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

hm_ResT_up = subset(x = hm_ResT, hm_ResT$NES > 1)

hm_ResT_up <- hm_ResT_up[order(hm_ResT_up$padj),]

hm_ResT_up_sig <- subset(hm_ResT_up,hm_ResT_up$padj < 0.05)

library(openxlsx)


write.xlsx(hm_ResT_up_sig,"HM_up_significant_GSEA.xlsx")

hm_ResT_dn = subset(x = hm_ResT, hm_ResT$NES < -1)


hm_ResT_dn <- hm_ResT_dn[order(hm_ResT_dn$padj),]


hm_ResT_dn_sig <- subset(hm_ResT_dn,hm_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(hm_ResT_dn_sig,"hm_dn_significant_GSEA.xlsx")


tophmPathwaysUp <- hm_Res[ES > 0][head(order(pval), n=20), pathway]

tophmPathwaysDown <- hm_Res[ES < 0][head(order(pval), n=20), pathway]


hm_ResT_up_FOR = subset(hm_ResT_up,hm_ResT_up$pathway %in% tophmPathwaysUp)


hm_ResT_dn_FOR = subset(hm_ResT_dn,hm_ResT_dn$pathway %in% tophmPathwaysDown)

hm_ResT_all <- rbind(hm_ResT_up_FOR,hm_ResT_dn_FOR)

hm_ResT_all_df <- rbind(hm_ResT_up,hm_ResT_dn)

#go_ResT_all_df <- subset(hm_ResT_all_df,go)


library(ggplot2)

ggplot(hm_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))
####################################################################

library(fgsea)

library(dplyr)

wp = gmtPathways("c2.cp.wikipathways.v2024.1.Hs.symbols.gmt")

wp_Res <- fgsea(pathways=wp, stats=ranks, nperm=1000, minSize = 5, maxSize = 500)


wp_ResT <- wp_Res %>%
  as_tibble() %>%
  arrange(desc(NES))

wp_ResT  %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

wp_ResT_up = subset(x = wp_ResT, wp_ResT$NES > 1)

wp_ResT_up <- wp_ResT_up[order(wp_ResT_up$padj),]

wp_ResT_up_sig <- subset(wp_ResT_up,wp_ResT_up$padj < 0.05)

library(openxlsx)


write.xlsx(hm_ResT_up_sig,"HM_up_significant_GSEA.xlsx")

wp_ResT_dn = subset(x = wp_ResT, wp_ResT$NES < -1)


wp_ResT_dn <- wp_ResT_dn[order(wp_ResT_dn$padj),]


wp_ResT_dn_sig <- subset(wp_ResT_dn,wp_ResT_dn$padj < 0.05)

library(openxlsx)

write.xlsx(hm_ResT_dn_sig,"hm_dn_significant_GSEA.xlsx")


topwpPathwaysUp <- wp_Res[ES > 0][head(order(pval), n=20), pathway]

topwpPathwaysDown <- wp_Res[ES < 0][head(order(pval), n=20), pathway]


wp_ResT_up_FOR = subset(wp_ResT_up,wp_ResT_up$pathway %in% topwpPathwaysUp)


wp_ResT_dn_FOR = subset(wp_ResT_dn,wp_ResT_dn$pathway %in% topwpPathwaysDown)

wp_ResT_all <- rbind(wp_ResT_up_FOR,wp_ResT_dn_FOR)

wp_ResT_all_df <- rbind(wp_ResT_up,wp_ResT_dn)

#go_ResT_all_df <- subset(hm_ResT_all_df,go)


library(ggplot2)

ggplot(wp_ResT_all , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways from GSEA") + 
  theme_minimal()+ scale_fill_manual(values = c("TRUE" = "blue","FALSE"= "red"))
################################################################

library(DESeq2)

rld <- DESeq2::rlog(dds, blind = FALSE)

rlogMat <- assay(rld)

low_df <- rlogMat[,4:6]


hypoxia_up <- gmtPathways("HALLMARK_HYPOXIA.v2023.2.Hs.gmt")


hypoxia_up<- hypoxia_up$HALLMARK_HYPOXIA


Gene_set= list(hypoxia_up)

names(Gene_set) = c("Hypoxia")


library(GSVA)

library(dplyr)

library(Hmisc)

library(rstatix)
library(dplyr)

library(psych)

library(future)

source("http://www.sthda.com/upload/rquery_cormat.r")


plan("multisession", workers = 8)

ssgsea = gsva(low_df, Gene_set,method = "ssgsea")

sig_high <- low_df[rownames(significant_results_up),]


sig_low <- low_df[rownames(significant_results_down),]


for_cor_high <- rbind(ssgsea,sig_high)

for_cor_low <- rbind(ssgsea,sig_low)

for_cor_high_t <- t(for_cor_high)


for_cor_low_t <- t(for_cor_low)

ccor_low_df_t <- corr.test(for_cor_low_t)

ccor_low_df_t$r

#ccor_low_df_t <- cor(for_cor_t , method = "pearson")

library(Hmisc)

library(rstatix)
library(dplyr)

library(psych)

library(future)

source("http://www.sthda.com/upload/rquery_cormat.r")


plan("multisession", workers = 8)
#######################################################################


#ccor2_low_df_t <- rquery.cormat(for_cor_t,type = "full" )

#ccor2_low_df_t$r

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

filter_cor_low = flattenCorrMatrix(ccor_low_df_t$r, ccor_low_df_t$p)
#####################################################################
########################################################################

library(Hmisc)

library(rstatix)
library(dplyr)

library(psych)

library(future)

library(tibble)

source("http://www.sthda.com/upload/rquery_cormat.r")


sig_low <- low_df[rownames(significant_results_EXP),]


for_cor_low <- rbind(ssgsea,sig_low)

for_cor_t <- t(for_cor_low)

ccor_low_df_t <- rquery.cormat(for_cor_t, type="flatten", graph=FALSE)

ccor_df <- ccor_low_df_t$r

ccor_df_sig <- subset(ccor_df, ccor_df$p < 0.05)

hp_ccor_df <- subset(ccor_df_sig,ccor_df_sig$column == "Hypoxia")

significant_results_EXP_df <- rownames_to_column(significant_results_EXP,"Gene")

colnames(hp_ccor_df)[1] <- "Gene"

hp_ccor_df_exp <- merge(hp_ccor_df,significant_results_EXP_df,"Gene")

library(openxlsx)

write.xlsx(hp_ccor_df_exp ,"hypoxia_cor_df.xlsx")


ccor2_low_df_t <- cor(for_cor_t , method = "pearson")

ccor2_low_df_t_df <- rownames_to_column(data.frame(ccor2_low_df_t),"row")

hypox_cor <- ccor2_low_df_t_df[,c(1:2)]

hypox_cor<- hypox_cor[order(hypox_cor$Hypoxia,decreasing = TRUE),]


write.xlsx(hypox_cor ,"hypoxia_only_cor_df.xlsx")

##################################################

library(enrichplot)

library(ggplot2)

library(tibble)

library(ggnewscale)

library(clusterProfiler)


ranks_df <- as.data.frame(ranks)

ranks_df <- rownames_to_column(ranks_df,"SYMBOL")



ranks_et= bitr(ranks_df$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

ranks_df_entrez <- merge(ranks_et,ranks_df,by= "SYMBOL")

ranks_df_entrez2 <- ranks_df_entrez[-which(duplicated(ranks_df_entrez$SYMBOL)),]

geneList = ranks_df_entrez2[,3]

names(geneList) = as.character(ranks_df_entrez2[,2])

geneList = sort(geneList, decreasing = TRUE)

geneList

geneList2 = ranks_df_entrez2[,3]

names(geneList2) = as.character(ranks_df_entrez2[,1])


geneList2 = sort(geneList2, decreasing = TRUE)

geneList2


gseWP<- gseWP(geneList     = geneList,
                 organism     = "Homo sapiens",
                 nPerm        = 1000,
                 minGSSize    = 5,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)




gseWP_df <- gseWP@result

###################################################

library(QuaternaryProd)

e2f3 <- system.file("extdata", "e2f3_sig.txt",
                    package = "QuaternaryProd")

e2f3 <- read.table(e2f3, sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE)

names(e2f3) <- c("entrez", "pvalue", "fc")
e2f3 <- e2f3[!duplicated(e2f3$entrez),]


quaternary_results <- RunCRE_HSAStringDB(e2f3, method = "Quaternary",
                                         fc.thresh = log2(1.3), pval.thresh = 0.05,
                                         only.significant.pvalues = TRUE,
                                         significance.level = 0.05,
                                         epsilon = 1e-16, progressBar = FALSE,
                                         relations = NULL, entities = NULL)
##################################################################################
library(pheatmap)

ferroptosis_related <- c( "UGT1A6", "ATP6V1C1", "MAFG", "NUDCD1", "PPP1R1A", "TSKU", "CTSB", "AIFM2", "CTSA",  "CTNND2")

#http://www.zhounan.org/ferrdb/current/operations/browsegene.html?genetype=marker
#FerrDBv2 
ferroptosis_marker <- c( "FTH1", "GPX4", "CHAC1", "HSPB1", "NFE2L2", "PTGS2", "SLC40A1", "TF", "TFRC","ACSL4")

annotation_col=data.frame(group=colData$condition)

row.names(annotation_col) <- colnames(rld_counts_sig_ALL)

colnames(annotation_col) <- "Group"

pheatmap(rlogMat[ferroptosis_marker,] ,
         annotation_col = colData,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         scale = "row",
         color=colorRampPalette(c("navy", "white", "red"))(50))

