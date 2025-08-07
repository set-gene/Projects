cells = read.table("rna_celline.tsv", sep = "\t",header = TRUE, stringsAsFactors = FALSE)


library(tidyverse)
cells %>% group_by(Cell.line) %>% summarize (num_genes = n())

library(DESeq2)


tpm <- reshape2::dcast (cells, Gene ~Cell.line, value.var = "TPM")


#tpm <- column_to_rownames(tpm, var = "Gene")

library(clusterProfiler)

library(org.Hs.eg.db)

genes <- bitr(tpm$Gene, fromType="ENSEMBL", 
           toType="SYMBOL", 
           OrgDb="org.Hs.eg.db")

colnames(genes)[1]= "Gene"

tpm2 = merge(genes,tpm)

exprs = aggregate(tpm2[,-c(1,2)],
                  by = list(Gene = tpm2$SYMBOL),
                  FUN = mean,
                  na.rm = TRUE)


exprs <- column_to_rownames(exprs, var = "Gene")

adherent = c("A-431", "A549", "AF22", "AN3-CA", "ASC diff", "ASC TERT1"
             , "BEWO"	,"BJ", "BJ hTERT+", "BJ hTERT+ SV40 Large T+", "BJ hTERT+ SV40 Large T+ RasG12V","CACO-2", "CAPAN-2",
             "EFO-21", "fHDF/TERT166", "HaCaT","HAP1", "HBEC3-KT", "HBF TERT88", "HEK 293", "HeLa", "Hep G2","HHSteC",
             "HSkMC", "hTCEpi", "hTEC/SVTERT24-B", "hTERT-HME1", "HUVEC TERT2", "LHCN-M2",
             "MCF7",	"NTERA-2", "PC-3", "RH-30", "RPTEC TERT1", "RT4",
             "SH-SY5Y", "SiHa", "SK-BR-3", "SK-MEL-30", "T-47d","TIME", "U-138 MG", 
             "U-2 OS", "U-2197", "U-251 MG", "U-87 MG")


suspension = c("WM-115","Daudi","HDLM-2","HEL"	,"HL-60"	,"HMC-1","K-562"	,"Karpas-707"
               ,"MOLT-4" ,"NB-4","REH","RPMI-8226",	"SCLC-21H","THP-1","U-266/70",
               "U-266/84"	,"U-698"	,"U-937")


adherent_tpm = exprs[,adherent]

suspension_tpm = exprs[,suspension]

adherent = data.frame(Sample = adherent)

adherent$Group = rep("Adherent",46)

suspension = data.frame(Sample = suspension)

suspension$Group = rep("Suspension",18)

exprs2 = cbind(adherent_tpm,suspension_tpm)

colData = rbind(adherent,suspension)

colData$Group = factor((colData$Group), levels = c("Adherent", "Suspension"))

colData$Group = relevel(colData$Group, ref = "Adherent")

library(DESeq2)

exprs2 = as.matrix(exprs2)

colData = column_to_rownames(colData, "Sample")

dds = DESeq2::DESeqDataSetFromMatrix(countData = round(exprs2 * 10),
                             colData = colData,
                             design = ~Group)

dds <- dds[ rowSums(counts(dds)) > 1, ]  

dds <- DESeq(dds)

res = DESeq2::results(dds  ,c("Group", "Suspension","Adherent"))


res2 = res[order(res$padj),]

res2

cut_pvalue <- 0.05

cut_lfc <- 1

significant_results <- res2[which(res2$padj < cut_pvalue),]

significant_results

significant_results_EXP <- res2[which(res2$padj < cut_pvalue & (res2$log2FoldChange<(-cut_lfc) | res2$log2FoldChange>cut_lfc)),]

significant_results_EXP = data.frame(significant_results_EXP)

significant_results_EXP_df = rownames_to_column(significant_results_EXP,"Gene")

rld <- DESeq2::rlog(dds, blind = FALSE)

rlogMat <- assay(rld)

ast_pla = c("IKZF1","NFE2","CDH1","VIM","SELP","SLC6A4","MYL9","SLC2A3","LRP8","FYB1","SCAMP2","PTPN1","FGB","FGG")

pla = c("SELP","SLC6A4","MYL9","SLC2A3","LRP8","FYB1","SCAMP2","PTPN1","FGA","FGB","FGG")

ast = c("IKZF1","NFE2","IRF8","BTG2","CD37")

AST_exp  = rlogMat[ast_pla,]

AST_1exp = rlogMat[ast,]

pla_exp = rlogMat[pla,]

AST_exp_ad = AST_exp[,c(1:46)]

AST_exp_sus = AST_exp[,c(47:64)]

NN_select_sus  = t(AST_exp_sus)

library(Hmisc)


NN_select_ad  = t(AST_exp_ad)

library(Hmisc)

ccor_select_ad <- rcorr(NN_select_ad, type ="pearson")


cor_mx_ad = as.matrix(ccor_select_ad$r)

cor_mx_ad_p = as.matrix(ccor_select_ad$P)


mycol_P <- ifelse(c(cor_mx_ad_p < 0.05), "black", "white")

library(corrplot)

cor_AST_ad = corrplot(cor_mx_ad ,        # Correlation matrix
                      method = "color", # Correlation plot method
                      type = "full",    # Correlation plot style (also "upper" and "lower")
                      diag = TRUE,      # If TRUE (default), adds the diagonal
                      # Labels color
                      bg = "white",     # Background color
                      title = "",       # Main title
                      col = NULL,
                      order = "original",
                      addCoef.col = mycol_P,
                      addCoefasPercent = TRUE
)

################################################################################

library(Hmisc)

ccor_select_sus <- rcorr(NN_select_sus, type ="pearson")


cor_mx_sus = as.matrix(ccor_select_sus$r)

cor_mx_sus_p = as.matrix(ccor_select_sus$P)


mycol <- ifelse(c(cor_mx_sus_p < 0.05), "black", "white")

library(corrplot)

cor_AST_sus = corrplot(cor_mx_sus ,        # Correlation matrix
                      method = "color", # Correlation plot method
                      type = "full",    # Correlation plot style (also "upper" and "lower")
                      diag = TRUE,      # If TRUE (default), adds the diagonal
                      # Labels color
                      bg = "white",     # Background color
                      title = "",       # Main title
                      col = NULL,
                      order = "original",
                      addCoef.col = mycol,
                      addCoefasPercent = TRUE
)







##############################################################################
library(pheatmap)

annotation_col=data.frame(group=colData$Group)

row.names(annotation_col) <- colnames(AST_exp)

annotation_col

AST_1expT = t(AST_1expT)


pheatmap(rlogMat[ast,] ,
         annotation_col = annotation_col,
         cluster_rows = T,
         show_rownames = T,
         show_colnames = T,
         border_color = NA,
         fontsize = 7,
         scale = "row",
         fontsize_row = 8,
         annotation_colors = list(group=c(Adherent = "blue",Suspension="red")))

selects = data.frame(t(AST_1exp))

selects$Group = c(rep("Adhesion",46), rep("Suspension",18))


library(ggplot2)

library(ggpubr)

norm_count = DESeq2::counts(dds)

norm_ast = norm_count[ast,]

norm_ast = data.frame(t(norm_ast))

norm_ast$Group = c(rep("Adhesion",46), rep("Suspension",18))

ggboxplot(selects,x = "Group",
          y = "CD37", color = "Group", add = "jitter",ylab = "CD37 normalization counts ") + stat_compare_means(method = "t.test")

#############################################3


ast = c("IKZF1","NFE2","IRF8","BTG2","TEAD2","HBA1","HBA2","GATA1")

select_genes = c("PECAM1","CD34","PTPRC","STAT3","NFE2","GATA1","MMP9","AKT1")

#######VOLCANO
topT <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(topT)[1] <- "Gene"

library(tibble)
# Adjusted P values
with(topT, plot(log2FoldChange, -log10(padj), main ="The volcano plot",pch=20,col='grey', cex=1.0, xlab=bquote(~Fold~change), ylab=bquote(~-log[10]~(adj~p~value))))

with(topT,legend("topright",c("Adherent","Suspension"),col = c("grey","purple"),pch = 20, title = "Group"))

with(subset(topT, padj<cut_pvalue & log2FoldChange>cut_lfc),points(log2FoldChange, -log10(padj), pch=20, col='purple', cex=1.5))


## Add lines for FC and P-value cut-off
abline(v=0, col='black', lty=4, lwd=2.0)


abline(h=-log10(max(topT$padj[topT$padj<cut_pvalue], na.rm=TRUE)), col='black', lty=4, lwd=2.0)

library(tibble)

topTT = column_to_rownames(topT,var = "Gene")

intentded_gene = c("IKZF1","NFE2","CDH1","VIM","SELPLG","MYL9","FYB1")

with(topTT[ast, ], {
  points(log2FoldChange, -log10(padj), col="red", cex=2, lwd=2)
  text(log2FoldChange, -log10(padj), ast, pos=2, col="red")
})
