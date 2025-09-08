exp_kirc <- read.table("TCGA.KIRC.sampleMap_HiSeqV2/HiSeqV2",sep = "\t",header = TRUE,check.names = FALSE)


clin <- read.table("TCGA.KIRC.sampleMap_KIRC_clinicalMatrix",sep = "\t",header = TRUE,check.names = FALSE)

survival<- read.table("survival_KIRC_survival.txt",sep = "\t",header = TRUE,check.names = FALSE)



library(tibble)

exp_kirc <- column_to_rownames(exp_kirc,"sample")

#################################################################selection samples
ls <- unlist(substr(colnames(exp_kirc),14,15))

table(ls) 

group_list = ifelse(as.numeric(substr(colnames(exp_kirc),14,15)) < 2, 'Tumor', 'normal')

table(group_list) 


count_matrix = na.omit(exp_kirc)

count_matrix_tu <-  count_matrix[,group_list=="Tumor"]

count_matrix_n <-  count_matrix[,group_list=="normal"]

count_matrix3 <- cbind(count_matrix_n,count_matrix_tu)


save(list = c("count_matrix_tu","count_matrix3","count_matrix"),file = "KIRC_tumor.RData")

#############################################################################################


intented_gene <- c("LCN2","ISCA2","OGT","ACSL3","GABPA","ENO2","HMOX1","OTUD4","KLF2",
                   "KLF11","MITD1","TRIB3","MDH2")

select_count <- data.frame(t(count_matrix_tu[intented_gene,]))


clin_tu <- subset(clin,clin$sampleID %in% rownames(select_count))

surv_tu <- subset(survival,survival$sample %in% rownames(select_count))


rownames(surv_tu) <- surv_tu$sample

surv_filtered <- surv_tu[match(rownames(select_count), surv_tu$sample), ]

library(survminer)
library(survival)

library(ggpubr)

library(rstatix)


for (gene in intented_gene) {
  

  Group <- ifelse(select_count[[gene]] > median(select_count[[gene]], na.rm = TRUE), "High", "Low")
  

  df <- data.frame(sample = surv_filtered$sample,
    time = surv_filtered$OS.time,
    status = surv_filtered$OS,
    Group = Group
  )
  
  
  # Surv 객체 생성
  surv_obj <- Surv(df$time, df$status)
  
  # 모델 적합
  fit <- survfit(surv_obj ~ Group, data = df)
  
  # 생존곡선 출력
  
  surv_plot <- ggsurvplot(
    fit,
    data = df,
    title = paste("Overall Survival Plot:", gene),
    legend = "top",
    legend.title = "Group",
    palette = c("red","blue"),
    xlab="Days",
    ylab="Survival probability",
    font.x = 11,
    font.y = 11,
    pval = TRUE,
    pval.method = TRUE,
    pval.size = 3,
    pval.coord = c(50,.24), 
    risk.table = FALSE,
    ggtheme = theme_classic2(base_size=11)
  )
  

  ggsave(
    filename = paste0("survival_plot_", gene, ".png"),
    plot = surv_plot$plot,
    width = 12, height = 10, dpi = 300, units = "cm"
  )
  
  # 리스크 테이블 포함된 전체를 저장하려면 다음을 사용
  # ggsave(paste0("survival_plot_", gene, "_full.png"), surv_plot, width = 7, height = 6, dpi = 300)
}

for (gene in intented_gene) {
  
  # 그룹 지정: 유전자 발현값 기준으로 High / Low 구분
  Group <- ifelse(select_count[[gene]] > median(select_count[[gene]], na.rm = TRUE), "High", "Low")
  
  # 생존 데이터 병합
  df2 <- data.frame(sample = surv_filtered$sample,
                   time = surv_filtered$PFI.time,
                   status = surv_filtered$PFI,
                   Group = Group
  )
  
  

  surv_obj2 <- Surv(df2$time, df2$status)
  

  fit2 <- survfit(surv_obj2 ~ Group, data = df2)
  

  surv_plot2 <- ggsurvplot(
    fit2,
    data = df2,
    title = paste("Progression-Free Interval Survival Plot:", gene),
    legend = "top",
    legend.title = "Group",
    palette = c("red","blue"),
    xlab="Days",
    ylab="Survival probability",
    font.x = 11,
    font.y = 11,
    pval = TRUE,
    pval.method = TRUE,
    pval.size = 3,
    pval.coord = c(50,.24), 
    risk.table = FALSE,
    ggtheme = theme_classic2(base_size=11)
  )
  
  ggsave(
    filename = paste0("PFI_survival_plot_", gene, ".png"),
    plot = surv_plot2$plot,
    width = 12, height = 10, dpi = 300, units = "cm"
  )
  
  # ggsave(paste0("survival_plot_", gene, "_full.png"), surv_plot, width = 7, height = 6, dpi = 300)
}
