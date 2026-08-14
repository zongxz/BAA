rm(list = ls())
library(survival)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)


# RCS====
library(rms)

data <- fread('data.csv')
data <- data.frame(data)

HR_all_diseases <- data.frame()

for (disease in diseases) {
  
  dis_time <- paste0(disease, '_time')
  print(dis_time)
  data.rcs <- data[, c(dis_time, disease, 'age_i2', 'age_acceleration', 
                       "sex", 'education_class','Mixed' , 'Asian', 'Black', 'Chinese', 'others')]
  
  colnames(data.rcs) <- c('time', 'status', 'age_i2', 'age_acceleration', "sex", 'education_class',
                          'Mixed' , 'Asian', 'Black', 'Chinese', 'others')
  
  dd <- datadist(data.rcs) 
  options(datadist = 'dd') 
  aic_values <- sapply(3:5, function(k) {
    model <- cph(Surv(time, status) ~ rcs(age_acceleration, k) + sex + age_i2 + education_class +
                   + Mixed + Asian + Black + Chinese + others, 
                 data = data.rcs)
    AIC(model)
  })
  
  names(aic_values) <- c("3", "4", "5")
  best_k <- as.numeric(names(which.min(aic_values)))
  dd$limits$age_acceleration[2] <- median(data.rcs$age_acceleration)
  fit1 <- cph(Surv(time, status) ~ rcs(age_acceleration, best_k) + sex + age_i2 + education_class +
                + Mixed + Asian + Black + Chinese + others, 
              data = data.rcs)
  
  HR1 <- Predict(fit1, age_acceleration, fun = exp, ref.zero = TRUE)
  HR1$disease <- disease
  HR_all_diseases <- rbind(HR_all_diseases, HR1)
  HR_all_diseases <- HR_all_diseases[,1:13]
}

# cox ====


for (disease in diseases) {
  dis_time <- paste0(disease, '_time')
  
  cox_formula <- as.formula(
    paste0('Surv(', dis_time, ', ', disease, ') ~ age_acceleration + age_i2 + sex + education_class + Mixed + Asian + Black + Chinese + others')
  )
  cox_model <- coxph(cox_formula, data = data)
  print(summary(cox_model))
  model_summary <- summary(cox_model)
  col <- c(colnames(model_summary$coefficients), colnames(model_summary$conf.int))
  results_matrix <- matrix(nrow = nrow(model_summary$coefficients), ncol = 9)
  
  results_matrix[,1:5] <- model_summary$coefficients
  results_matrix[,6:9] <- model_summary$conf.int
  colnames(results_matrix) <- col
  rownames(results_matrix) <- rownames(model_summary$coefficients)
  
  predictors <- model.matrix(cox_model)[, -1] 
  cor_matrix <- cor(predictors)
  
 
}

# logit ====
diseases <- c('cancer', 'CVD','digestive', 'respiratory','endocrine',
              'neurological',  'psychiatric','musculoskeleta', 'genitourinary', 
              'eye', 'ear', 'skin')
for (disease in diseases) {
  data.lo <- data[, c(disease,'age_i2',
                      'age_acceleration',"sex",'education_class',
                      'Mixed' , 'Asian', 'Black', 'Chinese', 'others')]
  colnames(data.lo) <- c('status','age_i2',
                         'age_acceleration',"sex",'education_class',
                         'Mixed' , 'Asian', 'Black', 'Chinese', 'others')
  logit_model <- glm(status ~ age_acceleration + age_i2 + 
                       sex + education_class +Mixed + 
                       Asian + Black + Chinese + others, 
                     data = data.lo, 
                     family = binomial(link = "logit"))
  summary(logit_model)
  result <- data.frame(summary(logit_model)$coefficients)
  

  coefficients <- coef(logit_model)
  std_errors <- sqrt(diag(vcov(logit_model)))
  
  result['OR'] <- exp(coefficients)
  result['lower.95'] <- exp(coefficients - 1.96 * std_errors)
  result['upper.95'] <- exp(coefficients + 1.96 * std_errors)

}


# ====
library(emmeans)

data.acc <- data
diseases <- c('cancer', 'CVD','digestive', 'neurological',  'psychiatric',
              'respiratory','endocrine', 'genitourinary', 'musculoskeleta',
              'eye', 'ear', 'skin')
for (disease in diseases) {
  threshold_age <- quantile(data.acc$age_i2[data.acc[disease] == 1], 0.10, na.rm = TRUE)
  
  data.acc <- data.acc %>%
    mutate(!!paste0(disease, "_group") := case_when(
      !!sym(disease) == 0 ~ "Control",
      !!sym(disease) == 1 & age_i2 <= threshold_age ~ "Early-onset",
      !!sym(disease) == 1 & age_i2 > threshold_age ~ "Other-onset"
    ))
}


library(emmeans)


all_results <- data.frame()

for (disease in diseases) {
  dis_group <- paste0(disease, '_group')
  
  model <- lm(age_acceleration ~ group + age_i2,data=data)
  summary(model)
  emm <- emmeans(model, ~group)
  summary_results <- summary(emm, infer = TRUE)
  summary_results <- summary_results %>%
    mutate(significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    disease = disease)
  
  pairwise_results <- data.frame(pairs(emm, adjust = "BH"))
  
  pairwise_results <- pairwise_results %>%
    mutate(
      significance = case_when(
        pairwise_results$p.value < 0.001 ~ "***",
        pairwise_results$p.value < 0.01 ~ "**",
        pairwise_results$p.value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      disease = disease
    )
  
  summary_results['const'] <- pairwise_results$contrast
  summary_results['const_estimate'] <-  pairwise_results$estimate
  summary_results['const_SE'] <-  pairwise_results$SE
  summary_results['const_p'] <-  pairwise_results$p.value
  summary_results['const_sig'] <- pairwise_results$significance
  all_results <- rbind(all_results, summary_results)
}


# correlation#####


library(ppcor)

# i0
for (cor in cors_i0) {
  data.cor <- data[,c("sex", "age_i2", "age_acceleration",cor,
                       "education_class", 'Mixed' , 'Asian', 'Black', 'Chinese', 
                       'others', 'p53_time')]
  data.cor <- na.omit(data.cor)
  
  pcor_result <- pcor.test(
    x = data.cor$age_acceleration, 
    y = data.cor[cor],
    z = data.cor[, c("sex", "age_i2", "education_class",
                      'Mixed' , 'Asian', 'Black', 'Chinese', 
                      'others', 'p53_time')],
    method = "pearson"  
  )
  pcor_result['R2'] <- pcor_result$estimate * pcor_result$estimate
  row.names(pcor_result) <- cor
  pcor_results <- rbind(pcor_results, pcor_result)
}
# i2
cors_i2 <- c("father_death_age","mother_death_age")
for (cor in cors_i2) {
  data.cor <- data[,c("sex", "age_i2", "age_acceleration",cor,
                       "education_class", 'Mixed' , 'Asian', 'Black', 'Chinese',
                       'others')]
  data.cor <- na.omit(data.cor)

  pcor_result <- pcor.test(
    x = data.cor$age_acceleration,
    y = data.cor[cor],
    z = data.cor[, c("sex", "age_i2", "education_class",
                      'Mixed' , 'Asian', 'Black', 'Chinese',
                      'others')],
    method = "pearson"
  )
  pcor_result['R2'] <- pcor_result$estimate * pcor_result$estimate
  row.names(pcor_result) <- cor
  pcor_results <- rbind(pcor_results, pcor_result)
}


# lifestyle-disease Mediation======

boot.med <- function(data, mediators, exposure){
  data <- data[sample(1:nrow(data), replace = T), ]
  
  lm_formula <- as.formula(paste(mediators, "~ ", value1, " + sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  alpha.temp <- coefficients(lm(lm_formula, data))[2] 
  
  cox_formula1 <- as.formula(paste("Surv(time, status == 1) ~", mediators, "+ ", exposure, " + 
                                    sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  beta.temp <- coefficients(coxph(cox_formula1, data, method = "breslow"))[1]
  
  cox_formula2 <- as.formula(paste("Surv(time, status == 1) ~ ", exposure, " +", mediators, "+  
                                   sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  c_prime.temp <- coefficients(coxph(cox_formula2, data, method = "breslow"))[1]
  
  cox_formula3 <- as.formula(paste("Surv(time, status == 1) ~ ", exposure, " +  
                                   sex + age_i2 + p53_time + Mixed + Asian + Black + Chinese + others + education_class"))
  c.temp <- coefficients(coxph(cox_formula3, data, method = "breslow"))[1]
  
  IE1.l <- alpha.temp * beta.temp
  IE2.l <- c.temp - c_prime.temp
  DE.l <- c_prime.temp
  TOT.l <- c.temp
  
  results <- c(IE1.l, IE2.l, DE.l, TOT.l)
  return(results)
}
for (life in lifes) {
  result.m <- matrix(nrow = length(diseases), ncol = 18)
  result.m <- data.frame(result.m)
  rownames(result.m) <- diseases
  
  for (disease in diseases) {
    time <- paste0(disease, '_time')
    
    # Bootstrap nums
    G <- 1000
    
    # Bootstrap
    med.boot.cox <- mclapply(1:G, FUN = function(i) boot.med(data.md, 'age_acceleration', life), mc.cores = num_cores)

    IE1.cox <- unlist(lapply(med.boot.cox, '[[', 1)) 
    IE2.cox <- unlist(lapply(med.boot.cox, '[[', 2))
    DE.cox <- unlist(lapply(med.boot.cox, '[[', 3)) 
    TOT.cox <- unlist(lapply(med.boot.cox, '[[', 4)) 
    
    calc_pval <- function(x) {
      pval <- mean(x < 0)
      if (pval > 0.5) {
        pval <- 1 - pval
      }
      return(2 * pval)
    }
    
    dir.cox <- mean(DE.cox) # IE c-c' method
    ci.dir.cox <- as.table(quantile(DE.cox, c(0.025, 0.975))) # 95% CI's c-c'
    
    tot.cox <- mean(TOT.cox) # IE c-c' method
    ci.tot.cox <- as.table(quantile(TOT.cox, c(0.025, 0.975)))
    
    ind.ab.cox <- mean(IE1.cox) # IE ab method
    ind.ccp.cox <- mean(IE2.cox) # IE c-c' method
    ci.ab.cox <- as.table(quantile(IE1.cox, c(0.025, 0.975))) # 95% CI's ab
    ci.ccp.cox <- as.table(quantile(IE2.cox, c(0.025, 0.975))) # 95% CI's c-c'
    
    PM_1 <- ind.ab.cox / tot.cox  
    PM_2 <- ind.ccp.cox / tot.cox
    
    pval_IE1 <- calc_pval(IE1.cox)
    pval_IE2 <- calc_pval(IE2.cox)
    pval_DE <- calc_pval(DE.cox)
    pval_TOT <- calc_pval(TOT.cox)
    
    result.m[which(rownames(result.m) == disease), ] <- c(dir.cox, ci.dir.cox[1], ci.dir.cox[2], pval_DE,
                                                          tot.cox, ci.tot.cox[1], ci.tot.cox[2],pval_TOT,
                                                          ind.ab.cox, ci.ab.cox[1], ci.ab.cox[2], pval_IE1,
                                                          ind.ccp.cox, ci.ccp.cox[1], ci.ccp.cox[2], pval_IE2,
                                                          PM_1, PM_2)
  }
  
  colnames(result.m) <- c('dir.cox', 'ci.dir.cox_2.5', 'ci.dir.cox_97.5','dir.p',
                          'tot.cox', 'ci.tot.cox_2.5', 'ci.tot.cox_97.5', 'tot.p',
                          'ind.ab.cox', 'ci.ab.cox_2.5', 'ci.ab.cox_97.5', 'ind.ab.p',
                          'ind.ccp.cox', 'ci.ccp.cox_2.5', 'ci.ccp.cox_97.5','ind.ccp.p',
                          'PM_1','PM_2')
}
# ====


for (life in lifes) {
  for (disease in diseases) {
    data.cox <- data[, c(disease, life, 'PRS', 'sex', 'age_i2','p53_time', 
                         'education_class')]
    data.cox <- na.omit(data.cox)
    
    data.cox$PRS <- as.factor(data.cox$PRS)
    data.cox[[life]] <- as.factor(data.cox[[life]])
    
    
    lm_formula <- as.formula(
      paste0('age_acceleration ~ ', life,'*PRS + age_i2 + p53_time + sex + education_class + Mixed + Asian + Black + Chinese + others')
    )
    
    model <- glm(lm_formula,data=data.cox)
    
    result <- data.frame(summary(model)$coefficients)
    
    coefficients <- coef(model)
    std_errors <- sqrt(diag(vcov(model)))
    
    result['OR'] <- exp(coefficients)
    result['lower.95'] <- exp(coefficients - 1.96 * std_errors)
    result['upper.95'] <- exp(coefficients + 1.96 * std_errors)

    
    #delta
    out <- interactionR(model,
                        #delta
                        exposure_names =c("PRS", life),
                        ci.type ="delta", ci.level = 0.95,em = F, recode = F)
    result <- out$dframe
    result$life <- life
    result$disease <- disease 

  }
}

# Multiplicative interaction stratified analysis
low_PRS <- subset(data, PRS == 0)
high_PRS <- subset(data, PRS == 1)

for (life in lifes) {
  lm_formula <- as.formula(
    paste0('age_acceleration ~ ',life,' + age_i2 + p53_time + sex + education_class + Mixed + Asian + Black + Chinese + others')
  )
  model <- glm(lm_formula,data=low_PRS)
  
  result <- data.frame(summary(model)$coefficients)
  coefficients <- coef(model)
  std_errors <- sqrt(diag(vcov(model)))
  
  
  result['lower.95'] <- coefficients - 1.96 * std_errors
  result['upper.95'] <- coefficients + 1.96 * std_errors


  
  model <- glm(lm_formula,data=high_PRS)
  
  result <- data.frame(summary(model)$coefficients)

  coefficients <- coef(model)
  std_errors <- sqrt(diag(vcov(model)))
  
  result['lower.95'] <- coefficients - 1.96 * std_errors
  result['upper.95'] <- coefficients + 1.96 * std_errors

}

# Additive interaction stratified analysis
for (life in lifes) {
  
  df <- subset(data, !is.na(data[[life]]))
  df$groups <- paste0(df$PRS, df[[life]], sep = "")
  df <- df %>%
    mutate(groups = case_when(
      groups == "00" ~ 1,
      groups == "01" ~ 2,
      groups == "10" ~ 3,
      groups == "11" ~ 4
    ))
  df$groups <- factor(df$groups)
  
  table(df[[life]], useNA='always')
  
  lm_formula <- as.formula(
    paste0('age_acceleration ~ ','groups + age_i2 + p53_time + sex + education_class + Mixed + Asian + Black + Chinese + others')
  )
  
  model <- glm(lm_formula,data=df)
  
  result <- data.frame(summary(model)$coefficients)
  
  coefficients <- coef(model)
  std_errors <- sqrt(diag(vcov(model)))
  
  result['lower.95'] <- exp(coefficients - 1.96 * std_errors)
  result['upper.95'] <- exp(coefficients + 1.96 * std_errors)
}

# scPagwas ====
library(scPagwas)
library(Seurat)

load('genes.by.pathway_kegg.RData')
load('chrom_ld_new.RData')
load('block_annotation_hg37_new.RData')

Pagwas <- scPagwas_main(Pagwas =NULL,
                        gwas_data="gwas.txt", # GWAS   
                        Single_data ='brain_umap.rds', # Seurat 
                        output.prefix="brain", 
                        output.dirs="scPagwastest_output",    
                        Pathway_list=genes.by.pathway_kegg,  
                        n.cores=10,
                        assay="RNA",    
                        singlecell=T,    
                        iters_singlecell = 100, 
                        celltype=T,    
                        block_annotation = block_annotation_hg37_new,    
                        chrom_ld = chrom_ld_new)

save(Pagwas,file="brain_scPagwas.RData")


# KEGG-GO====
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(reactome.db)

data <- fread('data.csv')
data <- data %>%
  filter(p_value < (0.05 / nrow(data)))
interactors <- data$protein 
id <- bitr(interactors, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# go
go <- enrichGO(gene          = id$ENTREZID,
               OrgDb         = "org.Hs.eg.db",
               ont           = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05)

GO <- go@result 

# kegg
R.utils::setOption("clusterProfiler.download.method","auto")

kegg <- enrichKEGG(
  gene = id$ENTREZID,
  pvalueCutoff = 0.05,
  organism = "hsa")  
KEGG <- kegg@result
