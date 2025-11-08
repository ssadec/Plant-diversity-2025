
prefix <- '/Users/yusha/'

filenames <- dir(path = prefix, pattern = '^0903Paired data for herbivore and crop.*\\.xlsx$')


target_cols_1 <- c("yi1", "vi1", "Trtmnt.s1", "Cntrl.s1", "cv2_trea_new1", "cv2_cont_new1")
target_cols_2 <- c("yi2", "vi2", "Trtmnt.s2", "Cntrl.s2", "cv2_trea_new2", "cv2_cont_new2")

dt_total0 <- NULL
nms0 <- NULL

for (i in seq_along(filenames)) {
  file_path <- file.path(prefix, filenames[i])

  
  matched <- regmatches(filenames[i], regexec("Paired data for ([A-Za-z0-9_]+) and ([A-Za-z0-9_]+)", filenames[i]))
  if (length(matched[[1]]) < 2) {
    warning(paste("file can not extract two categories，skip：", filenames[i]))
    next
  }
  
  category <- tolower(matched[[1]][2:3])
  sheets <- category
  
  dt_new0 <- NULL
  
  for (j in 1:2) {
   
    dt0 <- tryCatch({
      read_excel(file_path, sheet = sheets[j], range = "A1:AG417")
    }, error = function(e) {
      warning(paste("wrong：", sheets[j], "in file", filenames[i]))
      return(NULL)
    })
    if (is.null(dt0)) next
    
   
    colnames(dt0) <- make.names(abbreviate(colnames(dt0), minlength = 4))
    dt0 <- subset(dt0, !is.na(CONTROL_M))
    
    dt0$CONTROL_SD <- replace(dt0$CONTROL_SD, round(dt0$CONTROL_M / dt0$CONTROL_SD, 5) == 10, NA)
    dt0$TREATMENT_SD <- replace(dt0$TREATMENT_SD, is.na(dt0$CONTROL_SD), NA)
    
    
    ESData <- dt0 %>%
      mutate(cv_Control = na_if(CONTROL_SD / CONTROL_M, Inf),
             cv_Treatment = na_if(TREATMENT_SD / TREATMENT_M, Inf)) %>%
      cv_avg(x = TREATMENT_M, sd = TREATMENT_SD, n = TREATMENT_R, group = S_ID, label = "1") %>%
      cv_avg(x = CONTROL_M, sd = CONTROL_SD, n = TREATMENT_R, group = S_ID, label = "2")
    
  
    if (j == 1) {
      ESData <- ESData %>%
        mutate(cv2_trea_new1 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
               cv2_cont_new1 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2),
               yi1 = lnrr_laj(TREATMENT_M, CONTROL_M, b_CV2_1, b_CV2_2, TREATMENT_R, TREATMENT_R),
               vi1 = v_lnrr_laj(b_CV2_1, TREATMENT_R, b_CV2_2, TREATMENT_R),
               Trtmnt.s1 = if_else(is.na(TREATMENT_SD), cv2_trea_new1 * TREATMENT_M, TREATMENT_SD),
               Cntrl.s1  = if_else(is.na(CONTROL_SD), cv2_cont_new1 * CONTROL_M, CONTROL_SD))
    } else if (j == 2) {
      ESData <- ESData %>%
        mutate(cv2_trea_new2 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
               cv2_cont_new2 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2),
               yi2 = lnrr_laj(TREATMENT_M, CONTROL_M, b_CV2_1, b_CV2_2, TREATMENT_R, TREATMENT_R),
               vi2 = v_lnrr_laj(b_CV2_1, TREATMENT_R, b_CV2_2, TREATMENT_R),
               Trtmnt.s2 = if_else(is.na(TREATMENT_SD), cv2_trea_new2 * TREATMENT_M, TREATMENT_SD),
               Cntrl.s2  = if_else(is.na(CONTROL_SD), cv2_cont_new2 * CONTROL_M, CONTROL_SD))
    }
    
 
    if (is.null(dt_new0)) {
      dt_new0 <- ESData
    } else {
      if (j == 2) {
        dt_new0 <- cbind(dt_new0, ESData[, target_cols_2, drop = FALSE])
      }
    }
  }
  

  if (!is.null(dt_new0)) {
    S_ID <- table(dt_new0$S_ID)
    dt_new0$EsID <- unlist(lapply(S_ID, function(k) seq_len(k)))
    dt_new0$added_plant_species <- log2(dt_new0$NCP_N)
    dt_new0$antagonist <- category[1]
    
    if(all(c("yi1", "yi2") %in% colnames(dt_new0))) {
      dt_new0 <- dt_new0 %>%
        filter(is.finite(yi1), !is.na(yi1),
               is.finite(yi2), !is.na(yi2))
    } else {
      warning("not full data")
    }
    
    
    

    dt_new0 <- dt_new0 %>%
      filter(is.finite(yi1), !is.na(yi1),
             is.finite(yi2), !is.na(yi2))
    
    dt_total0 <- rbind(dt_total0, dt_new0)
    valid_cat <- category[!is.na(category)]
    nms0 <- rbind(nms0, c(paste(category, collapse = "_vs_"), category))
    
    
  }
}



View(dt_total0)
View(nms0)
View(dt_new0)





write.table(nms0, paste0(prefix, "paired_names.txt"), sep = "\t", row.names = T)
write.table(dt_total0, paste0(prefix, "herbivore_vs_crop.txt"), sep = "\t", row.names = T)



####Table S11####
res <- NULL
for (i in 1:nrow(paired_names)) {
  
  dat<-read.table('herbivore_vs_crop.txt', sep = "\t", header=TRUE)
  
  
  for (j in 1:2) {
    if (j == 1) {
      V1 <- tcrossprod(sqrt(dat$Cntrl.s1))/sqrt(outer(dat$CONTROL_R * (dat$CONTROL_M)^2, dat$CONTROL_R * (dat$CONTROL_M)^2))+
        diag((dat$Trtmnt.s1)/(dat$TREATMENT_R * (dat$TREATMENT_M)^2))
      V_ind <- pmax(as.matrix(crossprod(lFormula(yi1~1+(1|S_ID/CONTROL_M), dat)$reTrms$Zt)) - 1, 0)
      V1 <- V1*V_ind
      V1a<-make.positive.definite(V1)
      dat.fit <- try(rma.mv(yi1, V1a, random = list(~1|S_ID), data = dat, test = "t",control = list(optimizer = "optim")))
      
    } else if  (j == 2) {
      V2 <- tcrossprod(sqrt(dat$Cntrl.s2))/sqrt(outer(dat$CONTROL_R * (dat$CONTROL_M)^2, dat$CONTROL_R * (dat$CONTROL_M)^2))+
        diag((dat$Trtmnt.s2)/(dat$TREATMENT_R * (dat$TREATMENT_M)^2))
      V_ind <- pmax(as.matrix(crossprod(lFormula(yi2~1+(1|S_ID/CONTROL_M), dat)$reTrms$Zt)) - 1, 0)
      V2 <- V2*V_ind
      V2a<-make.positive.definite(V2)
      dat.fit <- try(rma.mv(yi2, V2a, random = list(~1|S_ID), data = dat, test = "t",control = list(optimizer = "optim")))
      
    }
    res <- rbind(res, c(paired_names[i, j+1], nrow(dat), length(unique(dat$S_ID)), dat.fit$b,
                        coef(summary(dat.fit))$tval, coef(summary(dat.fit))$pval, coef(summary(dat.fit))$df,
                        dat.fit$ci.lb, dat.fit$ci.ub))
  }
}
res <- data.frame(paired_names[, 1],
                  res)
colnames(res) <- c("Trophic interaction",
                   "Category",
                   "Number of observations",
                   "Number of studies",
                   "Effect size",
                   "t-value",
                   "P-value",
                   "df",
                   "Lower Bound of Confidence Interval",
                   "Upper Bound of Confidence  Interval")
View(res)







####Table S12####


library(nlme)
library(rsq)
library(openxlsx)

library(readxl)
library(nlme)
library(piecewiseSEM)   # rsq.lmm
library(xlsx)          

library(nlme)
library(piecewiseSEM)
library(dplyr)
library(readxl)


rm(list = ls()) 
.rs.restartR()   



prefix <- '/Users/yusha/'
paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE)

res <- NULL

for (i in 1:nrow(paired_names)) {
  
  
 dat <- read_xlsx(paste0(prefix, "0902-V15.xlsx"))
  View(dat)
  dat$S_ID <- as.factor(dat$S_ID)
  
#binary variable
  dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, ifelse(dat$NCP_N >= 2, 1, NA))
  

  safe_lme <- function(formula, data, random, weights = NULL) {
    tryCatch({
      lme(formula, data = data, random = random, weights = weights)
    }, error = function(e) NULL)
  }
  

  fit_yi2 <- safe_lme(yi2 ~ -1 + yi1, data = dat,
                      random = ~ -1 + yi1 | S_ID,
                      weights = varFixed(~vi2))
  dat$yi2_resid <- if(!is.null(fit_yi2)) residuals(fit_yi2, type = "pearson") else NA
  

  r1 <- tryCatch(residuals(safe_lme(yi1 ~ 1, data = dat, random = ~1 | S_ID, weights = varFixed(~vi1))), 
                 error = function(e) rep(NA, nrow(dat)))
  r2 <- tryCatch(residuals(safe_lme(yi2 ~ 1, data = dat, random = ~1 | S_ID, weights = varFixed(~vi2))), 
                 error = function(e) rep(NA, nrow(dat)))
  dat.residual <- data.frame(r1 = r1, r2 = r2, S_ID = dat$S_ID)

  
  mod_list <- list(
    list(name = "yi1_on_yi2", formula = r1 ~ -1 + r2, random = ~ -1 + r2 | S_ID, data = dat.residual, weights = NULL),
    list(name = "highdiv_on_yi1", formula = yi1 ~ plant_diversity, random = ~1 | S_ID, data = dat, weights = varFixed(~vi1)),
    list(name = "highdiv_on_yi2", formula = yi2 ~ plant_diversity, random = ~1 | S_ID, data = dat, weights = varFixed(~vi2)),
    list(name = "added_on_yi1", formula = yi1 ~ 1 + added_plant_species, random = ~1 | S_ID, data = dat, weights = varFixed(~vi1)),
    list(name = "added_on_yi2", formula = yi2 ~ 1 + added_plant_species, random = ~1 | S_ID, data = dat, weights = varFixed(~vi2))
  )
  

  for (mod in mod_list) {
    fit <- safe_lme(mod$formula, data = mod$data, random = mod$random, weights = mod$weights)
    if (is.null(fit)) next 
    
    summ <- summary(fit)$tTable

    idx <- intersect(c("plant_diversity", "added_plant_species", "r2", "(Intercept)"), rownames(summ))
    if (length(idx) == 0) next
    
    idx <- idx[1] 
    
    est  <- summ[idx, "Value"]
    se   <- summ[idx, "Std.Error"]
    tval <- summ[idx, "t-value"]
    pval <- summ[idx, "p-value"]
    
    if (any(!is.finite(c(est, se, tval, pval)))) next
    
    df <- tryCatch(nrow(mod$data) - length(fixef(fit)), error = function(e) NA)
    r2 <- tryCatch(rsq.lmm(fit)$model, error = function(e) NA)
    
    ci_l <- est - se * qt(0.975, df)
    ci_u <- est + se * qt(0.975, df)
    
    res <- rbind(res, data.frame(
      dataset = paired_names[i, 1],
      model = mod$name,
      predictor = idx,
      estimate = est,
      se = se,
      t = tval,
      p = pval,
      df = df,
      CI_lower = ci_l,
      CI_upper = ci_u,
      R2 = r2
    ))
  }
}


View(res)

 write.xlsx(res, paste0(prefix, "SuppTab12.xlsx"), row.names = FALSE)



 #####STable 13####
#####sub-categories####
res <- NULL

for (i in 1:nrow(paired_names)) {
  

  dat <- read_xlsx(paste0(prefix, "0902-V15.xlsx"))
  

  dat$S_ID <- as.factor(dat$S_ID)
  
  categories <- c("Herbaceous_plant", "Woody_plant") #CP_C3
  #categories <- c("Generalist", "Speciailist") #Mn__
  #categories <- c("Intercropping", "Sown_field_margins","Cover_cropping") #NCP_P
  #categories <- c("Food_crop", "Cash_crop") #CP_Ct
  #categories <- c("Temperate", "Tropic") #Zone
 #categories <- unique(na.omit(dat$Mn__))
  
  for (category in categories) {

    dat_cat <- dat[dat$CP_C3 == category, ]
    
    

    if (nrow(dat_cat) == 0) next
    
    # binary
    dat_cat$plant_diversity <- ifelse(dat_cat$NCP_N == 1, 0, ifelse(dat_cat$NCP_N >= 2, 1, NA))
    
  
    safe_lme <- function(formula, data, random, weights = NULL) {
      tryCatch({
        lme(formula, data = data, random = random, weights = weights)
      }, error = function(e) NULL)
    }
    
 
    fit_yi2 <- safe_lme(yi2 ~ -1 + yi1, data = dat_cat,
                        random = ~ -1 + yi1 | S_ID,
                        weights = varFixed(~vi2))
    dat_cat$yi2_resid <- if(!is.null(fit_yi2)) residuals(fit_yi2, type = "pearson") else NA
    
    
    r1 <- tryCatch(residuals(safe_lme(yi1 ~ 1, data = dat_cat, random = ~1 | S_ID, weights = varFixed(~vi1))), 
                   error = function(e) rep(NA, nrow(dat_cat)))
    r2 <- tryCatch(residuals(safe_lme(yi2 ~ 1, data = dat_cat, random = ~1 | S_ID, weights = varFixed(~vi2))), 
                   error = function(e) rep(NA, nrow(dat_cat)))
    dat.residual <- data.frame(r1 = r1, r2 = r2, S_ID = dat_cat$S_ID)
    

    mod_list <- list(
      list(name = "yi1_on_yi2", formula = r1 ~ -1 + r2, random = ~ -1 + r2 | S_ID, data = dat.residual, weights = NULL),
      list(name = "highdiv_on_yi1", formula = yi1 ~ plant_diversity, random = ~1 | S_ID, data = dat_cat, weights = varFixed(~vi1)),
      list(name = "highdiv_on_yi2", formula = yi2 ~ plant_diversity, random = ~1 | S_ID, data = dat_cat, weights = varFixed(~vi2)),
      list(name = "added_on_yi1", formula = yi1 ~ 1 + added_plant_species, random = ~1 | S_ID, data = dat_cat, weights = varFixed(~vi1)),
      list(name = "added_on_yi2", formula = yi2 ~ 1 + added_plant_species, random = ~1 | S_ID, data = dat_cat, weights = varFixed(~vi2))
    )
    

    for (mod in mod_list) {
      fit <- safe_lme(mod$formula, data = mod$data, random = mod$random, weights = mod$weights)
      if (is.null(fit)) next  
      
      summ <- summary(fit)$tTable
      # 选择目标系数
      idx <- intersect(c("plant_diversity", "added_plant_species", "r2", "(Intercept)"), rownames(summ))
      if (length(idx) == 0) next
      
      idx <- idx[1] 
      
      est  <- summ[idx, "Value"]
      se   <- summ[idx, "Std.Error"]
      tval <- summ[idx, "t-value"]
      pval <- summ[idx, "p-value"]
      
      if (any(!is.finite(c(est, se, tval, pval)))) next
      
      df <- tryCatch(nrow(mod$data) - length(fixef(fit)), error = function(e) NA)
      r2 <- tryCatch(rsq.lmm(fit)$model, error = function(e) NA)
      
      ci_l <- est - se * qt(0.975, df)
      ci_u <- est + se * qt(0.975, df)
      
      res <- rbind(res, data.frame(
        dataset = paired_names[i, 1],
        model = mod$name,
        predictor = idx,
        estimate = est,
        se = se,
        t = tval,
        p = pval,
        df = df,
        CI_lower = ci_l,
        CI_upper = ci_u,
        R2 = r2,
        category = category
      ))
    }
  }
}
View(res)
