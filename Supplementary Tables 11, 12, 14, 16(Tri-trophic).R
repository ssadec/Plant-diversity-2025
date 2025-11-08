

prefix <-'/Users/yusha/'

filenames <- dir(path = prefix, pattern = "^0901Paired data for enemies and herbivore and crop.*\\.xlsx$")


target_cols_1 <- c("yi1", "vi1", "Trtmnt.s1", "Cntrl.s1", "cv2_trea_new1", "cv2_cont_new1")
target_cols_2 <- c("yi2", "vi2", "Trtmnt.s2", "Cntrl.s2", "cv2_trea_new2", "cv2_cont_new2")
target_cols_3 <- c("yi3", "vi3", "Trtmnt.s3", "Cntrl.s3", "cv2_trea_new3", "cv2_cont_new3")

dt_total0 <- NULL
nms0 <- NULL

for (i in seq_along(filenames)) {
  file_path <- file.path(prefix, filenames[i])
  

  matched <- regmatches(filenames[i], regexec("Paired data for ([A-Za-z0-9_]+) and ([A-Za-z0-9_]+) and ([A-Za-z0-9_]+)", filenames[i]))
  if (length(matched[[1]]) < 4) {
    warning(paste("can not withdraw three different categories，skip：", filenames[i]))
    next
  }
  
  category <- tolower(matched[[1]][2:4])
  sheets <- category 
  
  dt_new0 <- NULL
  
  for (j in 1:3) {
   
    dt0 <- tryCatch({
      read_excel(file_path, sheet = sheets[j], range = "A1:AG417")
    }, error = function(e) {
      warning(paste("warning：", sheets[j], "in file", filenames[i]))
      return(NULL)
    })
    if (is.null(dt0)) next
    
    # preprosedure
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
    } else if (j == 3) {
      ESData <- ESData %>%
        mutate(cv2_trea_new3 = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
               cv2_cont_new3 = if_else(is.na(cv_Control), b_CV2_2, cv_Control^2),
               yi3 = lnrr_laj(TREATMENT_M, CONTROL_M, b_CV2_1, b_CV2_2, TREATMENT_R, TREATMENT_R),
               vi3 = v_lnrr_laj(b_CV2_1, TREATMENT_R, b_CV2_2, TREATMENT_R),
               Trtmnt.s3 = if_else(is.na(TREATMENT_SD), cv2_trea_new3 * TREATMENT_M, TREATMENT_SD),
               Cntrl.s3  = if_else(is.na(CONTROL_SD), cv2_cont_new3 * CONTROL_M, CONTROL_SD))
    }
    
   
    if (is.null(dt_new0)) {
      dt_new0 <- ESData
    } else {
      if (j == 2) {
        dt_new0 <- cbind(dt_new0, ESData[, target_cols_2, drop = FALSE])
      } else if (j == 3) {
        dt_new0 <- cbind(dt_new0, ESData[, target_cols_3, drop = FALSE])
      }
    }
  }
  
 
  if (!is.null(dt_new0)) {
    S_ID <- table(dt_new0$S_ID)
    dt_new0$EsID <- unlist(lapply(S_ID, function(k) seq_len(k)))
    dt_new0$added_plant_species <- log2(dt_new0$NCP_N)
    dt_new0$antagonist <- category[1]
    
    if(all(c("yi1", "yi2", "yi3") %in% colnames(dt_new0))) {
      dt_new0 <- dt_new0 %>%
        filter(is.finite(yi1), !is.na(yi1),
               is.finite(yi2), !is.na(yi2),
               is.finite(yi3), !is.na(yi3))
    } else {
      warning("effect sizes are not completed，skip")
    }
    
    
    
    
    dt_new0 <- dt_new0 %>%
      filter(is.finite(yi1), !is.na(yi1),
             is.finite(yi2), !is.na(yi2),
             is.finite(yi3), !is.na(yi3))
    
    dt_total0 <- rbind(dt_total0, dt_new0)
    nms0 <- rbind(nms0, c(paste(category, collapse = "_vs_"), category))
    
  }
}



View(dt_total0)
View(nms0)
View(dt_new0)





write.table(nms0, paste0(prefix, "paired_names.txt"), sep = "\t", row.names = T)
write.table(dt_total0, paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", row.names = T)



####Table S11####
paired_names <- read.table(paste(prefix, "paired_names.txt", sep = ""), sep = "\t", header = T)
View(paired_names)

res <- NULL
for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste(prefix, paired_names[i, 1], ".txt", sep = ""), sep = "\t", header= T)
  

  
  for (j in 1:3) {
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
      
    } else if  (j == 3){
      V3 <- tcrossprod(sqrt(dat$Cntrl.s3))/sqrt(outer(dat$CONTROL_R * (dat$CONTROL_M)^2, dat$CONTROL_R * (dat$CONTROL_M)^2))+
        diag((dat$Trtmnt.s3)/(dat$TREATMENT_R * (dat$TREATMENT_M)^2))
      V_ind <- pmax(as.matrix(crossprod(lFormula(yi3~1+(1|S_ID/CONTROL_M), dat)$reTrms$Zt)) - 1, 0)
      V3 <- V3*V_ind
      V3a<-make.positive.definite(V3)
      dat.fit <- try(rma.mv(yi3, V3a, random = list(~1|S_ID), data = dat, test = "t",control = list(optimizer = "optim")))
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
#xlsx::write.xlsx(res, paste(path, "SuppTab8b.xlsx", sep=""), row.names = F)
View(res)
write.csv(res, '/Users/yusha/SuppTab11.csv')


########Table S12########


library(nlme)
library(rsq)
library(openxlsx)

res <- NULL  

for (i in 1:nrow(paired_names)) {  
  dat <- read.table(paste(prefix, paired_names[i, 1], ".txt", sep = ""), sep = "\t", header= T)
  
   
  dat$S_ID <- as.factor(dat$S_ID)  
  
#binary
  dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, 
                               ifelse(dat$NCP_N >= 2, 1, NA))

  table(dat$plant_diversity, dat$NCP_N, useNA = "always")
  # ================
 
  dat$yi2_resid <- residuals(  
    lme(yi2 ~ -1 + yi1, data = dat,  
        random = ~ -1 + yi1 | S_ID,  
        weights = varFixed(~vi2),  
        control = lmeControl(sigma = 1)),  
    type = "pearson"  
  )  
  
  dat$yi3_resid <- residuals(  
    lme(yi3 ~ -1 + yi2, data = dat,  
        random = ~ -1 + yi2 | S_ID,  
        weights = varFixed(~vi3),  
        control = lmeControl(sigma = 1)),  
    type = "pearson"  
  )  
  
  dat.residual <- data.frame(  
    r1 = residuals(lme(yi1 ~ 1, data = dat, random = ~1 | S_ID, weights = varFixed(~vi1))),  
    r2 = residuals(lme(yi2 ~ 1, data = dat, random = ~1 | S_ID, weights = varFixed(~vi2))),  
    r3 = residuals(lme(yi3 ~ 1, data = dat, random = ~1 | S_ID, weights = varFixed(~vi3))),  
    S_ID = dat$S_ID  
  )  
  

  mod_list <- list(  
    list(name = "yi1_on_yi2", fit = lme(r1 ~ -1 + r2, data = dat.residual, random = ~ -1 + r2 | S_ID)),  
    list(name = "yi2_on_yi3", fit = lme(r2 ~ -1 + r3, data = dat.residual, random = ~ -1 + r3 | S_ID)),  
 
    list(name = "highdiv_on_yi1", fit = lme(yi1 ~ plant_diversity, data = dat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi1))),  
    list(name = "highdiv_on_yi2", fit = lme(yi2 ~ plant_diversity, data = dat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi2))),  
    list(name = "highdiv_on_yi3", fit = lme(yi3 ~ plant_diversity, data = dat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi3))),  
    list(name = "added_on_yi1", fit = lme(yi1 ~ 1 + added_plant_species, data = dat, 
                                          random = ~1 | S_ID, 
                                          weights = varFixed(~vi1))),  
    list(name = "added_on_yi2", fit = lme(yi2 ~ 1 + added_plant_species, data = dat, 
                                          random = ~1 | S_ID, 
                                          weights = varFixed(~vi2))),  
    list(name = "added_on_yi3", fit = lme(yi3 ~ 1 + added_plant_species, data = dat, 
                                          random = ~1 | S_ID, 
                                          weights = varFixed(~vi3)))  
  )  
  
  for (mod in mod_list) {  
    fit <- mod$fit  
    summ <- summary(fit)$tTable 
    
   
    if ("plant_diversity" %in% rownames(summ)) {  
      idx <- "plant_diversity"  
    } else if ("added_plant_species" %in% rownames(summ)) {  
      idx <- "added_plant_species"  
    } else if ("r2" %in% rownames(summ)) {  
      idx <- "r2"  
    } else if ("r3" %in% rownames(summ)) {  
      idx <- "r3"  
    } else if ("(Intercept)" %in% rownames(summ)) {  
      idx <- "(Intercept)"  
    } else {  
      next 
    }  
    

    est  <- summ[idx, "Value"]        
    se   <- summ[idx, "Std.Error"]    
    tval <- summ[idx, "t-value"]      
    pval <- summ[idx, "p-value"]      
    
  
    if (any(!is.finite(c(est, se, tval, pval)))) next  
    
 
    df <- tryCatch(nrow(dat) - length(fixef(fit)), error = function(e) NA)  
    

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



########Table S14########
####different subcategories####
res <- NULL  
for (i in 1:nrow(paired_names)) {
  

  dat <- read.table(paste(prefix, paired_names[i, 1], ".txt", sep = ""), sep = "\t", header= T)

  

  dat$S_ID <- as.factor(dat$S_ID)
  
  dat$plant_diversity <- ifelse(dat$NCP_N == 1, 0, 
                               ifelse(dat$NCP_N >= 2, 1, NA))
  
  

  categories <- c("Intercropping","Cover_cropping") #NCP_P,"Cover_cropping"
  #categories <- c("Annual", "Perennial") #CP_C2
  #categories <- c("Food_crop", "Cash_crop") #CP_Ct
  #categories <- c("Herbaceous_plant", "Woody_plant") #CP_C3
  #categories <- c("Temperate", "Tropic") #Zone
  #categories <- c("Generalist", "Specialist") #Mn__
  #categories <- c("Generalist", "Specialist") #NE_c
  
  #categories <- unique(na.omit(dat$NCP_P))
  
  for (category in categories) {

    dat_cat <- dat[!is.na(dat$NCP_P) & dat$NCP_P == category, ]
    
    
    
    if (nrow(dat_cat) < 3) next 
  
    
    yi2_resid <- residuals(  
      lme(yi2 ~ -1 + yi1, data = dat_cat,  
          random = ~ -1 + yi1 | S_ID,  
          weights = varFixed(~vi2),  
          control = lmeControl(sigma = 1)),  
      type = "pearson"  
    )  
    
  
    yi3_resid <- residuals(  
      lme(yi3 ~ -1 + yi2, data = dat_cat,  
          random = ~ -1 + yi2 | S_ID,  
          weights = varFixed(~vi3),  
          control = lmeControl(sigma = 1)),  
      type = "pearson"  
    )  
    

    dat.residual <- data.frame(  
      r1 = residuals(lme(yi1 ~ 1, data = dat_cat, random = ~1 | S_ID, weights = varFixed(~vi1))),  
      r2 = residuals(lme(yi2 ~ 1, data = dat_cat, random = ~1 | S_ID, weights = varFixed(~vi2))),  
      r3 = residuals(lme(yi3 ~ 1, data = dat_cat, random = ~1 | S_ID, weights = varFixed(~vi3))),  
      S_ID = dat_cat$S_ID  
    )  
    
 
    mod_list <- list(  
      list(name = "yi1_on_yi2", fit = lme(r1 ~ -1 + r2, data = dat.residual, random = ~ -1 + r2 | S_ID)),  
      list(name = "yi2_on_yi3", fit = lme(r2 ~ -1 + r3, data = dat.residual, random = ~ -1 + r3 | S_ID)),  

      list(name = "highdiv_on_yi1", fit = lme(yi1 ~ plant_diversity, data = dat_cat, 
                                              random = ~1 | S_ID, 
                                              weights = varFixed(~vi1))),  
      list(name = "highdiv_on_yi2", fit = lme(yi2 ~ plant_diversity, data = dat_cat, 
                                              random = ~1 | S_ID, 
                                              weights = varFixed(~vi2))),  
      list(name = "highdiv_on_yi3", fit = lme(yi3 ~ plant_diversity, data = dat_cat, 
                                              random = ~1 | S_ID, 
                                              weights = varFixed(~vi3))),  
      list(name = "added_on_yi1", fit = lme(yi1 ~ 1 + added_plant_species, data = dat_cat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi1))),  
      list(name = "added_on_yi2", fit = lme(yi2 ~ 1 + added_plant_species, data = dat_cat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi2))),  
      list(name = "added_on_yi3", fit = lme(yi3 ~ 1 + added_plant_species, data = dat_cat, 
                                            random = ~1 | S_ID, 
                                            weights = varFixed(~vi3)))  
    )  
    
    for (mod in mod_list) {  
      fit <- mod$fit  
      summ <- summary(fit)$tTable
      

      if ("plant_diversity" %in% rownames(summ)) {  
        idx <- "plant_diversity"  
      } else if ("added_plant_species" %in% rownames(summ)) {  
        idx <- "added_plant_species"  
      } else if ("r2" %in% rownames(summ)) {  
        idx <- "r2"  
      } else if ("r3" %in% rownames(summ)) {  
        idx <- "r3"  
      } else if ("(Intercept)" %in% rownames(summ)) {  
        idx <- "(Intercept)"  
      } else {  
        next 
      }  
      
     
      est  <- summ[idx, "Value"]        
      se   <- summ[idx, "Std.Error"]    
      tval <- summ[idx, "t-value"]      
      pval <- summ[idx, "p-value"]      
      
  
      if (any(!is.finite(c(est, se, tval, pval)))) next  
 
      df <- tryCatch(nrow(dat_cat) - length(fixef(fit)), error = function(e) NA)  
      
 
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


