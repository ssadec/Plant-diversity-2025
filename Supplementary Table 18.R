rm(list = ls())
library(nlme)
library(MuMIn)
library(rsq)

prefix <- '/Users/yusha/'


paired_names <- read.table(paste(prefix, "paired_names.txt", sep = ""), sep = "\t", header = T)

## SuppTab18a-1####
res <- NULL
for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste(prefix, paired_names[i, 1], ".txt", sep = ""), sep = "\t", header= T)
  dat.fit1 <- lme(yi1~1, dat, list(~1|S_ID), 
                  weights = varFixed(~vi1), 
                  control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2~1, dat, list(~1|S_ID),
                  weights = varFixed(~vi2), 
                  control = lmeControl(sigma = 1))
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"),
                             r2 = residuals(dat.fit2, type = "pearson"),
                             vi1 = dat$vi1,
                             vi2 = dat$vi2,
                             S_ID = dat$S_ID)
  yi1_on_yi2 <- lme(r2~-1+r1, dat.residual, list(~-1+r1|S_ID))
  res <- rbind(res, data.frame(paired_names[i, 1],
                               AIC(yi1_on_yi2),
                               AICc(yi1_on_yi2),
                               BIC(yi1_on_yi2),
                               logLik(yi1_on_yi2),
                               rsq.lmm(yi1_on_yi2),
                               fix.empty.names = F))
}

Strat <- c("Global",
           "Herbaceous_plant",
           "Woody_plant",
           "Food_crop",
           "Cash_crop",
           "Generalist",
           "Specialist",
           "Intercropping",
           "Cover_cropping",
           "Sown_field_margins",
           "Natural",
           "Temperate",
           "Tropic")

for (i in Strat) {
  dat <- read.table(paste(prefix, "herbivore_vs_crop.txt", sep = ""), sep = "\t", header= T)
  if (i != "Global") {
    dat <- subset(dat,
                  (CP_C3 == i)| (CP_Ct == i)|(Mn__ == i)|(NCP_P == i)|(Std_ == i)|(Zone == i))
  }
  dat.fit1 <- lme(yi1~1, dat, list(~1|S_ID), 
                  weights = varFixed(~vi1), 
                  control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2~1, dat, list(~1|S_ID),
                  weights = varFixed(~vi2), 
                  control = lmeControl(sigma = 1))
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"),
                             r2 = residuals(dat.fit2, type = "pearson"),
                             vi1 = dat$vi1,
                             vi2 = dat$vi2,
                             S_ID = dat$S_ID)
  yi1_on_yi2 <- lme(r2~-1+r1, dat.residual, list(~-1+r1|S_ID))
  res <- rbind(res, data.frame(i,
                               AIC(yi1_on_yi2),
                               AICc(yi1_on_yi2),
                               BIC(yi1_on_yi2),
                               logLik(yi1_on_yi2),
                               rsq.lmm(yi1_on_yi2),
                               fix.empty.names = F))
}
res <- data.frame(res)
colnames(res) <- c("Item",
                   "AIC",
                   "AICc",
                   "BIC",
                   "logLik",
                   "Model R-square",
                   "Fixed R-square",
                   "Random R-square")
View(res)

write.csv(res, '/Users/yusha/Stable18a-1.csv')




## SuppTab18a-2####
prefix <- '/Users/yusha/'
paired_names <- read.table(paste(prefix, "paired_names.txt", sep = ""), sep = "\t", header = T)
res <- NULL
for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste(prefix, paired_names[i, 1], ".txt", sep = ""), sep = "\t", header= T)
  dat.fit1 <- lme(yi1~1+added_plant_species, dat, list(~1|S_ID), 
                  weights = varFixed(~vi1), 
                  control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2~1+added_plant_species, dat, list(~1|S_ID),
                  weights = varFixed(~vi2), 
                  control = lmeControl(sigma = 1))
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"),
                             r2 = residuals(dat.fit2, type = "pearson"),
                             vi1 = dat$vi1,
                             vi2 = dat$vi2,
                             S_ID = dat$S_ID)
  yi1_on_yi2 <- lme(r2~-1+r1, dat.residual, list(~-1+r1|S_ID))
  res <- rbind(res, data.frame(paired_names[i, 1],
                               AIC(yi1_on_yi2),
                               AICc(yi1_on_yi2),
                               BIC(yi1_on_yi2),
                               logLik(yi1_on_yi2),
                               rsq.lmm(yi1_on_yi2),
                               fix.empty.names = F))
}

Strat <- c("Global",
           "Herbaceous_plant",
           "Woody_plant",
           "Food_crop",
           "Cash_crop",
           "Generalist",
           "Specialist",
           "Intercropping",
           "Cover_cropping",
           "Sown_field_margins",
           "Natural",
           "Temperate",
           "Tropic")


for (i in Strat) {
  dat <- read.table(paste(prefix, "herbivore_vs_crop.txt", sep = ""), sep = "\t", header= T)
  if (i != "Global") {
    dat <- subset(dat,
                  (CP_C3 == i)| (CP_Ct == i)|(Mn__ == i)|(NCP_P == i)|(Std_ == i)|(Zone == i))
  }
  dat.fit1 <- lme(yi1~1+added_plant_species, dat, list(~1|S_ID), 
                  weights = varFixed(~vi1), 
                  control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2~1+added_plant_species, dat, list(~1|S_ID),
                  weights = varFixed(~vi2), 
                  control = lmeControl(sigma = 1))
  dat.residual <- data.frame(r1 = residuals(dat.fit1, type = "pearson"),
                             r2 = residuals(dat.fit2, type = "pearson"),
                             vi1 = dat$vi1,
                             vi2 = dat$vi2,
                             S_ID = dat$S_ID)
  yi1_on_yi2 <- lme(r2~-1+r1, dat.residual, list(~-1+r1|S_ID))
  res <- rbind(res, data.frame(i,
                               AIC(yi1_on_yi2),
                               AICc(yi1_on_yi2),
                               BIC(yi1_on_yi2),
                               logLik(yi1_on_yi2),
                               rsq.lmm(yi1_on_yi2),
                               fix.empty.names = F))
}
res <- data.frame(res)
colnames(res) <- c("Item",
                   "AIC",
                   "AICc",
                   "BIC",
                   "logLik",
                   "Model R-square",
                   "Fixed R-square",
                   "Random R-square")
View(res)

write.csv(res, '/Users/yusha/Stable18a-2.csv')








####STable 18b-1####
library(nlme)
library(MuMIn)
library(rsq)

prefix <- '/Users/yusha/'
paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE)

res <- NULL


for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste0(prefix, paired_names[i, 1], ".txt"), sep = "\t", header = TRUE)
  
  dat.fit1 <- lme(yi1 ~ 1, dat, random = ~1 | S_ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1, dat, random = ~1 | S_ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1, dat, random = ~1 | S_ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(
    r1 = residuals(dat.fit1, type = "pearson"),
    r2 = residuals(dat.fit2, type = "pearson"),
    r3 = residuals(dat.fit3, type = "pearson"),
    S_ID = dat$S_ID
  )
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID)
  
  rsq1 <- rsq.lmm(yi1_on_yi2)
  rsq2 <- rsq.lmm(yi2_on_yi3)
  
  res <- rbind(res, data.frame(
    Item = paired_names[i, 1],
    AIC_1 = AIC(yi1_on_yi2),
    AIC_2 = AIC(yi2_on_yi3),
    AICc_1 = AICc(yi1_on_yi2),
    AICc_2 = AICc(yi2_on_yi3),
    BIC_1 = BIC(yi1_on_yi2),
    BIC_2 = BIC(yi2_on_yi3),
    logLik_1 = logLik(yi1_on_yi2),
    logLik_2 = logLik(yi2_on_yi3),
    R2_model_1 = rsq1$model,
    R2_fixed_1 = rsq1$fixed,
    R2_random_1 = rsq1$random,
    R2_model_2 = rsq2$model,
    R2_fixed_2 = rsq2$fixed,
    R2_random_2 = rsq2$random,
    fix.empty.names = FALSE
  ))
}


Strat <- c("Global", "Herbaceous_plant", "Food_crop", "Cash_crop","Mn__-Generalist","NE_c-Generalist",
           "Intercropping", "Cover_cropping", "Natural", "Temperate", "Tropic")

dat <- read.table(paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", header = TRUE)

for (i in Strat) {
  if (i == "Global") {
    dat_sub <- dat
  } else if (grepl("Mn__-", i)) {
    level <- sub("Mn__-", "", i)
    dat_sub <- subset(dat, Mn__ == level)
  } else if (grepl("NE_c-", i)) {
    level <- sub("NE_c-", "", i)
    dat_sub <- subset(dat, NE_c == level)
  } else {
    dat_sub <- subset(dat, CP_C3 == i | CP_Ct == i | Mn__ == i | NE_c == i | NCP_P == i | Std_ == i | Zone == i)
  }
  
  if (nrow(dat_sub) < 2) next
  
  dat.fit1 <- lme(yi1 ~ 1, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(
    r1 = residuals(dat.fit1, type = "pearson"),
    r2 = residuals(dat.fit2, type = "pearson"),
    r3 = residuals(dat.fit3, type = "pearson"),
    S_ID = dat_sub$S_ID
  )
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID)
  
  rsq1 <- rsq.lmm(yi1_on_yi2)
  rsq2 <- rsq.lmm(yi2_on_yi3)
  
  res <- rbind(res, data.frame(
    Item = i,
    AIC_1 = AIC(yi1_on_yi2),
    AIC_2 = AIC(yi2_on_yi3),
    AICc_1 = AICc(yi1_on_yi2),
    AICc_2 = AICc(yi2_on_yi3),
    BIC_1 = BIC(yi1_on_yi2),
    BIC_2 = BIC(yi2_on_yi3),
    logLik_1 = logLik(yi1_on_yi2),
    logLik_2 = logLik(yi2_on_yi3),
    R2_model_1 = rsq1$model,
    R2_fixed_1 = rsq1$fixed,
    R2_random_1 = rsq1$random,
    R2_model_2 = rsq2$model,
    R2_fixed_2 = rsq2$fixed,
    R2_random_2 = rsq2$random,
    fix.empty.names = FALSE
  ))
}


View(res)
write.csv(res, '/Users/yusha/SuppTab18b-1.csv')





 ####Stable18b-2#### 

rm(list = ls())
library(nlme)
library(MuMIn)
library(rsq)


prefix <- '/Users/yusha/'
paired_names <- read.table(paste0(prefix, "paired_names.txt"), sep = "\t", header = TRUE)


res <- NULL


for (i in 1:nrow(paired_names)) {
  dat <- read.table(paste0(prefix, paired_names[i, 1], ".txt"), sep = "\t", header = TRUE)
  
  dat.fit1 <- lme(yi1 ~ 1 + added_plant_species, dat, random = ~1 | S_ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1 + added_plant_species, dat, random = ~1 | S_ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1 + added_plant_species, dat, random = ~1 | S_ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(
    r1 = residuals(dat.fit1, type = "pearson"),
    r2 = residuals(dat.fit2, type = "pearson"),
    r3 = residuals(dat.fit3, type = "pearson"),
    S_ID = dat$S_ID
  )
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID)
  
  rsq1 <- rsq.lmm(yi1_on_yi2)
  rsq2 <- rsq.lmm(yi2_on_yi3)
  
  res <- rbind(res, data.frame(
    Item = paired_names[i, 1],
    AIC_1 = AIC(yi1_on_yi2),
    AIC_2 = AIC(yi2_on_yi3),
    AICc_1 = AICc(yi1_on_yi2),
    AICc_2 = AICc(yi2_on_yi3),
    BIC_1 = BIC(yi1_on_yi2),
    BIC_2 = BIC(yi2_on_yi3),
    logLik_1 = logLik(yi1_on_yi2),
    logLik_2 = logLik(yi2_on_yi3),
    R2_model_1 = rsq1$model,
    R2_fixed_1 = rsq1$fixed,
    R2_random_1 = rsq1$random,
    R2_model_2 = rsq2$model,
    R2_fixed_2 = rsq2$fixed,
    R2_random_2 = rsq2$random,
    fix.empty.names = FALSE
  ))

}


Strat <- c("Global", "Herbaceous_plant", "Food_crop", "Cash_crop","Mn__-Generalist","NE_c-Generalist",
           "Intercropping", "Cover_cropping", "Natural", "Temperate", "Tropic")

dat <- read.table(paste0(prefix, "enemies_vs_herbivore_vs_crop.txt"), sep = "\t", header = TRUE)

for (i in Strat) {
  if (i == "Global") {
    dat_sub <- dat
  } else if (grepl("Mn__-", i)) {
    level <- sub("Mn__-", "", i)
    dat_sub <- subset(dat, Mn__ == level)
  } else if (grepl("NE_c-", i)) {
    level <- sub("NE_c-", "", i)
    dat_sub <- subset(dat, NE_c == level)
  } else {
    dat_sub <- subset(dat, CP_C3 == i | CP_Ct == i | Mn__ == i | NE_c == i | NCP_P == i | Std_ == i | Zone == i)
  }
  
  if (nrow(dat_sub) < 2) next
  
  dat.fit1 <- lme(yi1 ~ 1 + added_plant_species, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi1), control = lmeControl(sigma = 1))
  dat.fit2 <- lme(yi2 ~ 1 + added_plant_species, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi2), control = lmeControl(sigma = 1))
  dat.fit3 <- lme(yi3 ~ 1 + added_plant_species, dat_sub, random = ~1 | S_ID, weights = varFixed(~vi3), control = lmeControl(sigma = 1))
  
  dat.residual <- data.frame(
    r1 = residuals(dat.fit1, type = "pearson"),
    r2 = residuals(dat.fit2, type = "pearson"),
    r3 = residuals(dat.fit3, type = "pearson"),
    S_ID = dat_sub$S_ID
  )
  
  yi1_on_yi2 <- lme(r2 ~ -1 + r1, dat.residual, random = ~ -1 + r1 | S_ID)
  yi2_on_yi3 <- lme(r3 ~ -1 + r2, dat.residual, random = ~ -1 + r2 | S_ID)
  
  rsq2 <- rsq.lmm(yi2_on_yi3)
  
  res <- rbind(res, data.frame(
    Item = paired_names[i, 1],
    AIC_1 = AIC(yi1_on_yi2),
    AICc_1 = AICc(yi1_on_yi2),
    BIC_1 = BIC(yi1_on_yi2),
    logLik_1 = logLik(yi1_on_yi2),
    R2_model_1 = rsq1$model,
    R2_fixed_1 = rsq1$fixed,
    R2_random_1 = rsq1$random,
    
    AIC_2 = AIC(yi2_on_yi3),
    AICc_2 = AICc(yi2_on_yi3),
    BIC_2 = BIC(yi2_on_yi3),
    logLik_2 = logLik(yi2_on_yi3),
    R2_model_2 = rsq2$model,
    R2_fixed_2 = rsq2$fixed,
    R2_random_2 = rsq2$random,
    
    fix.empty.names = FALSE
  ))
  
}


View(res)

write.csv(res, '/Users/yusha/SuppTab18b-2.csv')








