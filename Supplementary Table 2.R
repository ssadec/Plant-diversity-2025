##### Supplementary Table 2.1####
library(readxl)
library(metafor)
library(dplyr)
library(ape)
library(rotl)
library(doParallel)
library(foreach)
library(stringr)
library(parallel)
library(lme4)


Data <- read.table("Alldata.txt", sep = "\t", header = TRUE)
Whole_data <- Data[!is.na(Data$yi), ]
Whole_data$S_ID <- factor(Whole_data$S_ID)
Whole_data$Trophic_group <- factor(Whole_data$Trophic_group)
Whole_data$Trophic_group_response <- factor(Whole_data$Trophic_group_response)
Whole_data$SqrtNum <- sqrt((Whole_data$CONTROL_R + Whole_data$TREATMENT_R) /
                             (Whole_data$CONTROL_R * Whole_data$TREATMENT_R))


Whole_data <- Whole_data[order(Whole_data$S_ID), ]
Cod <- table(Whole_data$S_ID)
esid <- unlist(lapply(Cod, function(x) seq(1, x, 1)))
Whole_data$EsID <- esid


V <- tcrossprod(sqrt(Whole_data$CONTROL_SD)) /
  sqrt(outer(Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2,
             Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2)) +
  diag((Whole_data$TREATMENT_SD) / (Whole_data$TREATMENT_R * (Whole_data$TREATMENT_M)^2))

V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), Whole_data)$reTrms$Zt)) - 1, 0)
V <- V * V_ind

####part1####
# -------------------- Trophic_group--------------------
model.null <- rma.mv(yi, V, mods = ~Trophic_group - 1,
                   random = list(~1 | S_ID, ~1 | EsID),
                   data = Whole_data, method = "ML")

egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null)) * sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r ~ s, data = egger_dt)

res_null<- data.frame(
  Predictor = "Trophic_group",
  N = nrow(egger_dt),
  test_value = coef(summary(egg_test))[1, 3],
  P_value = coef(summary(egg_test))[1, 4]
)
View(res_null)

# -------------------- recycle Pre_vars without NA category--------------------
Pre_vars <- c("CP_Ct","CP_C3", "NCP_P" ,"Std_", "added_plant_species")

num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res_loop <- foreach(varname = Pre_vars, .combine = "rbind",
                    .packages = c("metafor","lme4")) %dopar% {
                      
                      Tem_data <- Whole_data[!is.na(Whole_data[[varname]]), ]
                      
                      formula.A <- as.formula(paste("~Trophic_group+", varname, "+SqrtNum -1"))
                      formula.B <- as.formula(paste("~Trophic_group+", varname, "+SqrtNum+", 
                                                    "Trophic_group:", varname, "-1"))
                      
                      # model A
                      model.A <- rma.mv(yi, V, mods = formula.A,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_A <- data.frame(r = residuals(model.A),
                                               s = sign(residuals(model.A)) * sqrt(diag(vcov(model.A, type="resid"))))
                      egg_test_A <- lm(r ~ s, data = egger_dt_A)
                      
                      res_A <- data.frame(
                        Predictor = varname,
                        Model = "A",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_A))[1,3],
                        P_value = coef(summary(egg_test_A))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                      # model B）
                      model.B <- rma.mv(yi, V, mods = formula.B,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_B <- data.frame(r = residuals(model.B),
                                               s = sign(residuals(model.B)) * sqrt(diag(vcov(model.B, type="resid"))))
                      egg_test_B <- lm(r ~ s, data = egger_dt_B)
                      
                      res_B <- data.frame(
                        Predictor = varname,
                        Model = "B_interaction",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_B))[1,3],
                        P_value = coef(summary(egg_test_B))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                     
                      rbind(res_A, res_B)
                    }

stopCluster(cl)
View(res_loop)
write.csv(res, paste(path,"SupplementTable2.1-1.csv",sep = ""), row.names = F)




# -------------------- recycle Pre_vars with NA category--------------------
##### Supplementary Table 2.1-part 2 ####
library(readxl)
library(metafor)
library(dplyr)
library(ape)
library(rotl)
library(doParallel)
library(foreach)
library(stringr)
library(parallel)
library(lme4)


Data <- read.table("Alldata.txt", sep = "\t", header = TRUE)
Whole_data <- Data[!is.na(Data$yi), ]
Whole_data$S_ID <- factor(Whole_data$S_ID)
Whole_data$Trophic_group <- factor(Whole_data$Trophic_group)
Whole_data$Trophic_group_response <- factor(Whole_data$Trophic_group_response)
Whole_data$SqrtNum <- sqrt((Whole_data$CONTROL_R + Whole_data$TREATMENT_R) /
                             (Whole_data$CONTROL_R * Whole_data$TREATMENT_R))


Whole_data <- Whole_data[order(Whole_data$S_ID), ]
Cod <- table(Whole_data$S_ID)
esid <- unlist(lapply(Cod, function(x) seq(1, x, 1)))
Whole_data$EsID <- esid


V_full <- tcrossprod(sqrt(Whole_data$CONTROL_SD)) /
  sqrt(outer(Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2,
             Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2)) +
  diag((Whole_data$TREATMENT_SD) / (Whole_data$TREATMENT_R * (Whole_data$TREATMENT_M)^2))

V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), Whole_data)$reTrms$Zt)) - 1, 0)
V_full <- V_full * V_ind


Pre_vars <- c("Zone")#"Mn__","NE_c"

num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res_loop <- foreach(varname = Pre_vars, .combine = "rbind",
                    .packages = c("metafor","lme4")) %dopar% {
                      
                 
                      Tem_data <- Whole_data[!is.na(Whole_data[[varname]]), ]
                      idx <- as.numeric(rownames(Tem_data)) 
                      V <- V_full[idx, idx] 
                      
                      formula.A <- as.formula(paste("~Trophic_group+", varname, "+SqrtNum -1"))
                      formula.B <- as.formula(paste("~Trophic_group+", varname, "+SqrtNum+", 
                                                    "Trophic_group:", varname, "-1"))
                      
                     
                      model.A <- rma.mv(yi, V, mods = formula.A,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_A <- data.frame(r = residuals(model.A),
                                               s = sign(residuals(model.A)) * sqrt(diag(vcov(model.A, type="resid"))))
                      egg_test_A <- lm(r ~ s, data = egger_dt_A)
                      
                      res_A <- data.frame(
                        Predictor = varname,
                        Model = "A",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_A))[1,3],
                        P_value = coef(summary(egg_test_A))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                     
                      model.B <- rma.mv(yi, V, mods = formula.B,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_B <- data.frame(r = residuals(model.B),
                                               s = sign(residuals(model.B)) * sqrt(diag(vcov(model.B, type="resid"))))
                      egg_test_B <- lm(r ~ s, data = egger_dt_B)
                      
                      res_B <- data.frame(
                        Predictor = varname,
                        Model = "B_interaction",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_B))[1,3],
                        P_value = coef(summary(egg_test_B))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                  
                      rbind(res_A, res_B)
                    }

stopCluster(cl)

View(res_loop)


write.csv(res_loop, "SupplementTable2.1-1.csv", row.names = FALSE)




####part2####
# -------------------- Trophic_group_response  --------------------
model.null <- rma.mv(yi, V, mods = ~Trophic_group_response - 1,
                   random = list(~1 | S_ID, ~1 | EsID),
                   data = Whole_data, method = "ML")

egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null)) * sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r ~ s, data = egger_dt)

res_TG <- data.frame(
  Predictor = "Trophic_group_response",
  N = nrow(egger_dt),
  test_value = coef(summary(egg_test))[1, 3],
  P_value = coef(summary(egg_test))[1, 4]
)
View(res_TG)

# -------------------- recycle Pre_vars without NA category--------------------
Pre_vars <- c("CP_Ct","CP_C3", "NCP_P" ,"Std_", "added_plant_species")

num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res_loop <- foreach(varname = Pre_vars, .combine = "rbind",
                    .packages = c("metafor","lme4")) %dopar% {
                      
                      Tem_data <- Whole_data[!is.na(Whole_data[[varname]]), ]
                      
                      formula.A <- as.formula(paste("~Trophic_group_response+", varname, "+SqrtNum -1"))
                      formula.B <- as.formula(paste("~Trophic_group_response+", varname, "+SqrtNum+", 
                                                    "Trophic_group_response:", varname, "-1"))
                      
                      
                      model.A <- rma.mv(yi, V, mods = formula.A,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_A <- data.frame(r = residuals(model.A),
                                               s = sign(residuals(model.A)) * sqrt(diag(vcov(model.A, type="resid"))))
                      egg_test_A <- lm(r ~ s, data = egger_dt_A)
                      
                      res_A <- data.frame(
                        Predictor = varname,
                        Model = "A",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_A))[1,3],
                        P_value = coef(summary(egg_test_A))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                     
                      model.B <- rma.mv(yi, V, mods = formula.B,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_B <- data.frame(r = residuals(model.B),
                                               s = sign(residuals(model.B)) * sqrt(diag(vcov(model.B, type="resid"))))
                      egg_test_B <- lm(r ~ s, data = egger_dt_B)
                      
                      res_B <- data.frame(
                        Predictor = varname,
                        Model = "B_interaction",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_B))[1,3],
                        P_value = coef(summary(egg_test_B))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                      
                      rbind(res_A, res_B)
                    }

stopCluster(cl)
View(res_loop)
write.csv(res, paste(path,"SupplementTable2.1-1.csv",sep = ""), row.names = F)




# -------------------- recycle Pre_vars with NA category--------------------
##### Supplementary Table 2.1-part 2 ####
library(readxl)
library(metafor)
library(dplyr)
library(ape)
library(rotl)
library(doParallel)
library(foreach)
library(stringr)
library(parallel)
library(lme4)


Data <- read.table("Alldata.txt", sep = "\t", header = TRUE)
Whole_data <- Data[!is.na(Data$yi), ]
Whole_data$S_ID <- factor(Whole_data$S_ID)
Whole_data$Trophic_group <- factor(Whole_data$Trophic_group)
Whole_data$Trophic_group_response <- factor(Whole_data$Trophic_group_response)
Whole_data$SqrtNum <- sqrt((Whole_data$CONTROL_R + Whole_data$TREATMENT_R) /
                             (Whole_data$CONTROL_R * Whole_data$TREATMENT_R))


Whole_data <- Whole_data[order(Whole_data$S_ID), ]
Cod <- table(Whole_data$S_ID)
esid <- unlist(lapply(Cod, function(x) seq(1, x, 1)))
Whole_data$EsID <- esid


V_full <- tcrossprod(sqrt(Whole_data$CONTROL_SD)) /
  sqrt(outer(Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2,
             Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2)) +
  diag((Whole_data$TREATMENT_SD) / (Whole_data$TREATMENT_R * (Whole_data$TREATMENT_M)^2))

V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), Whole_data)$reTrms$Zt)) - 1, 0)
V_full <- V_full * V_ind


Pre_vars <- c("NE_c")#"Mn__","Zone"

num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res_loop <- foreach(varname = Pre_vars, .combine = "rbind",
                    .packages = c("metafor","lme4")) %dopar% {
                      
                    
                      Tem_data <- Whole_data[!is.na(Whole_data[[varname]]), ]
                      idx <- as.numeric(rownames(Tem_data))
                      V <- V_full[idx, idx]
                      
                      formula.A <- as.formula(paste("~Trophic_group_response+", varname, "+SqrtNum -1"))
                      formula.B <- as.formula(paste("~Trophic_group_response+", varname, "+SqrtNum+", 
                                                    "Trophic_group_response:", varname, "-1"))
                      
                  
                      model.A <- rma.mv(yi, V, mods = formula.A,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_A <- data.frame(r = residuals(model.A),
                                               s = sign(residuals(model.A)) * sqrt(diag(vcov(model.A, type="resid"))))
                      egg_test_A <- lm(r ~ s, data = egger_dt_A)
                      
                      res_A <- data.frame(
                        Predictor = varname,
                        Model = "A",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_A))[1,3],
                        P_value = coef(summary(egg_test_A))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                  
                      model.B <- rma.mv(yi, V, mods = formula.B,
                                        random = list(~1|S_ID, ~1|EsID),
                                        data = Tem_data, method = "ML")
                      
                      egger_dt_B <- data.frame(r = residuals(model.B),
                                               s = sign(residuals(model.B)) * sqrt(diag(vcov(model.B, type="resid"))))
                      egg_test_B <- lm(r ~ s, data = egger_dt_B)
                      
                      res_B <- data.frame(
                        Predictor = varname,
                        Model = "B_interaction",
                        N = nrow(Tem_data),
                        test_value = coef(summary(egg_test_B))[1,3],
                        P_value = coef(summary(egg_test_B))[1,4],
                        stringsAsFactors = FALSE
                      )
                      
                   
                      rbind(res_A, res_B)
                    }

stopCluster(cl)

View(res_loop)
write.csv(res_loop, "SupplementTable2.1-1.csv", row.names = FALSE)
