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




                      ##### Supplementary Table 2.2 ####
Plant_data <- Data[!is.na(Data$yi),]
Plant_data$S_ID <- factor(Plant_data$S_ID)
Plant_data$CP_C3 <- factor(Plant_data$CP_C3)
Plant_data$CP_Ct <- factor(Plant_data$CP_Ct)
Plant_data$Mn__ <- factor(Plant_data$Mn__)
Plant_data$NE_c <- factor(Plant_data$NE_c)
Plant_data$NCP_P <- factor(Plant_data$NCP_P)
Plant_data$Std_ <- factor(Plant_data$Std_)
Plant_data$Zone <- factor(Plant_data$Zone)
Plant_data$Trophic_group <- factor(Plant_data$Trophic_group)
Plant_data$Trophic_group_response <- factor(Plant_data$Trophic_group_response)


Plant_data <- subset(Plant_data,!is.na(Plant_data$CP_Ltn_nm_f_))
Species <- sort(unique(Plant_data$CP_Ltn_nm_f_))
Taxa_Plant <- tnrs_match_names(Species)
Plant_InTree <- is_in_tree(ott_id(Taxa_Plant))
Plant_Tree <- tol_induced_subtree(ott_ids = ott_id(Taxa_Plant)[Plant_InTree],
                                  label_format = "id")
Plant_brlen <- compute.brlen(Plant_Tree)
Plant_Cor <- vcv(Plant_brlen, corr = T)

### make the colnames of Cor become id
for (i in 1:nrow(Plant_Cor)){
  colnames(Plant_Cor)[i] <- substr(colnames(Plant_Cor)[i],4,100)
  rownames(Plant_Cor)[i] <- substr(rownames(Plant_Cor)[i], 4, 100)
}
Colnam_Cor <- colnames(Plant_Cor)
Match_Plant_ootid <- match(Colnam_Cor, Taxa_Plant$ott_id)
Matched_Plant_Sp <- Taxa_Plant[Match_Plant_ootid,"search_string"]
colnames(Plant_Cor) <- Matched_Plant_Sp
rownames(Plant_Cor) <- Matched_Plant_Sp
Species <- stringr::str_to_lower(Species)
Matched_Plant_InCor <- match(Species, colnames(Plant_Cor))
Plant_data <- Plant_data[!is.na(Matched_Plant_InCor),]
Plant_Cor <- Plant_Cor[order(dimnames(Plant_Cor)[[1]]), order(dimnames(Plant_Cor)[[1]])]
Plant_data$Plant_Species <- factor(Plant_data$CP_Ltn_nm_f_)
levels(Plant_data$Plant_Species) <- sort(dimnames(Plant_Cor)[[1]])

Plant_data$Plant_Species <- tolower(as.character(Plant_data$Plant_Species))


common_species <- intersect(Plant_data$Plant_Species, dimnames(Plant_Cor)[[1]])
Plant_data <- subset(Plant_data, Plant_Species %in% common_species)
Plant_Cor <- Plant_Cor[common_species, common_species]


Plant_data$Plant_Species <- factor(Plant_data$Plant_Species, levels = dimnames(Plant_Cor)[[1]])


length(levels(Plant_data$Plant_Species))
length(dimnames(Plant_Cor)[[1]])
setequal(levels(Plant_data$Plant_Species), dimnames(Plant_Cor)[[1]])



Plant_data$Plant_Species2 <- Plant_data$Plant_Species

Pre_vars <- c("CP_Ct","NCP_P","Std_","added_species_richness")#, "CP_Ct", "Mn__", "NE_c", "NCP_P","Std_","Zone","added_species_richness"，"CP_C3"

V <- tcrossprod(sqrt(Plant_data$CONTROL_SD)) /sqrt(outer(Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2,Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2)) +
  diag((Plant_data$TREATMENT_SD)/(Plant_data$TREATMENT_R * (Plant_data$TREATMENT_M)^2))

V_ind <- pmax(as.matrix(crossprod(lFormula(yi~1+(1|S_ID/CONTROL_M), Plant_data)$reTrms$Zt)) - 1, 0)
V <- V*V_ind


library(foreach)
library(doParallel)

cl <- makeCluster(7)        
registerDoParallel(cl)
model.null <- rma.mv(yi, V,mods = ~Trophic_group+SqrtNum-1,
                     random =list(~ 1 | S_ID,
                                  ~ 1 | Plant_Species,
                                  ~ 1 | EsID,
                                  ~ 1 | Plant_Species2),
                     R=list(Plant_Species=Plant_Cor),
                     control = list(optimizer = "optim", optmethod = "BFGS"),
                     data = Plant_data,
                     method = "ML")


egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null))*
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
 


cl <- makeCluster(7)         
registerDoParallel(cl)
res <- data.frame("Trophic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = F)
res <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
               .packages = c("metafor","ape","lme4")) %dopar%{
                 formula.A <- as.formula(paste("~",paste(paste("Trophic_group",Pre_vars[i],sep="+"),"SqrtNum",sep="+"), "-1", sep = ""))
                 formula.B <- as.formula(paste("~", paste(paste("Trophic_group",Pre_vars[i], sep = "+"),"SqrtNum",sep = "+"),
                                               paste("+Trophic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
                 TemData <- Plant_data
                 
                 if (Pre_vars[i]=="Clmz"){
                   TemData.Clmz <- subset(TemData,Etioge=="Outdoor")
                   V <- tcrossprod(sqrt(Plant_data$CONTROL_SD)) /sqrt(outer(Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2,Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2)) +
                     diag((Plant_data$TREATMENT_SD)/(Plant_data$TREATMENT_R * (Plant_data$TREATMENT_M)^2))
                   V_ind <- pmax(as.matrix(crossprod(lFormula(yi~1+(1|S_ID/CONTROL_M), Plant_data)$reTrms$Zt)) - 1, 0)
                   V.Clmz <- V.Clmz*V_ind
                   model.A <- rma.mv(yi, V, mods = formula.A,
                                     random =list(~ 1 | S_ID,
                                                  ~ 1 | Plant_Species,
                                                  ~ 1 | EsID,
                                                  ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData.Clmz,
                                     method = "ML")
                   model.B <- rma.mv(yi, V, mods = formula.B,
                                     random =list( ~ 1 | S_ID,
                                                   ~ 1 | Plant_Species,
                                                   ~ 1 | EsID,
                                                   ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData.Clmz,
                                     method = "ML")
                   egger_dt <- data.frame(r = residuals(model.A),
                                          s = sign(residuals(model.A))*
                                            sqrt(diag(vcov(model.A, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.A)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   egger_dt <- data.frame(r = residuals(model.B),
                                          s = sign(residuals(model.B))*
                                            sqrt(diag(vcov(model.B, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.B)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   
                 }else{
                   model.A <- rma.mv(yi, V,mods = formula.A,
                                     random =list(  ~ 1 | S_ID,
                                                    ~ 1 | Plant_Species,
                                                    ~ 1 | EsID,
                                                    ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData,
                                     method = "ML")
                   model.B <- rma.mv(yi, V, mods = formula.B,
                                     random =list( ~ 1 | S_ID,
                                                   ~ 1 | Plant_Species,
                                                   ~ 1 | EsID,
                                                   ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData,
                                     method = "ML")
                   egger_dt <- data.frame(r = residuals(model.A),
                                          s = sign(residuals(model.A))*
                                            sqrt(diag(vcov(model.A, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.A)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   egger_dt <- data.frame(r = residuals(model.B),
                                          s = sign(residuals(model.B))*
                                            sqrt(diag(vcov(model.B, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.B)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                 }
               }

# V <- vcalc(vi,cluster = S_ID, obs = EsID, data = Plant_data, rho = 0.5)
       
model.null2 <- rma.mv(yi, V,mods = ~Trophic_group_response+SqrtNum-1,
                      random =list( ~ 1 | S_ID,
                                    ~ 1 | Plant_Species,
                                    ~ 1 | EsID,
                                    ~ 1 | Plant_Species2),
                      R=list(Plant_Species=Plant_Cor),
                      control = list(optimizer = "optim", optmethod = "BFGS"),
                      data = Plant_data,
                      method = "ML")
egger_dt <- data.frame(r = residuals(model.null2),
                       s = sign(residuals(model.null2))*
                         sqrt(diag(vcov(model.null2, type = "resid"))))
egg_test <- lm(r~s, egger_dt)
res <- rbind(res,data.frame("Trophic Group Response",
                            nrow(egg_test$model),
                            coef(summary(egg_test))[1, 3],
                            coef(summary(egg_test))[1, 4],
                            fix.empty.names = F))
res <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
               .packages = c("metafor","ape","lme4")) %dopar%{
                 formula.A <- as.formula(paste("~",paste(paste("Trophic_group_response",Pre_vars[i],sep="+"),"SqrtNum",sep="+"), "-1", sep = ""))
                 formula.B <- as.formula(paste("~", paste(paste("Trophic_group_response",Pre_vars[i], sep = "+"),"SqrtNum",sep='+'),
                                               paste("+Trophic_group_response",Pre_vars[i],sep = ":"), "-1", sep = ""))
                 TemData <- Plant_data
                 
                 if (Pre_vars[i]=="Clmz"){
                   TemData.Clmz <- subset(TemData,Etioge=="Outdoor")
                   V <- tcrossprod(sqrt(Plant_data$CONTROL_SD)) /sqrt(outer(Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2,Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2)) +
                     diag((Plant_data$TREATMENT_SD)/(Plant_data$TREATMENT_R * (Plant_data$TREATMENT_M)^2))
                   V_ind <- pmax(as.matrix(crossprod(lFormula(yi~1+(1|S_ID/CONTROL_M), Plant_data)$reTrms$Zt)) - 1, 0)
                   V.Clmz <- V.Clmz*V_ind
                   model.A <- rma.mv(yi, V, mods = formula.A,
                                     random =list( ~ 1 | S_ID,
                                                   ~ 1 | Plant_Species,
                                                   ~ 1 | EsID,
                                                   ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData.Clmz,
                                     method = "ML")
                   model.B <- rma.mv(yi, V, mods = formula.B,
                                     random =list( ~ 1 | S_ID,
                                                   ~ 1 | Plant_Species,
                                                   ~ 1 | EsID,
                                                   ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData.Clmz,
                                     method = "ML")
                   egger_dt <- data.frame(r = residuals(model.A),
                                          s = sign(residuals(model.A))*
                                            sqrt(diag(vcov(model.A, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.A)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   egger_dt <- data.frame(r = residuals(model.B),
                                          s = sign(residuals(model.B))*
                                            sqrt(diag(vcov(model.B, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.B)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   
                 }else{
                   model.A <- rma.mv(yi, V,mods = formula.A,
                                     random =list(  ~ 1 | S_ID,
                                                    ~ 1 | Plant_Species,
                                                    ~ 1 | EsID,
                                                    ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData,
                                     method = "ML")
                   model.B <- rma.mv(yi, V, mods = formula.B,
                                     random =list( ~ 1 | S_ID,
                                                   ~ 1 | Plant_Species,
                                                   ~ 1 | EsID,
                                                   ~ 1 | Plant_Species2),
                                     R=list(Plant_Species=Plant_Cor),
                                     control = list(optimizer = "optim", optmethod = "BFGS"),
                                     data = TemData,
                                     method = "ML")
                   egger_dt <- data.frame(r = residuals(model.A),
                                          s = sign(residuals(model.A))*
                                            sqrt(diag(vcov(model.A, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.A)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                   egger_dt <- data.frame(r = residuals(model.B),
                                          s = sign(residuals(model.B))*
                                            sqrt(diag(vcov(model.B, type = "resid"))))
                   egg_test <- lm(r~s, egger_dt)
                   res <- rbind(res, data.frame(as.character(formula.B)[-1],
                                                nrow(egg_test$model),
                                                coef(summary(egg_test))[1, 3],
                                                coef(summary(egg_test))[1, 4],
                                                fix.empty.names = F))
                 }
               }
colnames(res) <- c("Predictor variables",
                   "N",
                   "test-value",
                   "P-value")
res <- res[!duplicated(res[,1]),]
stopCluster(cl)
View(res)
write.csv(res, paste('SupplementTable2_2.2.csv'), row.names = F)

                      
##### Supplementary Table 2.2 with NA ####


Plant_data <- Data[!is.na(Data$yi), ]
Plant_data$S_ID <- factor(Plant_data$S_ID)
Plant_data$CP_C3 <- factor(Plant_data$CP_C3)
Plant_data$CP_Ct <- factor(Plant_data$CP_Ct)
Plant_data$Mn__ <- factor(Plant_data$Mn__)
Plant_data$NE_c <- factor(Plant_data$NE_c)
Plant_data$NCP_P <- factor(Plant_data$NCP_P)
Plant_data$Std_ <- factor(Plant_data$Std_)
Plant_data$Zone <- factor(Plant_data$Zone)
Plant_data$Trophic_group <- factor(Plant_data$Trophic_group)
Plant_data$Trophic_group_response <- factor(Plant_data$Trophic_group_response)

Plant_data <- subset(Plant_data, !is.na(Plant_data$CP_Ltn_nm_f_))
Species <- sort(unique(Plant_data$CP_Ltn_nm_f_))
Taxa_Plant <- tnrs_match_names(Species)
Plant_InTree <- is_in_tree(ott_id(Taxa_Plant))
Plant_Tree <- tol_induced_subtree(ott_ids = ott_id(Taxa_Plant)[Plant_InTree],
                                  label_format = "id")
Plant_brlen <- compute.brlen(Plant_Tree)
Plant_Cor <- vcv(Plant_brlen, corr = TRUE)

for (i in 1:nrow(Plant_Cor)) {
  colnames(Plant_Cor)[i] <- substr(colnames(Plant_Cor)[i], 4, 100)
  rownames(Plant_Cor)[i] <- substr(rownames(Plant_Cor)[i], 4, 100)
}

Colnam_Cor <- colnames(Plant_Cor)
Match_Plant_ootid <- match(Colnam_Cor, Taxa_Plant$ott_id)
Matched_Plant_Sp <- Taxa_Plant[Match_Plant_ootid, "search_string"]
colnames(Plant_Cor) <- Matched_Plant_Sp
rownames(Plant_Cor) <- Matched_Plant_Sp

Species <- stringr::str_to_lower(Species)
Matched_Plant_InCor <- match(Species, colnames(Plant_Cor))
Plant_data <- Plant_data[!is.na(Matched_Plant_InCor), ]

Plant_Cor <- Plant_Cor[order(dimnames(Plant_Cor)[[1]]),
                       order(dimnames(Plant_Cor)[[1]])]
Plant_data$Plant_Species <- tolower(as.character(Plant_data$CP_Ltn_nm_f_))
common_species <- intersect(Plant_data$Plant_Species, dimnames(Plant_Cor)[[1]])
Plant_data <- subset(Plant_data, Plant_Species %in% common_species)
Plant_Cor <- Plant_Cor[common_species, common_species]
Plant_data$Plant_Species <- factor(Plant_data$Plant_Species,
                                   levels = dimnames(Plant_Cor)[[1]])
Plant_data$Plant_Species2 <- Plant_data$Plant_Species

### ==== VCV ====
V <- tcrossprod(sqrt(Plant_data$CONTROL_SD)) /
  sqrt(outer(Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2,
             Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2)) +
  diag((Plant_data$TREATMENT_SD) /
         (Plant_data$TREATMENT_R * (Plant_data$TREATMENT_M)^2))


V_ind <- tryCatch({
  if (nrow(Plant_data) > 1 &&
      length(unique(Plant_data$S_ID)) > 1 &&
      length(unique(Plant_data$CONTROL_M)) > 1) {
    pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID/CONTROL_M),
                                      Plant_data)$reTrms$Zt)) - 1, 0)
  } else {
    message("⚠️ not enough sample sizes")
    diag(nrow(Plant_data))
  }
}, error = function(e) {
  message("⚠️ VCV failed：", e$message)
  diag(nrow(Plant_data))
})
V <- V * V_ind

### ==== 基础模型 ====
library(metafor)
library(ape)
library(lme4)
library(foreach)
library(doParallel)

cl <- makeCluster(6)
registerDoParallel(cl)

model.null <- rma.mv(yi, V,
                     mods = ~ Trophic_group + SqrtNum - 1,
                     random = list(~1 | S_ID,
                                   ~1 | Plant_Species,
                                   ~1 | EsID,
                                   ~1 | Plant_Species2),
                     R = list(Plant_Species = Plant_Cor),
                     control = list(optimizer = "optim",
                                    optmethod = "BFGS"),
                     data = Plant_data, method = "ML")

egger_dt <- data.frame(r = residuals(model.null),
                       s = sign(residuals(model.null)) *
                         sqrt(diag(vcov(model.null, type = "resid"))))
egg_test <- lm(r ~ s, egger_dt)

res <- data.frame("Trophic Group",
                  nrow(egg_test$model),
                  coef(summary(egg_test))[1, 3],
                  coef(summary(egg_test))[1, 4],
                  fix.empty.names = FALSE)


cl <- makeCluster(6)  
registerDoParallel(cl)

Pre_vars <- c("Mn__", "NE_c", "Zone")
Group_vars <- c("Trophic_group", "Trophic_group_response") 

res <- foreach(g = Group_vars, .combine = "rbind",
               .packages = c("metafor", "ape", "lme4")) %:%
  foreach(i = 1:length(Pre_vars), .combine = "rbind") %dopar% {
    
    this_var <- Pre_vars[i]
    group_var <- g
    
    formula.A <- as.formula(paste("~", paste(paste(group_var, this_var, sep = "+"),
                                             "SqrtNum", sep = "+"), "-1"))
    formula.B <- as.formula(
      paste("~", paste(group_var, this_var, sep = "+"),
            "+ SqrtNum",
            "+", paste(group_var, this_var, sep = ":"),
            "-1")
    )
    
    # delete NA
    TemData <- Plant_data[!is.na(Plant_data[[this_var]]), ]
    if (nrow(TemData) <= 3) {
      message("⚠️ ", group_var, "-", this_var, " not enough data，skip")
      return(data.frame(Group = group_var,
                        Predictor = this_var,
                        N = nrow(TemData),
                        tval_A = NA, pval_A = NA,
                        tval_B = NA, pval_B = NA))
    }
    
    # VCV construction
    V_tmp <- tryCatch({
      if (nrow(TemData) > 1 &&
          length(unique(TemData$S_ID)) > 1 &&
          length(unique(TemData$CONTROL_M)) > 1) {
        Vbase <- tcrossprod(sqrt(TemData$CONTROL_SD)) /
          sqrt(outer(TemData$CONTROL_R * (TemData$CONTROL_M)^2,
                     TemData$CONTROL_R * (TemData$CONTROL_M)^2)) +
          diag((TemData$TREATMENT_SD) /
                 (TemData$TREATMENT_R * (TemData$TREATMENT_M)^2))
        Vind <- pmax(as.matrix(crossprod(
          lFormula(yi ~ 1 + (1 | S_ID/CONTROL_M), TemData)$reTrms$Zt)) - 1, 0)
        Vbase * Vind
      } else diag(nrow(TemData))
    }, error = function(e) {
      message("⚠️ ", group_var, "-", this_var, " VCV falied：", e$message)
      diag(nrow(TemData))
    })
    
    # model
    model.A <- rma.mv(yi, V_tmp, mods = formula.A,
                      random = list(~1 | S_ID,
                                    ~1 | Plant_Species,
                                    ~1 | EsID,
                                    ~1 | Plant_Species2),
                      R = list(Plant_Species = Plant_Cor),
                      control = list(optimizer = "optim", optmethod = "BFGS"),
                      data = TemData, method = "ML")
    
    model.B <- rma.mv(yi, V_tmp, mods = formula.B,
                      random = list(~1 | S_ID,
                                    ~1 | Plant_Species,
                                    ~1 | EsID,
                                    ~1 | Plant_Species2),
                      R = list(Plant_Species = Plant_Cor),
                      control = list(optimizer = "optim", optmethod = "BFGS"),
                      data = TemData, method = "ML")
    
    # Egger test
    egger_dt.A <- data.frame(r = residuals(model.A),
                             s = sign(residuals(model.A)) *
                               sqrt(diag(vcov(model.A, type = "resid"))))
    egg_test.A <- lm(r ~ s, egger_dt.A)
    
    egger_dt.B <- data.frame(r = residuals(model.B),
                             s = sign(residuals(model.B)) *
                               sqrt(diag(vcov(model.B, type = "resid"))))
    egg_test.B <- lm(r ~ s, egger_dt.B)
    
    # Results
    data.frame(Group = group_var,
               Predictor = this_var,
               N = nrow(TemData),
               tval_A = coef(summary(egg_test.A))[1, 3],
               pval_A = coef(summary(egg_test.A))[1, 4],
               tval_B = coef(summary(egg_test.B))[1, 3],
               pval_B = coef(summary(egg_test.B))[1, 4])
  }

stopCluster(cl)

View(res)
write.csv(res, paste('SupplementTable2_2.2-with NA.csv'), row.names = F)

