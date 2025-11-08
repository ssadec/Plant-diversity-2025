rm(list = ls())
library(metafor)
library(rotl)
library(ape)
library(dplyr)
library(stringr)
library(lme4)
library(orchaRd)
library(parallel)
library(doParallel)
library(lqmm)


prefix <- "/Users/yusha/"
Data <- read.table("Alldata.txt", sep = "\t", header = TRUE)

####Supplementary Table 1.1####
Data <- Data[order(Data$S_ID), ]
Cod <- table(Data$S_ID)
esid <- unlist(lapply(Cod, function(x) seq(1, x, 1)))
Data$EsID <- esid


Whole_data <- Data[!is.na(Data$yi), ]
Whole_data$S_ID <- factor(Whole_data$S_ID)
Whole_data$CP_C3 <- factor(Whole_data$CP_C3)
Whole_data$CP_Ct <- factor(Whole_data$CP_Ct)
Whole_data$Mn__ <- factor(Whole_data$Mn__)
Whole_data$NE_c <- factor(Whole_data$NE_c)
Whole_data$NCP_P <- factor(Whole_data$NCP_P)
Whole_data$Std_ <- factor(Whole_data$Std_)
Whole_data$Zone <- factor(Whole_data$Zone)
Whole_data$Trophic_group <- factor(Whole_data$Trophic_group)
Whole_data$Trophic_group_response <- factor(Whole_data$Trophic_group_response)


###part1####
V <- tcrossprod(sqrt(Whole_data$CONTROL_SD))/sqrt(outer(Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2, Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2))+
  diag((Whole_data$TREATMENT_SD)/(Whole_data$TREATMENT_R * (Whole_data$TREATMENT_M)^2))
V_ind <- pmax(as.matrix(crossprod(lFormula(yi~1+(1|S_ID/CONTROL_M), Whole_data)$reTrms$Zt)) - 1, 0)
V <- V*V_ind

model.null <- rma.mv(yi, V,
                     mods = ~Trophic_group-1, 
                     random = list(~ 1|S_ID,~1|EsID), 
                     data = Whole_data, control = list(optimize = "BFGS"),
                     method = "ML")
I2 <- i2_ml(model.null)
Results <- data.frame("Trophic_group",
                      NA,
                      AIC(model.null),
                      logLik(model.null),
                      NA,
                      summary(model.null)$p,
                      NA,
                      summary(model.null)$k,
                      I2,
                      fix.empty.names = F)
Pre_vars <- c("CP_C3","CP_Ct","NCP_P","Std_","added_plant_species")


num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~",paste("Trophic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group",Pre_vars[i], sep = "+"),
                                                   paste("+Trophic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
                     
                     Tem_data <- Whole_data
                     Tem_data <- Tem_data[!is.na(Tem_data[[Pre_vars[i]]]), ]  
                     
                     model.A <- rma.mv(yi, V, mods = formula.A,
                                       random = list(~1|S_ID,~1|EsID),
                                       data = Tem_data,control = list(optimize = "BFGS"),
                                       method = "ML")
                     
                     model.B <- rma.mv(yi, V, mods = formula.B,
                                       random = list(~1|S_ID,~1|EsID),control = list(optimize = "BFGS"),
                                       data = Tem_data,
                                       method = "ML")
                     
                     I2.A <- i2_ml(model.A)
                     I2.B <- i2_ml(model.B)
                     
                     # Compare model B vs A
                     if (length(model.A$beta) != length(model.B$beta)){
                       BA <- anova(model.B, model.A)
                       BA_res <- data.frame(
                         formula = as.character(formula.B)[-1],
                         Model = paste("1.", i, "A", sep = ""),
                         AIC = BA[["fit.stats.f"]][["AIC"]],
                         logLik = BA[["fit.stats.f"]][["ll"]],
                         LRT = BA[["LRT"]],
                         df = BA[["parms.f"]],
                         pval = BA[["pval"]],
                         k = summary(model.B)$k,
                         I2 = I2.B,
                         stringsAsFactors = FALSE
                       )
                     } else {
                       BA_res <- data.frame(
                         formula = as.character(formula.B)[-1],
                         Model = paste("1.", i, "A", sep = ""),
                         AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA,
                         k = summary(model.B)$k,
                         I2 = I2.B,
                         stringsAsFactors = FALSE
                       )
                     }
                     
                     # Also record model A info
                     A1 <- anova(model.A, model.null)
                     A_res <- data.frame(
                       formula = as.character(formula.A)[-1],
                       Model = "1",
                       AIC = A1[["fit.stats.f"]][["AIC"]],
                       logLik = A1[["fit.stats.f"]][["ll"]],
                       LRT = A1[["LRT"]],
                       df = A1[["parms.f"]],
                       pval = A1[["pval"]],
                       k = summary(model.A)$k,
                       I2 = I2.A,
                       stringsAsFactors = FALSE
                     )
                     
                     # 合并 A 模型和 B 模型的结果
                     rbind(A_res, BA_res)
                   }

View(Results)

write.csv(Results,'/Users/yusha/Stable1a.csv')




###part2####

library(metafor)
library(ape)
library(orchaRd)
library(lme4)
library(doParallel)

Pre_vars <- c("Mn__") #,"Zone","Mn__","NE_c"
num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


colnames_final <- c("Model", "Step", "AIC", "logLik", "LRT", "df", "pval", "k", "I2")

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~ Trophic_group +", Pre_vars[i], "-1"))
                     
                     formula.B <- as.formula(paste("~ Trophic_group +", Pre_vars[i],
                                                   "+ Trophic_group:", Pre_vars[i], "-1"))
                     
                     Tem_data <- Whole_data
                     
                   
                     if (Pre_vars[i] == "Mn__") {  #"Zone","NE_c"
                       
                  
                       Tem_data.Clmz <- Tem_data[!is.na(Tem_data$Mn__) & !is.na(Tem_data$Trophic_group), ]  #"Zone","NE_c"
                       
                      
                       cat(">>> subset distribution：\n")
                       print(table(Tem_data.Clmz$Trophic_group, Tem_data.Clmz$Mn__))  #"Zone","NE_c"
                       
                       
                       V.Clmz <- tcrossprod(sqrt(Tem_data.Clmz$CONTROL_SD)) / 
                         sqrt(outer(Tem_data.Clmz$CONTROL_R * (Tem_data.Clmz$CONTROL_M)^2,
                                    Tem_data.Clmz$CONTROL_R * (Tem_data.Clmz$CONTROL_M)^2)) +
                         diag((Tem_data.Clmz$TREATMENT_SD) / (Tem_data.Clmz$TREATMENT_R * (Tem_data.Clmz$TREATMENT_M)^2))
                       
                       V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1|S_ID/CONTROL_M), Tem_data.Clmz)$reTrms$Zt)) - 1, 0)
                       V.Clmz <- V.Clmz * V_ind
                       
                      
                       model.A <- rma.mv(yi, V.Clmz, mods = formula.A,
                                         random = list(~1|S_ID, ~1|EsID),
                                         data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                         method = "ML")
                       
                 
                       model.B <- rma.mv(yi, V.Clmz, mods = formula.B,
                                         random = list(~1|S_ID, ~1|EsID),
                                         data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                         method = "ML")
                       
                      
                       model.Clmz <- rma.mv(yi, V.Clmz, mods = ~Trophic_group - 1,
                                            random = list(~1|S_ID, ~1|EsID),
                                            data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                            method = "ML")
                       
                   
                       A1 <- anova(model.A, model.Clmz)
                       I2.A <- i2_ml(model.A)
                       I2.B <- i2_ml(model.B)
                       
                   
                       if (length(model.A$beta) != length(model.B$beta)) {
                         BA <- anova(model.B, model.A)
                         BA_res <- data.frame("Model" = as.character(formula.B)[-1],
                                              "Step" = paste("1.", i, "A", sep = ""),
                                              "AIC" = BA[["fit.stats.f"]][["AIC"]],
                                              "logLik" = BA[["fit.stats.f"]][["ll"]],
                                              "LRT" = BA[["LRT"]],
                                              "df" = BA[["parms.f"]],
                                              "pval" = BA[["pval"]],
                                              "k" = summary(model.B)$k,
                                              "I2" = I2.B,
                                              stringsAsFactors = FALSE)
                       } else {
                         BA_res <- data.frame("Model" = as.character(formula.B)[-1],
                                              "Step" = paste("1.", i, "A", sep = ""),
                                              "AIC" = NA,
                                              "logLik" = NA,
                                              "LRT" = NA,
                                              "df" = NA,
                                              "pval" = NA,
                                              "k" = summary(model.B)$k,
                                              "I2" = I2.B,
                                              stringsAsFactors = FALSE)
                       }
                       
                       A_res <- data.frame("Model" = as.character(formula.A)[-1],
                                           "Step" = "1",
                                           "AIC" = A1[["fit.stats.f"]][["AIC"]],
                                           "logLik" = A1[["fit.stats.f"]][["ll"]],
                                           "LRT" = A1[["LRT"]],
                                           "df" = A1[["parms.f"]],
                                           "pval" = A1[["pval"]],
                                           "k" = summary(model.A)$k,
                                           "I2" = I2.A,
                                           stringsAsFactors = FALSE)
                       
                     
                       rbind(A_res, BA_res)[, colnames_final]
                     }
                   }

stopCluster(cl)

View(Results)

write.csv(Results,'/Users/yusha/Stable1b1-herbivore_Mn__.csv') #"Zone","NE_c"




###part3####
prefix <- "/Users/yusha/"
Data <- read.table("Alldata.txt", sep = "\t", header = TRUE)


Data <- Data[order(Data$S_ID), ]
Cod <- table(Data$S_ID)
esid <- unlist(lapply(Cod, function(x) seq(1, x, 1)))
Data$EsID <- esid


Whole_data <- Data[!is.na(Data$yi), ]




V <- tcrossprod(sqrt(Whole_data$CONTROL_SD))/sqrt(outer(Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2, Whole_data$CONTROL_R * (Whole_data$CONTROL_M)^2))+
  diag((Whole_data$TREATMENT_SD)/(Whole_data$TREATMENT_R * (Whole_data$TREATMENT_M)^2))
V_ind <- pmax(as.matrix(crossprod(lFormula(yi~1+(1|S_ID/CONTROL_M), Whole_data)$reTrms$Zt)) - 1, 0)
V <- V*V_ind

model.null <- rma.mv(yi, V,
                     mods = ~Trophic_group_response-1, 
                     random = list(~ 1|S_ID,~1|EsID), 
                     data = Whole_data, control = list(optimize = "BFGS"),
                     method = "ML")
I2 <- i2_ml(model.null)
Results <- data.frame("Trophic_group_response",
                      NA,
                      AIC(model.null),
                      logLik(model.null),
                      NA,
                      summary(model.null)$p,
                      NA,
                      summary(model.null)$k,
                      I2,
                      fix.empty.names = F)
Pre_vars <- c("CP_C3","CP_Ct","NCP_P","Std_","added_plant_species")

num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~",paste("Trophic_group_response",Pre_vars[i],sep="+"), "-1", sep = ""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group_response",Pre_vars[i], sep = "+"),
                                                   paste("+Trophic_group_response",Pre_vars[i],sep = ":"), "-1", sep = ""))
                     
                     Tem_data <- Whole_data
                     Tem_data <- Tem_data[!is.na(Tem_data[[Pre_vars[i]]]), ]  
                     
                     model.A <- rma.mv(yi, V, mods = formula.A,
                                       random = list(~1|S_ID,~1|EsID),
                                       data = Tem_data,control = list(optimize = "BFGS"),
                                       method = "ML")
                     
                     model.B <- rma.mv(yi, V, mods = formula.B,
                                       random = list(~1|S_ID,~1|EsID),control = list(optimize = "BFGS"),
                                       data = Tem_data,
                                       method = "ML")
                     
                     I2.A <- i2_ml(model.A)
                     I2.B <- i2_ml(model.B)
                     
                     # Compare model B vs A
                     if (length(model.A$beta) != length(model.B$beta)){
                       BA <- anova(model.B, model.A)
                       BA_res <- data.frame(
                         formula = as.character(formula.B)[-1],
                         Model = paste("1.", i, "A", sep = ""),
                         AIC = BA[["fit.stats.f"]][["AIC"]],
                         logLik = BA[["fit.stats.f"]][["ll"]],
                         LRT = BA[["LRT"]],
                         df = BA[["parms.f"]],
                         pval = BA[["pval"]],
                         k = summary(model.B)$k,
                         I2 = I2.B,
                         stringsAsFactors = FALSE
                       )
                     } else {
                       BA_res <- data.frame(
                         formula = as.character(formula.B)[-1],
                         Model = paste("1.", i, "A", sep = ""),
                         AIC = NA, logLik = NA, LRT = NA, df = NA, pval = NA,
                         k = summary(model.B)$k,
                         I2 = I2.B,
                         stringsAsFactors = FALSE
                       )
                     }
                     
                     # Also record model A info
                     A1 <- anova(model.A, model.null)
                     A_res <- data.frame(
                       formula = as.character(formula.A)[-1],
                       Model = "1",
                       AIC = A1[["fit.stats.f"]][["AIC"]],
                       logLik = A1[["fit.stats.f"]][["ll"]],
                       LRT = A1[["LRT"]],
                       df = A1[["parms.f"]],
                       pval = A1[["pval"]],
                       k = summary(model.A)$k,
                       I2 = I2.A,
                       stringsAsFactors = FALSE
                     )
                     
                     # 合并 A 模型和 B 模型的结果
                     rbind(A_res, BA_res)
                   }

View(Results)

write.csv(Results,'/Users/yusha/Stable1c.csv')




###part4####
Pre_vars <- c("Zone") #"Mn__","NE_c"
num_cores <- 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


colnames_final <- c("Model", "Step", "AIC", "logLik", "LRT", "df", "pval", "k", "I2")

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~ Trophic_group_response +", Pre_vars[i], "-1"))
                     
                     formula.B <- as.formula(paste("~ Trophic_group_response +", Pre_vars[i],
                                                   "+ Trophic_group_response:", Pre_vars[i], "-1"))
                     
                     Tem_data <- Whole_data
                     
                     if (Pre_vars[i] == "Zone") {
                       
                       
                       Tem_data.Clmz <- Tem_data[!is.na(Tem_data$Zone) & !is.na(Tem_data$Trophic_group_response), ]
                       
                      
                       cat(">>> subset data distribution：\n")
                       print(table(Tem_data.Clmz$Trophic_group_response, Tem_data.Clmz$Zone))
                       
                     
                       V.Clmz <- tcrossprod(sqrt(Tem_data.Clmz$CONTROL_SD)) / 
                         sqrt(outer(Tem_data.Clmz$CONTROL_R * (Tem_data.Clmz$CONTROL_M)^2,
                                    Tem_data.Clmz$CONTROL_R * (Tem_data.Clmz$CONTROL_M)^2)) +
                         diag((Tem_data.Clmz$TREATMENT_SD) / (Tem_data.Clmz$TREATMENT_R * (Tem_data.Clmz$TREATMENT_M)^2))
                       
                       V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1|S_ID/CONTROL_M), Tem_data.Clmz)$reTrms$Zt)) - 1, 0)
                       V.Clmz <- V.Clmz * V_ind
                       
              
                       model.A <- rma.mv(yi, V.Clmz, mods = formula.A,
                                         random = list(~1|S_ID, ~1|EsID),
                                         data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                         method = "ML")
                       
                 
                       model.B <- rma.mv(yi, V.Clmz, mods = formula.B,
                                         random = list(~1|S_ID, ~1|EsID),
                                         data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                         method = "ML")
                       
                     
                       model.Clmz <- rma.mv(yi, V.Clmz, mods = ~Trophic_group_response - 1,
                                            random = list(~1|S_ID, ~1|EsID),
                                            data = Tem_data.Clmz, control = list(optimize = "BFGS"),
                                            method = "ML")
                       
                     
                       A1 <- anova(model.A, model.Clmz)
                       I2.A <- i2_ml(model.A)
                       I2.B <- i2_ml(model.B)
                       
                     
                       if (length(model.A$beta) != length(model.B$beta)) {
                         BA <- anova(model.B, model.A)
                         BA_res <- data.frame("Model" = as.character(formula.B)[-1],
                                              "Step" = paste("1.", i, "A", sep = ""),
                                              "AIC" = BA[["fit.stats.f"]][["AIC"]],
                                              "logLik" = BA[["fit.stats.f"]][["ll"]],
                                              "LRT" = BA[["LRT"]],
                                              "df" = BA[["parms.f"]],
                                              "pval" = BA[["pval"]],
                                              "k" = summary(model.B)$k,
                                              "I2" = I2.B,
                                              stringsAsFactors = FALSE)
                       } else {
                         BA_res <- data.frame("Model" = as.character(formula.B)[-1],
                                              "Step" = paste("1.", i, "A", sep = ""),
                                              "AIC" = NA,
                                              "logLik" = NA,
                                              "LRT" = NA,
                                              "df" = NA,
                                              "pval" = NA,
                                              "k" = summary(model.B)$k,
                                              "I2" = I2.B,
                                              stringsAsFactors = FALSE)
                       }
                       
                       A_res <- data.frame("Model" = as.character(formula.A)[-1],
                                           "Step" = "1",
                                           "AIC" = A1[["fit.stats.f"]][["AIC"]],
                                           "logLik" = A1[["fit.stats.f"]][["ll"]],
                                           "LRT" = A1[["LRT"]],
                                           "df" = A1[["parms.f"]],
                                           "pval" = A1[["pval"]],
                                           "k" = summary(model.A)$k,
                                           "I2" = I2.A,
                                           stringsAsFactors = FALSE)
                       
                   
                       rbind(A_res, BA_res)[, colnames_final]
                     }
                   }

stopCluster(cl)

View(Results)

write.csv(Results,'/Users/yusha/Stable1d-Zone.csv')



####Supplementary Table 1.2####

library(dplyr)
library(stringr)
library(rotl)
library(ape)
library(Matrix)    #  bdiag
library(metafor)


Data <- read.table("Alldata.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


Plant_data <- Data[!is.na(Data$yi), ]


factor_cols <- c("S_ID","CP_C3","CP_Ct","Mn__","NE_c","NCP_P","Std_",
                 "Zone","Trophic_group","Trophic_group_response")

factor_cols_existing <- intersect(factor_cols, names(Plant_data))
Plant_data[factor_cols_existing] <- lapply(Plant_data[factor_cols_existing], factor)


num_cols <- c("CONTROL_SD", "CONTROL_R", "CONTROL_M",
              "TREATMENT_SD", "TREATMENT_R", "TREATMENT_M")
num_cols <- intersect(num_cols, names(Plant_data))
Plant_data[num_cols] <- lapply(Plant_data[num_cols], function(x) as.numeric(as.character(x)))



Plant_data$CP_Ltn_nm_f_ <- as.character(Plant_data$CP_Ltn_nm_f_)
Plant_data$Plant_Species_clean <- Plant_data$CP_Ltn_nm_f_ %>%
  str_replace_all("_", " ") %>%
  str_replace_all("×", "x") %>%
  str_replace_all("X", "x") %>%
  str_replace_all("[^A-Za-z\\sx]", " ") %>% 
  str_squish() %>%
  str_to_lower()


Plant_data$Plant_Species_clean <- recode(Plant_data$Plant_Species_clean,
                                         "brassica chinensis" = "brassica rapa subsp. chinensis",
                                         "citrus x aurantium" = "citrus aurantium",
                                         "fragaria x ananassa" = "fragaria vesca")


Species <- sort(unique(Plant_data$Plant_Species_clean))
message("Unique species after cleaning: ", length(Species))


# tnrs_match_names 
Taxa_Plant <- tnrs_match_names(Species)


ott_ids_all <- ott_id(Taxa_Plant)

in_tree_logical <- is_in_tree(ott_ids_all)


ott_used <- ott_ids_all[in_tree_logical]

if (length(ott_used) == 0) stop("No  OTT id of plant species can be found in OpenTree")

Plant_Tree <- tol_induced_subtree(ott_ids = ott_used, label_format = "id")

Plant_brlen <- compute.brlen(Plant_Tree)
Plant_Cor <- vcv(Plant_brlen, corr = TRUE)


coln <- colnames(Plant_Cor)
extract_ott <- function(label) {
num <- str_extract(label, "\\d+")
  return(num)
}
ott_labels <- sapply(coln, extract_ott)


tp_df <- as.data.frame(Taxa_Plant)

if (!"search_string" %in% colnames(tp_df)) {
  tp_df$search_string <- NA
}

if (!"ott_id" %in% colnames(tp_df)) {
  tp_df$ott_id <- ott_id(Taxa_Plant)
}


matched_search <- tp_df$search_string[match(ott_labels, tp_df$ott_id)]

colnames(Plant_Cor) <- ifelse(is.na(matched_search), coln, matched_search)
rownames(Plant_Cor) <- colnames(Plant_Cor)


missing_species <- setdiff(Species, colnames(Plant_Cor))
if (length(missing_species) > 0) {
  message("missing species：")
  print(missing_species)

  if (nrow(Plant_Cor) > 0) {
    Plant_Cor_extended <- as.matrix(bdiag(Plant_Cor, diag(1, length(missing_species))))
    old_names <- colnames(Plant_Cor)
    new_names <- c(old_names, missing_species)
    colnames(Plant_Cor_extended) <- new_names
    rownames(Plant_Cor_extended) <- new_names
    Plant_Cor <- as.matrix(Plant_Cor_extended)
  } else {

    Plant_Cor <- diag(1, length(missing_species))
    colnames(Plant_Cor) <- rownames(Plant_Cor) <- missing_species
  }
} else {
  message("All plant species are matched")
}


Plant_Cor <- Plant_Cor[order(colnames(Plant_Cor)), order(colnames(Plant_Cor))]


Plant_data$Plant_Species_clean <- str_to_lower(Plant_data$Plant_Species_clean)

Plant_data <- Plant_data %>% mutate(Plant_Species_clean = as.character(Plant_Species_clean))

not_in_cor <- unique(Plant_data$Plant_Species_clean[!(Plant_data$Plant_Species_clean %in% colnames(Plant_Cor))])
if (length(not_in_cor) > 0) {
  warning("the species below are not matched in Plant_Cor.please check the name or find other ways")
  print(not_in_cor)
  
  Plant_data <- Plant_data[Plant_data$Plant_Species_clean %in% colnames(Plant_Cor), ]
}

Plant_data$Plant_Species <- factor(Plant_data$Plant_Species_clean, levels = colnames(Plant_Cor))
Plant_data$Plant_Species2 <- Plant_data$Plant_Species


stopifnot(all(levels(Plant_data$Plant_Species) %in% colnames(Plant_Cor)))

stopifnot(sum(is.na(Plant_data$Plant_Species)) == 0)

message("Final no. of rows of Plant_Cor: ", nrow(Plant_Cor))
message("Plant_data: species number: ", length(levels(Plant_data$Plant_Species)))





V <- tcrossprod(sqrt(Plant_data$CONTROL_SD)) /
  sqrt(outer(Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2,
             Plant_data$CONTROL_R * (Plant_data$CONTROL_M)^2)) +
  diag((Plant_data$TREATMENT_SD)/(Plant_data$TREATMENT_R * (Plant_data$TREATMENT_M)^2))


V_ind <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1|S_ID/CONTROL_M), Plant_data)$reTrms$Zt)) - 1, 0)
V <- V * V_ind

###part 1####
model.null <- rma.mv(yi, V,
                     mods = ~Trophic_group - 1,
                     random = list(~1|S_ID, ~1|Plant_Species, ~1|EsID, ~1|Plant_Species2),
                     R = list(Plant_Species = Plant_Cor[as.character(Plant_data$Plant_Species), as.character(Plant_data$Plant_Species)]),
                     data = Plant_data,
                     control = list(optimize = "BFGS", stepadj = 0.5, maxiter = 10000),
                     method = "ML")

message("Null model done.")
print(summary(model.null))


num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# -------------------- Pre_vars category without NA --------------------#
Pre_vars <- c("CP_C3","CP_Ct","Std_","added_plant_species")#"CP_C3","CP_Ct","Std_","added_plant_species"


Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~",paste("Trophic_group",Pre_vars[i],sep="+"), "-1", sep = ""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group",Pre_vars[i], sep = "+"),
                                                   paste("+Trophic_group",Pre_vars[i],sep = ":"), "-1", sep = ""))
                     
                     Tem_data <- Plant_data
                     
                    
                     model.A <- rma.mv(yi, V, mods = formula.A,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species,
                                                                          Plant_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     model.B <- rma.mv(yi, V, mods = formula.B,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species,
                                                                          Plant_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     
                     I2.A <- tryCatch(i2_ml(model.A), error = function(e) NA)
                     I2.B <- tryCatch(i2_ml(model.B), error = function(e) NA)
                     
                    
                     logLik.A <- as.numeric(logLik(model.A))
                     logLik.B <- as.numeric(logLik(model.B))
                     logLik.null <- as.numeric(logLik(model.null))
                     
                 
                     LRT_A <- 2 * (logLik.A - logLik.null)
                     df_A <- length(coef(model.A)) - length(coef(model.null))
                     pval_A <- 1 - pchisq(LRT_A, df_A)
                     
                    
                     LRT_B <- 2 * (logLik.B - logLik.A)
                     df_B <- length(coef(model.B)) - length(coef(model.A))
                     pval_B <- 1 - pchisq(LRT_B, df_B)
                     
            
                     extract_i2_values <- function(i2_obj) {
                       if (is.list(i2_obj) && "random" %in% names(i2_obj)) {
                         random_effects <- i2_obj[["random"]]
                         list(
                           total = ifelse("total" %in% names(i2_obj), i2_obj[["total"]], NA),
                           S_ID = ifelse("S_ID" %in% names(random_effects), random_effects[["S_ID"]], NA),
                           Plant_Species = ifelse("Plant_Species" %in% names(random_effects), random_effects[["Plant_Species"]], NA),
                           EsID = ifelse("EsID" %in% names(random_effects), random_effects[["EsID"]], NA),
                           Plant_Species2 = ifelse("Plant_Species2" %in% names(random_effects), random_effects[["Plant_Species2"]], NA)
                         )
                       } else {
                         list(
                           total = NA,
                           S_ID = NA,
                           Plant_Species = NA,
                           EsID = NA,
                           Plant_Species2 = NA
                         )
                       }
                     }
                     
                     i2_A_values <- extract_i2_values(I2.A)
                     i2_B_values <- extract_i2_values(I2.B)
                     
                 
                     resA <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.A),
                       Model = "A",
                       AIC = AIC(model.A),
                       logLik = logLik.A,
                       LRT = LRT_A,  
                       df = df_A,   
                       pval = pval_A,
                       k = model.A$k,
                       I2_Total = i2_A_values$total,
                       I2_S_ID = i2_A_values$S_ID,
                       I2_Plant_Species = i2_A_values$Plant_Species,
                       I2_EsID = i2_A_values$EsID,
                       I2_Plant_Species2 = i2_A_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     resB <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.B),
                       Model = "B",
                       AIC = AIC(model.B),
                       logLik = logLik.B,
                       LRT = LRT_B,
                       df = df_B,
                       pval = pval_B,
                       k = model.B$k,
                       I2_Total = i2_B_values$total,
                       I2_S_ID = i2_B_values$S_ID,
                       I2_Plant_Species = i2_B_values$Plant_Species,
                       I2_EsID = i2_B_values$EsID,
                       I2_Plant_Species2 = i2_B_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     rbind(resA, resB)
                   }
stopCluster(cl)
View(Results)

write.csv(Results,'/Users/yusha/Desktop/文章提交/NC提交_Round 2/last_round/new_code/check/系统发育树/Stable1.2-Trophic_group-NCP_P.csv')#"CP_C3","CP_Ct","Std_","added_plant_species"





###part2####
num_cores <- 7
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# --------------------  Pre_vars category with NA--------------------#
Pre_vars <- c("Zone")#"NE_c","Mn__"

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~", paste("Trophic_group", Pre_vars[i], sep="+"), "-1", sep=""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group", Pre_vars[i], sep="+"),
                                                   paste("+Trophic_group", Pre_vars[i], sep=":"), "-1", sep=""))
                     
                   
                     Tem_data <- Plant_data[!is.na(Plant_data[[Pre_vars[i]]]), ]
                     
                    
                     V_temp <- tcrossprod(sqrt(Tem_data$CONTROL_SD)) /
                       sqrt(outer(Tem_data$CONTROL_R * (Tem_data$CONTROL_M)^2,
                                  Tem_data$CONTROL_R * (Tem_data$CONTROL_M)^2)) +
                       diag((Tem_data$TREATMENT_SD)/(Tem_data$TREATMENT_R * (Tem_data$TREATMENT_M)^2))
                     
                     V_ind_temp <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1|S_ID/CONTROL_M), 
                                                                     Tem_data)$reTrms$Zt)) - 1, 0)
                     V_temp <- V_temp * V_ind_temp
                     
                  
                  
                     model.A <- rma.mv(yi, V_temp, mods = formula.A,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species,
                                                                          Tem_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     model.B <- rma.mv(yi, V_temp, mods = formula.B,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species,
                                                                          Tem_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                    
                     I2.A <- tryCatch(i2_ml(model.A), error = function(e) NA)
                     I2.B <- tryCatch(i2_ml(model.B), error = function(e) NA)
                     I2.null <- tryCatch(i2_ml(model.null.temp), error = function(e) NA)
                     
                  
                     logLik.A <- as.numeric(logLik(model.A))
                     logLik.B <- as.numeric(logLik(model.B))
                     logLik.null <- as.numeric(logLik(model.null.temp))
                     
                 
                     LRT_A <- 2 * (logLik.A - logLik.null)
                     df_A <- length(coef(model.A)) - length(coef(model.null.temp))
                     pval_A <- 1 - pchisq(LRT_A, df_A)
                     
                   
                     LRT_B <- 2 * (logLik.B - logLik.A)
                     df_B <- length(coef(model.B)) - length(coef(model.A))
                     pval_B <- 1 - pchisq(LRT_B, df_B)
                     
                 
                     extract_i2_values <- function(i2_obj) {
                       if (is.list(i2_obj) && "random" %in% names(i2_obj)) {
                         random_effects <- i2_obj[["random"]]
                         list(
                           total = ifelse("total" %in% names(i2_obj), i2_obj[["total"]], NA),
                           S_ID = ifelse("S_ID" %in% names(random_effects), random_effects[["S_ID"]], NA),
                           Plant_Species = ifelse("Plant_Species" %in% names(random_effects), random_effects[["Plant_Species"]], NA),
                           EsID = ifelse("EsID" %in% names(random_effects), random_effects[["EsID"]], NA),
                           Plant_Species2 = ifelse("Plant_Species2" %in% names(random_effects), random_effects[["Plant_Species2"]], NA)
                         )
                       } else {
                         list(
                           total = NA,
                           S_ID = NA,
                           Plant_Species = NA,
                           EsID = NA,
                           Plant_Species2 = NA
                         )
                       }
                     }
                     
                     i2_A_values <- extract_i2_values(I2.A)
                     i2_B_values <- extract_i2_values(I2.B)
                     i2_null_values <- extract_i2_values(I2.null)
                     
                   
                     resA <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.A),
                       Model = "A",
                       AIC = AIC(model.A),
                       logLik = logLik.A,
                       LRT = LRT_A, 
                       df = df_A, 
                       pval = pval_A,
                       k = model.A$k,
                       I2_Total = i2_A_values$total,
                       I2_S_ID = i2_A_values$S_ID,
                       I2_Plant_Species = i2_A_values$Plant_Species,
                       I2_EsID = i2_A_values$EsID,
                       I2_Plant_Species2 = i2_A_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     resB <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.B),
                       Model = "B",
                       AIC = AIC(model.B),
                       logLik = logLik.B,
                       LRT = LRT_B,
                       df = df_B,
                       pval = pval_B,
                       k = model.B$k,
                       I2_Total = i2_B_values$total,
                       I2_S_ID = i2_B_values$S_ID,
                       I2_Plant_Species = i2_B_values$Plant_Species,
                       I2_EsID = i2_B_values$EsID,
                       I2_Plant_Species2 = i2_B_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     rbind(resA, resB)
                   }


stopCluster(cl)

View(Results)

write.csv(Results,'/Users/yusha/Stable1.2-Trophic_group-Zone.csv')#"NE_c","Mn__"



###part 3####
model.null <- rma.mv(yi, V,
                     mods = ~Trophic_group_response - 1,
                     random = list(~1|S_ID, ~1|Plant_Species, ~1|EsID, ~1|Plant_Species2),
                     R = list(Plant_Species = Plant_Cor[as.character(Plant_data$Plant_Species), as.character(Plant_data$Plant_Species)]),
                     data = Plant_data,
                     control = list(optimize = "BFGS", stepadj = 0.5, maxiter = 10000),
                     method = "ML")

message("Null model done.")
print(summary(model.null))


num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# -------------------- Pre_vars category without NA --------------------#
Pre_vars <- c("NCP_P")#"CP_C3","CP_Ct","Std_","added_plant_species"


Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~",paste("Trophic_group_response",Pre_vars[i],sep="+"), "-1", sep = ""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group_response",Pre_vars[i], sep = "+"),
                                                   paste("+Trophic_group_response",Pre_vars[i],sep = ":"), "-1", sep = ""))
                     
                     Tem_data <- Plant_data
                     
                     
                     model.A <- rma.mv(yi, V, mods = formula.A,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species,
                                                                          Plant_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     model.B <- rma.mv(yi, V, mods = formula.B,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Plant_data$Plant_Species,
                                                                          Plant_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     
                     I2.A <- tryCatch(i2_ml(model.A), error = function(e) NA)
                     I2.B <- tryCatch(i2_ml(model.B), error = function(e) NA)
                     
                     
                     logLik.A <- as.numeric(logLik(model.A))
                     logLik.B <- as.numeric(logLik(model.B))
                     logLik.null <- as.numeric(logLik(model.null))
                     
                     
                     LRT_A <- 2 * (logLik.A - logLik.null)
                     df_A <- length(coef(model.A)) - length(coef(model.null))
                     pval_A <- 1 - pchisq(LRT_A, df_A)
                     
                     
                     LRT_B <- 2 * (logLik.B - logLik.A)
                     df_B <- length(coef(model.B)) - length(coef(model.A))
                     pval_B <- 1 - pchisq(LRT_B, df_B)
                     
                     
                     extract_i2_values <- function(i2_obj) {
                       if (is.list(i2_obj) && "random" %in% names(i2_obj)) {
                         random_effects <- i2_obj[["random"]]
                         list(
                           total = ifelse("total" %in% names(i2_obj), i2_obj[["total"]], NA),
                           S_ID = ifelse("S_ID" %in% names(random_effects), random_effects[["S_ID"]], NA),
                           Plant_Species = ifelse("Plant_Species" %in% names(random_effects), random_effects[["Plant_Species"]], NA),
                           EsID = ifelse("EsID" %in% names(random_effects), random_effects[["EsID"]], NA),
                           Plant_Species2 = ifelse("Plant_Species2" %in% names(random_effects), random_effects[["Plant_Species2"]], NA)
                         )
                       } else {
                         list(
                           total = NA,
                           S_ID = NA,
                           Plant_Species = NA,
                           EsID = NA,
                           Plant_Species2 = NA
                         )
                       }
                     }
                     
                     i2_A_values <- extract_i2_values(I2.A)
                     i2_B_values <- extract_i2_values(I2.B)
                     
                     
                     resA <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.A),
                       Model = "A",
                       AIC = AIC(model.A),
                       logLik = logLik.A,
                       LRT = LRT_A,  
                       df = df_A,   
                       pval = pval_A,
                       k = model.A$k,
                       I2_Total = i2_A_values$total,
                       I2_S_ID = i2_A_values$S_ID,
                       I2_Plant_Species = i2_A_values$Plant_Species,
                       I2_EsID = i2_A_values$EsID,
                       I2_Plant_Species2 = i2_A_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     resB <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.B),
                       Model = "B",
                       AIC = AIC(model.B),
                       logLik = logLik.B,
                       LRT = LRT_B,
                       df = df_B,
                       pval = pval_B,
                       k = model.B$k,
                       I2_Total = i2_B_values$total,
                       I2_S_ID = i2_B_values$S_ID,
                       I2_Plant_Species = i2_B_values$Plant_Species,
                       I2_EsID = i2_B_values$EsID,
                       I2_Plant_Species2 = i2_B_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     rbind(resA, resB)
                   }
stopCluster(cl)
View(Results)

write.csv(Results,'/Users/yusha/Stable1.2-Trophic_group_response-category without NA-NCP_P.csv')#"CP_C3","CP_Ct","Std_","added_plant_species"





###part4####
num_cores <- 7
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# --------------------  Pre_vars category with NA--------------------#
Pre_vars <- c("Zone")#"NE_c","Mn__"

Results <- foreach(i = 1:length(Pre_vars), .combine = "rbind",
                   .packages = c("metafor","ape","orchaRd","lme4")) %dopar% {
                     
                     formula.A <- as.formula(paste("~", paste("Trophic_group_response", Pre_vars[i], sep="+"), "-1", sep=""))
                     formula.B <- as.formula(paste("~", paste("Trophic_group_response", Pre_vars[i], sep="+"),
                                                   paste("+Trophic_group_response", Pre_vars[i], sep=":"), "-1", sep=""))
                     
                     
                     Tem_data <- Plant_data[!is.na(Plant_data[[Pre_vars[i]]]), ]
                     
                     
                     V_temp <- tcrossprod(sqrt(Tem_data$CONTROL_SD)) /
                       sqrt(outer(Tem_data$CONTROL_R * (Tem_data$CONTROL_M)^2,
                                  Tem_data$CONTROL_R * (Tem_data$CONTROL_M)^2)) +
                       diag((Tem_data$TREATMENT_SD)/(Tem_data$TREATMENT_R * (Tem_data$TREATMENT_M)^2))
                     
                     V_ind_temp <- pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1|S_ID/CONTROL_M), 
                                                                     Tem_data)$reTrms$Zt)) - 1, 0)
                     V_temp <- V_temp * V_ind_temp
                     
                     
                     
                     model.A <- rma.mv(yi, V_temp, mods = formula.A,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species,
                                                                          Tem_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     model.B <- rma.mv(yi, V_temp, mods = formula.B,
                                       random = list(~1|S_ID,~1|Plant_Species,~1|EsID,~1|Plant_Species2),
                                       R = list(Plant_Species = Plant_Cor[Tem_data$Plant_Species,
                                                                          Tem_data$Plant_Species]),
                                       data = Tem_data, control = list(optimize = "BFGS"), method = "ML")
                     
                     
                     I2.A <- tryCatch(i2_ml(model.A), error = function(e) NA)
                     I2.B <- tryCatch(i2_ml(model.B), error = function(e) NA)
                     I2.null <- tryCatch(i2_ml(model.null.temp), error = function(e) NA)
                     
                     
                     logLik.A <- as.numeric(logLik(model.A))
                     logLik.B <- as.numeric(logLik(model.B))
                     logLik.null <- as.numeric(logLik(model.null.temp))
                     
                     
                     LRT_A <- 2 * (logLik.A - logLik.null)
                     df_A <- length(coef(model.A)) - length(coef(model.null.temp))
                     pval_A <- 1 - pchisq(LRT_A, df_A)
                     
                     
                     LRT_B <- 2 * (logLik.B - logLik.A)
                     df_B <- length(coef(model.B)) - length(coef(model.A))
                     pval_B <- 1 - pchisq(LRT_B, df_B)
                     
                     
                     extract_i2_values <- function(i2_obj) {
                       if (is.list(i2_obj) && "random" %in% names(i2_obj)) {
                         random_effects <- i2_obj[["random"]]
                         list(
                           total = ifelse("total" %in% names(i2_obj), i2_obj[["total"]], NA),
                           S_ID = ifelse("S_ID" %in% names(random_effects), random_effects[["S_ID"]], NA),
                           Plant_Species = ifelse("Plant_Species" %in% names(random_effects), random_effects[["Plant_Species"]], NA),
                           EsID = ifelse("EsID" %in% names(random_effects), random_effects[["EsID"]], NA),
                           Plant_Species2 = ifelse("Plant_Species2" %in% names(random_effects), random_effects[["Plant_Species2"]], NA)
                         )
                       } else {
                         list(
                           total = NA,
                           S_ID = NA,
                           Plant_Species = NA,
                           EsID = NA,
                           Plant_Species2 = NA
                         )
                       }
                     }
                     
                     i2_A_values <- extract_i2_values(I2.A)
                     i2_B_values <- extract_i2_values(I2.B)
                     i2_null_values <- extract_i2_values(I2.null)
                     
                     
                     resA <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.A),
                       Model = "A",
                       AIC = AIC(model.A),
                       logLik = logLik.A,
                       LRT = LRT_A, 
                       df = df_A, 
                       pval = pval_A,
                       k = model.A$k,
                       I2_Total = i2_A_values$total,
                       I2_S_ID = i2_A_values$S_ID,
                       I2_Plant_Species = i2_A_values$Plant_Species,
                       I2_EsID = i2_A_values$EsID,
                       I2_Plant_Species2 = i2_A_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     resB <- data.frame(
                       Variable = Pre_vars[i],
                       Formula = deparse(formula.B),
                       Model = "B",
                       AIC = AIC(model.B),
                       logLik = logLik.B,
                       LRT = LRT_B,
                       df = df_B,
                       pval = pval_B,
                       k = model.B$k,
                       I2_Total = i2_B_values$total,
                       I2_S_ID = i2_B_values$S_ID,
                       I2_Plant_Species = i2_B_values$Plant_Species,
                       I2_EsID = i2_B_values$EsID,
                       I2_Plant_Species2 = i2_B_values$Plant_Species2,
                       stringsAsFactors = FALSE
                     )
                     
                     rbind(resA, resB)
                   }


stopCluster(cl)

View(Results)

write.csv(Results,'/Users/yusha/Stable1.2-Trophic_group_response-Zone.csv')#"NE_c","Mn__"