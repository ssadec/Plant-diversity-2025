rm(list = ls())
library(readxl)
library(dplyr)
library(metafor)



prefix <- "/Users/yusha/"

sheets <- c(
  "Predator_reproduction_response","Predator_diversity_response","Predator_predation_response",
  "Parasitoid_growth_response","Parasitoid_reproduction_response","Parasitoid_diversity_response",
  "Parasitoid_parasitism_response","NEs_reproduction_response","NEs_diversity_response",
  "Herbivore_growth_response","Herbivore_reproduction_response","Herbivore_damage_response",
  "Crop_growth_response","Crop_reproduction_response","Crop_quality_response"
)


geary_eq13 <- function(mean, sd, n){
  ifelse(is.na(mean) | is.na(sd) | is.na(n) | sd == 0, NA,
         (mean / sd) * (4 * n^(3/2) / (1 + 4*n)))
}


all_data <- list()
for (i in seq_along(sheets)) {
  df <- read_xlsx(paste0(prefix, "data of all trophic groups_2025_1001.xlsx"),
                  sheet = i, n_max = 1200)[,1:41]
  colnames(df) <- make.names(colnames(df))
  all_data[[sheets[i]]] <- df
}
colnames_abbr <- abbreviate(colnames(all_data[["Crop_growth_response"]]), minlength = 4)


Alldata_full  <- NULL
Alldata_geary <- NULL
Antadata_full <- NULL
Antadata_geary <- NULL
NEsdata_full  <- NULL
NEsdata_geary <- NULL


for (i in sheets) {
  
  data_i <- all_data[[i]]
  colnames(data_i) <- colnames_abbr
  data_i <- subset(data_i, !is.na(CONTROL_M))
  

  colnames(data_i)[colnames(data_i) == "TREATMENT_MEAN"] <- "TREATMENT_M"
  colnames(data_i)[colnames(data_i) == "TREATMENT_SD"]   <- "TREATMENT_SD"
  colnames(data_i)[colnames(data_i) == "TREATMENT_REPLICATION_NUMBER"] <- "TREATMENT_R"
  colnames(data_i)[colnames(data_i) == "CONTROL_MEAN"] <- "CONTROL_M"
  colnames(data_i)[colnames(data_i) == "CONTROL_SD"]   <- "CONTROL_SD"
  colnames(data_i)[colnames(data_i) == "CONTROL_REPLICATION_NUMBER"] <- "CONTROL_R"
  

  ESData <- data_i %>%
    mutate(
      cv_Control   = na_if(CONTROL_SD / CONTROL_M, Inf),
      cv_Treatment = na_if(TREATMENT_SD / TREATMENT_M, Inf)
    )
  
  ESData <- cv_avg(x = TREATMENT_M, sd = TREATMENT_SD,
                   n = TREATMENT_R, group = S_ID, label="1", data=ESData)
  ESData <- cv_avg(x = CONTROL_M, sd = CONTROL_SD,
                   n = CONTROL_R, group = S_ID, label="2", data=ESData)
  
  ESData <- ESData %>%
    mutate(
      cv2_trea_new = if_else(is.na(cv_Treatment), b_CV2_1, cv_Treatment^2),
      cv2_cont_new = if_else(is.na(cv_Control),   b_CV2_2, cv_Control^2)
    )
  
  ESData$TREATMENT_SD[is.na(ESData$TREATMENT_SD)] <- 
    ESData$cv2_trea_new[is.na(ESData$TREATMENT_SD)] * ESData$TREATMENT_M[is.na(ESData$TREATMENT_SD)]
  ESData$CONTROL_SD[is.na(ESData$CONTROL_SD)] <- 
    ESData$cv2_cont_new[is.na(ESData$CONTROL_SD)] * ESData$CONTROL_M[is.na(ESData$CONTROL_SD)]
  

  ESData <- ESData %>%
    mutate(
      yi = lnrr_laj(m1=TREATMENT_M, m2=CONTROL_M, cv1_2=b_CV2_1, cv2_2=b_CV2_2, n1=TREATMENT_R, n2=CONTROL_R),
      vi = v_lnrr_laj(cv1_2=b_CV2_1, n1=TREATMENT_R, cv2_2=b_CV2_2, n2=CONTROL_R)
    )
  
  # ---- Geary test ----
  ESData <- ESData %>%
    mutate(
      geary_control   = geary_eq13(CONTROL_M, CONTROL_SD, CONTROL_R),
      geary_treatment = geary_eq13(TREATMENT_M, TREATMENT_SD, TREATMENT_R),
      pass_control    = ifelse(!is.na(geary_control) & geary_control >= 3, TRUE, FALSE),
      pass_treatment  = ifelse(!is.na(geary_treatment) & geary_treatment >= 3, TRUE, FALSE),
      pass_geary      = pass_control & pass_treatment
    )
  

  S_IDtab <- table(ESData$S_ID)
  esid <- unlist(lapply(S_IDtab, function(n) seq_len(n)))
  ESData$EsID <- esid
  

  trophic <- regmatches(i, regexec("(Predator)?(Parasitoid)?(NEs)?(Herbivore)?(Crop)?", i))[[1]][1]
  ESData$Trophic_group <- trophic
  ESData$Trophic_group_response <- i
  

  if (trophic == "Herbivore") {
    Antadata_full <- rbind(Antadata_full, ESData)
  } else if (trophic %in% c("Predator", "Parasitoid", "NEs")) {
    NEsdata_full <- rbind(NEsdata_full, ESData)
  }
  
  
  

  

  Alldata_full <- rbind(Alldata_full, ESData)
  if (grepl("^Herbivore", i)) Antadata_full <- rbind(Antadata_full, ESData)
  if (grepl("Predator|Parasitoid|NEs", i)) NEsdata_full <- rbind(NEsdata_full, ESData)
  

  ESData_geary <- ESData %>% filter(pass_geary)
  if ("NCP_N" %in% colnames(ESData_geary)) {
    ESData_geary <- ESData_geary %>%
      mutate(added_plant_species = ifelse(!is.na(NCP_N) & NCP_N > 0, log2(NCP_N), NA))
  }
  
  Alldata_geary <- rbind(Alldata_geary, ESData_geary)
  if (grepl("^Herbivore", i)) Antadata_geary <- rbind(Antadata_geary, ESData_geary)
  if (grepl("Predator|Parasitoid|NEs", i)) NEsdata_geary <- rbind(NEsdata_geary, ESData_geary)
  

  write.table(ESData, paste0(i, "_full.txt"), sep="\t", row.names=FALSE)
  write.table(ESData_geary, paste0(i, ".txt"), sep="\t", row.names=FALSE)
}



write.table(Alldata_full,  "Alldata_full.txt",  sep="\t", row.names=FALSE)
write.table(Antadata_full, "Antadata_full.txt", sep="\t", row.names=FALSE)
write.table(NEsdata_full,  "NEsdata_full.txt",  sep="\t", row.names=FALSE)

write.table(Alldata_geary,  "Alldata.txt",  sep="\t", row.names=FALSE)
write.table(Antadata_geary, "Antadata.txt", sep="\t", row.names=FALSE)
write.table(NEsdata_geary,  "NEsdata.txt",  sep="\t", row.names=FALSE)



data_names <- sheets
Total_response_list <- list(
  Predator   = list(),
  Parasitoid = list(),
  NEs        = list(),
  Herbivore  = list(),
  Crop       = list()
)

for (i in data_names) {
  geary_file_path <- paste0(prefix, i, ".txt")
  if (!file.exists(geary_file_path)) {
    warning("Geary file not found: ", geary_file_path)
    next
  }
  
  ESData_geary <- read.table(geary_file_path, sep = "\t", header = TRUE)
  trophic <- strsplit(i, "_")[[1]][1]
  
  if (trophic == "Crop") {
    Total_response_list$Crop <- bind_rows(Total_response_list$Crop, ESData_geary)
  } else if (trophic == "Herbivore") {
    Total_response_list$Herbivore <- bind_rows(Total_response_list$Herbivore, ESData_geary)
  } else if (trophic == "Predator") {
    Total_response_list$Predator <- bind_rows(Total_response_list$Predator, ESData_geary)
  } else if (trophic == "Parasitoid") {
    Total_response_list$Parasitoid <- bind_rows(Total_response_list$Parasitoid, ESData_geary)
  } else if (trophic == "NEs") {
    Total_response_list$NEs <- bind_rows(Total_response_list$NEs, ESData_geary)
  } else {
    warning("Unrecognized trophic group in dataset: ", i)
  }
}


write.table(Total_response_list$Predator,   paste0(prefix, "Total_response_of_predator.txt"), sep="\t", row.names=FALSE)
write.table(Total_response_list$Parasitoid, paste0(prefix, "Total_response_of_parasitoid.txt"), sep="\t", row.names=FALSE)
write.table(Total_response_list$NEs,        paste0(prefix, "Total_response_of_NEs.txt"), sep="\t", row.names=FALSE)
write.table(Total_response_list$Herbivore,  paste0(prefix, "Total_response_of_herbivore.txt"), sep="\t", row.names=FALSE)
write.table(Total_response_list$Crop,       paste0(prefix, "Total_response_of_crop.txt"), sep="\t", row.names=FALSE)
