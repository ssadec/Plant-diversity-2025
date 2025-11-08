library(nlme)
library(piecewiseSEM)
library(readr)  # safer reading for tsv
library(rsq)
prefix <- "/Users/yusha/"


total_names_df <- read_tsv(paste0(prefix, "total_names.txt"), col_names = TRUE)
total_names <- total_names_df$total_name 


Epdt <- read_tsv(paste0(prefix, "Alldata.txt"), col_names = TRUE)
Epdt$S_ID <- paste0("Athr", 1:nrow(Epdt))


res <- NULL


Category <- c("Predator performance", "Parasitoid performance", 
              "NEs performance", "Herbivore performance", "Crop performance")
Strat <- c("Global", "Herbaceous_plant", "Woody_plant", "Food_crop", "Cash_crop",
           "Mn__Generalist", "Mn__Specialist",
           "NE_c_Generalist", "NE_c_Specialist",
           "Intercropping", "Cover_cropping", "Sown_field_margins",
           "Natural","Semi_natural", "Controlled", "Temperate", "Tropic")
Fignum <- paste(rep(1:17, each = length(total_names)), rep(letters[1:length(total_names)], 17), sep = "")
Figsource <- paste("Supplementary Data Fig.", Fignum, sep = "")


col_mapping <- list(
  S_ID = "S_ID",
  yi = "yi",
  vi = "vi",
  NCP_N = "NCP_N" 
)

k <- 3
for (j in Strat) {
  for (i in seq_along(total_names)) {
    dat_path <- paste0(prefix, total_names[i], ".txt")
    if (!file.exists(dat_path)) {
      warning("file does not exist：", dat_path)
      next
    }
    
  
    dat <- read_tsv(dat_path, col_types = cols(.default = "c"))
    

    for (std_col in names(col_mapping)) {
      if (!(col_mapping[[std_col]] %in% colnames(dat))) {
        warning("lack of important column, skip：", dat_path)
        next
      } else {
        dat[[std_col]] <- dat[[ col_mapping[[std_col]] ]]
      }
    }
    

    dat$yi <- as.numeric(dat$yi)
    dat$vi <- as.numeric(dat$vi)
    dat$NCP_N <- as.numeric(dat$NCP_N)
    

    dat$added_plant_species <- ifelse(!is.na(dat$NCP_N) & dat$NCP_N > 0, log2(dat$NCP_N), NA)
    

    dat$Epdt <- Epdt[match(dat$S_ID, Epdt$S_ID), 2]
    

    dat <- dat[!is.na(dat$vi), ]
    
  
    if (j == "Global") {
      dat_sub <- dat
    } else if (j %in% c("Mn__Generalist","Mn__Specialist")) {
      if (Category[i] == "Herbivore performance") {
        target <- sub("Mn__","", j)
        dat_sub <- subset(dat, Mn__ == target)
      } else {
        dat_sub <- dat[0,]
      }
    } else if (j %in% c("NE_c_Generalist","NE_c_Specialist")) {
      if (Category[i] %in% c("Predator performance","Parasitoid performance","NEs performance")) {
        target <- sub("NE_c_","", j)
        dat_sub <- subset(dat, NE_c == target)
      } else {
        dat_sub <- dat[0,]
      }
    } else {
    
      dat_sub <- subset(dat,
                        CP_C3 == j | CP_Ct == j | Mn__ == j | NE_c == j |
                          NCP_P == j | Std_ == j | Zone == j)
    }
    

    if (nrow(dat_sub) <= 1 || length(unique(dat_sub$added_plant_species)) == 1) {
      res <- rbind(res, data.frame(Figsource=Figsource[k], Category=Category[i], Strat=j,
                                   n=nrow(dat_sub), n_study=length(unique(dat_sub$S_ID)),
                                   regression_eq=NA, slope=NA, lower=NA, upper=NA,
                                   se=NA, df=NA, t_val=NA, p_val=NA, R2=NA))
    } else if (length(unique(dat_sub$S_ID)) == 1) {
      fit <- lm(yi ~ added_plant_species, data=dat_sub, weights=1/vi)
      coefs <- coef(summary(fit))
      intercept <- coefs[1,1]; slope <- coefs[2,1]
      regression_eq <- paste0("Y = ", round(intercept,4), ifelse(slope>=0," + "," - "), abs(round(slope,4)), "X")
      res <- rbind(res, data.frame(Figsource=Figsource[k], Category=Category[i], Strat=j,
                                   n=nrow(dat_sub), n_study=length(unique(dat_sub$S_ID)),
                                   regression_eq=regression_eq,
                                   slope=slope,
                                   lower=slope - coefs[2,2]*qt(0.975, fit$df.residual),
                                   upper=slope - coefs[2,2]*qt(0.025, fit$df.residual),
                                   se=coefs[2,2], df=fit$df.residual,
                                   t_val=coefs[2,3], p_val=coefs[2,4], R2=NA))
    } else {
      fit <- lme(yi ~ added_plant_species, random=~1|S_ID, data=dat_sub, weights=varFixed(~vi))
      coefs <- coef(summary(fit))
      intercept <- coefs[1,1]; slope <- coefs[2,1]
      regression_eq <- paste0("Y = ", round(intercept,4), ifelse(slope>=0," + "," - "), abs(round(slope,4)), "X")
      res <- rbind(res, data.frame(Figsource=Figsource[k], Category=Category[i], Strat=j,
                                   n=nrow(dat_sub), n_study=length(unique(dat_sub$S_ID)),
                                   regression_eq=regression_eq,
                                   slope=slope,
                                   lower=slope - coefs[2,2]*qt(0.975, coefs[2,3]),
                                   upper=slope - coefs[2,2]*qt(0.025, coefs[2,3]),
                                   se=coefs[2,2], df=coefs[2,3],
                                   t_val=coefs[2,4], p_val=coefs[2,5],
                                   R2=rsq.lmm(fit)$model))
    }
    k <- k+1
  }
}



View(res)
write.csv(res,'/Users/yusha/Stable17.csv')
