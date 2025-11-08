
rm(list = ls())

library(nlme)
library(effects)
library(scales)

# 📁 path
prefix <- "/Users/yusha/"
prefix2 <- "/Users/yusha/Desktop/check/"

# 📊 data
total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[, 1]
Epdt <- read.table(paste0(prefix, "Alldata.txt"), sep = "\t", header = TRUE)
Epdt$S_ID <- paste0("Athr", 1:nrow(Epdt))

res <- read.csv(paste0(prefix2, "Stable17.csv"), na.strings = c("NA", ""))
colnames(res) <- gsub("\\.", "_", trimws(colnames(res)))
View(res)
# labs
ylabs <- c("Predator performance", "Parasitoid performance", "NEs performance",
           "Herbivore performance", "Crop performance")

cat_colors <- c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73")
names(cat_colors) <- c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop")

# strat
Strat <- c("Global")

k <- 1
all_plots <- list()
plot_labels <- letters
label_index <- 1

# --- MODIFIED ---
safe_condition <- function(pval, dat_sub) {
  if (is.null(pval)) return(FALSE)
  if (length(pval) != 1) return(FALSE)
  if (is.na(pval) || !is.finite(pval)) return(FALSE)
  if (pval >= 0.05) return(FALSE)
  if (length(unique(na.omit(dat_sub$added_plant_species))) <= 1) return(FALSE)
  return(TRUE)
}
# --- END MODIFIED ---

for (j in seq_along(Strat)) {
  for (i in seq_along(ylabs)) {
    # read file
    dat_path <- paste0(prefix, total_names[i], ".txt")
    if (!file.exists(dat_path)) {
      warning("file doesnot exit, skip: ", dat_path)
      next
    }
    
    dat <- read.table(dat_path, sep = "\t", header = TRUE)
    
    if ("S_ID" %in% colnames(dat)) {
      dat$Epdt <- Epdt[match(dat$S_ID, Epdt$S_ID), 2]
    } else {
      dat$Epdt <- NA
    }
    
    dat <- dat[!is.na(dat$vi), ]
    
   
    dat_sub <- dat
    
    # 🔹added_plant_species = log2(NCP_N + 1)
    # --- MODIFIED: use log2(NCP_N + 1) to avoid log2(0) issues ---
    if ("NCP_N" %in% colnames(dat_sub)) {
      # note: if NCP_N can be non-numeric strings, coerce safely
      dat_sub$NCP_N <- suppressWarnings(as.numeric(as.character(dat_sub$NCP_N)))
      dat_sub$added_plant_species <- ifelse(!is.na(dat_sub$NCP_N),
                                               log2(dat_sub$NCP_N),
                                               NA)
    } else {
      dat_sub$added_plant_species <- NA
    }
    # --- END MODIFIED ---
    
    # --- MODIFIED: only keep rows that have added_plant_species, yi, vi (valid_rows)
    valid_rows <- dat_sub[complete.cases(dat_sub[, c("added_plant_species", "yi", "vi")]), ]
    if (nrow(valid_rows) == 0) {
      k <- k + 1
      next
    }
    # replace dat_sub with valid_rows so subsequent operations align
    dat_sub <- valid_rows
    # --- END MODIFIED ---
    
    # point size（dat_sub$vi）
    wi <- 1 / sqrt(dat_sub$vi)
    size <- if (length(wi) == 1 || max(wi) == min(wi)) {
      rep(3, length(wi))
    } else {
      1.5 + 5 * (wi - min(wi)) / (max(wi) - min(wi))
    }
    dat_sub$size <- size
    
   
    # if res has a column named "Stratification" use it; otherwise match only by Category
    if ("Strat" %in% colnames(res)) {
      res_row <- subset(res, Strat == Strat[j] & Category == ylabs[i])
    } else {
      res_row <- subset(res, Category == ylabs[i])
    }
    # --- END MODIFIED ---
    
    if (nrow(res_row) == 1) {
      eq <- as.character(res_row$Regression_equation)
      pval <- as.numeric(as.character(res_row$P_value))
      obs <- res_row$n
      stud <- res_row$n_study
    } else {
      eq <- NA; pval <- NA; obs <- NA; stud <- NA
    }
    obs_label <- paste0("(", obs, "/", stud, ")")
    
    # Category-color
    dat_sub$Category_short <- factor(ylabs[i],
                                     levels = ylabs,
                                     labels = c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop"))
    
    # save results
    all_plots[[length(all_plots) + 1]] <- list(
      data = dat_sub,
      eq = eq,
      pval = pval,
      obs_label = obs_label,
      label = plot_labels[label_index],
      ylab = ylabs[i]
    )
    
    label_index <- label_index + 1
    k <- k + 1
  }
}

#pics

panel_w <- 8   
panel_h <- 10  
nrow <- 3    
ncol <- 5     

png("Global_All.png", units = "cm", width = 40, height = 30, res = 300)

par(mai = c(1, 0.8, 0.15, 0.3), 
    mfrow = c(nrow, ncol), 
    mgp = c(3, 0.6, 0), 
    family = "Arial")

for (i in seq_along(all_plots)) {
  p <- all_plots[[i]]
  dat_sub <- p$data
  
  plot(dat_sub$added_plant_species, dat_sub$yi, type = "n",
       xlab = "",
       ylab = p$ylab,
       ylim = c(-8, 8), xlim = c(0, 6),
       cex.lab = 1.5, font.lab = 1, cex.axis = 1.2)
  mtext(expression(atop("Log"[2]*" Number of added", "plant species")),
        side = 1, line = 5, cex = 1, family = "Arial") 
  
  # --- MODIFIED: use safe_condition() to avoid NA/length issues in if ---
  if (safe_condition(p$pval, dat_sub)) {
    
    # choose fit type depending on number of studies
    if (length(unique(dat_sub$S_ID)) == 1) {
      fit <- lm(yi ~ added_plant_species, data = dat_sub, weights = 1 / vi)
    } else {
      # wrap lme in tryCatch to avoid singularity crashes
      fit <- tryCatch(
        lme(yi ~ added_plant_species, data = dat_sub,
            random = ~1 | S_ID, weights = varFixed(~vi),
            control = lmeControl(sigma = 1)),
        error = function(e) {
          message("lme failed（skip）: ", conditionMessage(e))
          return(NULL)
        }
      )
    }
    
    if (!is.null(fit)) {
      # compute effect prediction safely
      pred <- tryCatch(
        Effect("added_plant_species", fit,
               xlevels = list(added_plant_species = seq(0, 6, 0.01))),
        error = function(e) {
          message("Effect failed: ", conditionMessage(e))
          return(NULL)
        }
      )
      
      if (!is.null(pred)) {
        polygon(c(rev(pred$x[, 1]), pred$x[, 1]),
                c(rev(pred$upper[, 1]), pred$lower[, 1]),
                col = alpha("grey", 0.5), border = NA)
        lines(pred$x[, 1], pred$fit[, 1], col = "black", lwd = 3)
        
        if (!is.na(p$eq)) {
          text(5.5, 6.8, labels = p$eq, cex = 1.2, adj = 1)
        }
        text(5.5, 6.0, labels = paste("P =", format(signif(p$pval, 4), scientific = TRUE)), cex = 1.2, adj = 1)
      }
    }
  } else {
    
    if (!is.null(p$pval) && length(p$pval) == 1 && !is.na(p$pval)) {
      text(5.5, 6.0, labels = paste("P =", format(signif(p$pval, 4), scientific = TRUE)), cex = 1.0, adj = 1)
    }
  }
  # --- END MODIFIED ---
  
  points(dat_sub$added_plant_species, dat_sub$yi,
         bg = cat_colors[as.character(dat_sub$Category_short)], col = "black", pch = 21, cex = dat_sub$size)
  
  text(0.4, 7, labels = p$label, cex = 1.5, font = 2)
  
  text(5.5, 5.1, labels = p$obs_label, cex = 1.1, adj = 1)
}


dev.off()

