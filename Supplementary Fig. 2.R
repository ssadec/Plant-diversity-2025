rm(list = ls())
library(nlme)
library(effects)
library(scales)
library(readxl)


prefix <- "/Users/yusha/"
prefix2 <- "/Users/yusha/Desktop/check/"


total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[, 1]
Epdt <- read.table(paste0(prefix, "Alldata.txt"), sep = "\t", header = TRUE)
Epdt$S_ID <- paste0("Athr", 1:nrow(Epdt))


res <- read.csv('Stable17_Herbaceous&Woody.csv')
print(res)
colnames(res) <- gsub("\\.", "_", trimws(colnames(res)))


res$Strat <- as.character(as.character(res$Strat))
res$Category <- as.character(res$Category)



ylabs <- c("Predator performance", "Parasitoid performance", "NEs performance",
           "Herbivore performance", "Crop performance")


ylabs_short <- c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop")

cat_colors <- c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73")
names(cat_colors) <- c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop")


Strat <- c(1, 2)

Strat_labels <- c("Herbaceous_plant", "Woody_plant")


safe_condition <- function(pval, dat_sub) {
  if (is.null(pval)) return(FALSE)
  if (length(pval) != 1) return(FALSE)
  if (is.na(pval) || !is.finite(pval)) return(FALSE)
  if (pval >= 0.05) return(FALSE)
  if (length(unique(na.omit(dat_sub$added_plant_species))) <= 1) return(FALSE)
  return(TRUE)
}

k <- 1
all_plots <- list()
plot_labels <- letters
label_index <- 1



for (j in seq_along(Strat)) {
  for (i in seq_along(ylabs)) {
    
   
    dat_path <- paste0(prefix, total_names[i], ".txt")
    if (!file.exists(dat_path)) {
      cat("file does not exsist:", dat_path, "\n")
      next
    }
    dat <- read.table(dat_path, sep = "\t", header = TRUE)
    dat$Epdt <- Epdt[match(dat$S_ID, Epdt$S_ID), 2]
    dat <- dat[!is.na(dat$vi), ]
    
    cat("\ntreatment:", Strat_labels[j], "-", ylabs[i], "\n")
    cat("number:", nrow(dat), "\n")
    

    strat_cols <- c("CP_Ct", "CP_C3", "Mn__", "NE_c", "NCP_P", "Std_", "Zone")
    available_cols <- strat_cols[strat_cols %in% colnames(dat)]
    
    if (length(available_cols) == 0) {
      cat("warning: no strat\n")
      next
    }
    
    has_strat <- FALSE
    for (col in available_cols) {
      if (Strat_labels[j] %in% unique(dat[[col]])) {
        has_strat <- TRUE
        cat("in column", col, "find the data", Strat_labels[j], "\n")
        break
      }
    }
    
    if (!has_strat) {
      cat("warning: no data", Strat_labels[j], "\n")
      next
    }
    
    condition <- paste(available_cols, "=='", Strat_labels[j], "'", sep = "", collapse = " | ")
    dat_sub <- subset(dat, eval(parse(text = condition)))
    
    cat("after screening:", nrow(dat_sub), "\n")
    if (nrow(dat_sub) == 0) {
      cat("warning: after screening 0\n")
      k <- k + 1
      next
    }
    
    # --- added_plant_species ---
    if ("NCP_N" %in% colnames(dat_sub)) {
      dat_sub$NCP_N <- suppressWarnings(as.numeric(as.character(dat_sub$NCP_N)))
      dat_sub$added_plant_species <- ifelse(!is.na(dat_sub$NCP_N),
                                               log2(dat_sub$NCP_N), NA)
    } else {
      dat_sub$added_plant_species <- NA
    }
    
    valid_rows <- dat_sub[complete.cases(dat_sub[, c("added_plant_species", "yi", "vi")]), ]
    if (nrow(valid_rows) == 0) { 
      cat("warning: after screening 0\n")
      k <- k + 1
      next 
    }
    dat_sub <- valid_rows
    
    # 点大小
    wi <- 1 / sqrt(dat_sub$vi)
    size <- if (length(wi) == 1 || max(wi) == min(wi)) rep(3, length(wi)) else 1.5 + 5 * (wi - min(wi)) / (max(wi) - min(wi))
    dat_sub$size <- size
    
    
    res_row <- subset(res, Strat == Strat[j] & Category == ylabs_short[i])
    cat("from res the strat results:", Strat_labels[j], "-", ylabs_short[i], "\n")
    cat("number of rows:", nrow(res_row), "\n")
    
    if (!is.null(res_row) && nrow(res_row) == 1) {
      eq <- as.character(res_row$Regression_equation)
      pval <- as.numeric(as.character(res_row$P_value))
      obs <- res_row$n
      stud <- res_row$n_study
      cat("results: obs =", obs, "stud =", stud, "\n")
    } else {
      eq <- NA
      pval <- NA
      obs <- nrow(dat_sub)
      stud <- length(unique(dat_sub$S_ID))
      cat("no matched data，calculation: obs =", obs, "stud =", stud, "\n")
    }
    
    obs_label <- paste0("(", obs, "/", stud, ")")
    
    dat_sub$Category_short <- factor(ylabs[i],
                                     levels = ylabs,
                                     labels = c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop"))
    
    all_plots[[length(all_plots) + 1]] <- list(
      data = dat_sub,
      eq = eq,
      pval = pval,
      obs_label = obs_label,
      label = plot_labels[label_index],
      ylab = ylabs[i],
      strat = Strat_labels[j]   # 
    )
    
    label_index <- label_index + 1
    k <- k + 1
    cat("-----------------------------\n")
  }
}

cat("\n produced", length(all_plots), "pics\n")


panel_w <- 8   
panel_h <- 10  
nrow <- 3     
ncol <- 5     
if (length(all_plots) > 0) {
  png("Herbaceous&woody_res.png", units = "cm", width = 40, height = 30, res = 300)
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
          side = 1, line = 5, cex = 1, family = "Arial",col="black") 
    
    
    if (safe_condition(p$pval, dat_sub)) {
      if (length(unique(dat_sub$S_ID)) == 1) {
        fit <- lm(yi ~ added_plant_species, data = dat_sub, weights = 1 / vi)
      } else {
        fit <- tryCatch(
          lme(yi ~ added_plant_species, data = dat_sub,
              random = ~1 | S_ID, weights = varFixed(~vi),
              control = lmeControl(sigma = 1)),
          error = function(e) { return(NULL) }
        )
      }
      
      if (!is.null(fit)) {
        pred <- tryCatch(
          Effect("added_plant_species", fit, xlevels = list(added_plant_species = seq(0, 6, 0.01))),
          error = function(e) { return(NULL) }
        )
        
        if (!is.null(pred)) {
          polygon(c(rev(pred$x[,1]), pred$x[,1]), c(rev(pred$upper[,1]), pred$lower[,1]),
                  col = alpha("grey", 0.5), border = NA)
          lines(pred$x[,1], pred$fit[,1], col = "black", lwd = 3)
          if (!is.na(p$eq)) text(5.5, 6.8, labels = p$eq, cex = 1.2, adj = 1)
          text(5.5, 6.0, labels = paste("P =", format(signif(p$pval, 4), scientific = TRUE)), cex = 1.2, adj = 1)
        }
      }
    }
    
    points(dat_sub$added_plant_species, dat_sub$yi,
           bg = cat_colors[as.character(dat_sub$Category_short)], col = "black", pch = 21, cex = dat_sub$size)
    
    text(0.4, 7, labels = p$label, cex = 1.5, font = 2)
    text(5.5, 5.1, labels = p$obs_label, cex = 1.1, adj = 1)
  }
  
  dev.off()
  cat("pics saved Herbaceous&woody_res.png\n")
}

