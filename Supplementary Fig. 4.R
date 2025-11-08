library(nlme)
library(effects)
library(scales)


prefix <- "/Users/yusha/"
prefix2 <- "/Users/yusha/Desktop/check/"




total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[, 1]
Epdt <- read.table(paste0(prefix, "Alldata.txt"))
Epdt$S_ID <- paste0("Athr", 1:nrow(Epdt))
res <- read.csv(paste0(prefix2, "Stable17_Herbivore_G&S.csv"), na.strings = c("NA", ""))
colnames(res) <- gsub("\\.", "_", trimws(colnames(res)))


ylabs <- c("Predator performance",
           "Parasitoid performance",
           "NEs performance",
           "Herbivore performance",
           "Crop performance")


cat_colors <- c("#BC3A24", "#EFB882", "#06798F", "#7A9A01", "#384D73")
names(cat_colors) <- c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop")


Strat <- c("Generalist", "Specialist")


all_plots <- list()
plot_labels <- letters
label_index <- 1
k <- 1

for (j in seq_along(Strat)) {
  for (i in seq_along(ylabs)) {
    dat <- read.table(paste0(prefix, total_names[i], ".txt"), sep = "\t", header = TRUE)
    dat$Epdt <- Epdt[match(dat$S_ID, Epdt$S_ID), 2]
    dat <- dat[!is.na(dat$vi), ]
    
    dat_sub <- subset(dat, Mn__ == Strat[j])
    
  
    valid_rows <- na.omit(dat_sub[, c("added_plant_species", "yi", "vi")])
    if (nrow(valid_rows) == 0) {
      k <- k + 1
      next
    }
    

    wi <- 1 / sqrt(valid_rows$vi)
    size <- if (length(wi) == 1 || max(wi) == min(wi)) {
      rep(3, length(wi))
    } else {
      1.5 + 5 * (wi - min(wi)) / (max(wi) - min(wi))
    }
    dat_sub$size <- size
    
    eq <- if (k <= nrow(res)) as.character(res$Regression_equation[k]) else NA
    pval <- if (k <= nrow(res)) as.numeric(as.character(res$P_value[k])) else NA
    obs <- if (k <= nrow(res)) res$Number_of_observations[k] else NA
    stud <- if (k <= nrow(res)) res$Number_of_studies[k] else NA
    obs_label <- paste0("(", obs, "/", stud, ")")
    
    dat_sub$Category <- factor(ylabs[i],
                               levels = ylabs,
                               labels = c("Predator", "Parasitoid", "NEs", "Herbivore", "Crop"))
    
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


subplot_labels <- letters[1:length(all_plots)]


png("Herbivore_G&S.png", units = "cm", width = 40, height = 30, res = 300)

panel_w <- 8   
panel_h <- 10  
nrow <- 3    
ncol <- 5     

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
  

  pval_ok <- !is.null(p$pval) && !is.na(p$pval) && is.finite(p$pval) && p$pval < 0.05
  richness_ok <- "added_plant_species" %in% colnames(dat_sub) &&
    sum(!is.na(dat_sub$added_plant_species)) > 1
  

  if (pval_ok && richness_ok) {
    if (length(unique(na.omit(dat_sub$S_ID))) == 1) {
      fit <- lm(yi ~ added_plant_species, data = dat_sub, weights = 1 / vi)
    } else {
      fit <- lme(yi ~ added_plant_species, data = dat_sub,
                 random = ~1 | S_ID, weights = varFixed(~vi),
                 control = lmeControl(sigma = 1))
    }
    
    pred <- Effect("added_plant_species", fit,
                   xlevels = list(added_plant_species = seq(0, 6, 0.01)))
    
    polygon(c(rev(pred$x[, 1]), pred$x[, 1]),
            c(rev(pred$upper[, 1]), pred$lower[, 1]),
            col = alpha("grey", 0.5), border = NA)
    lines(pred$x[, 1], pred$fit[, 1], col = "black", lwd = 3)
    
    if (!is.null(p$eq) && !is.na(p$eq)) {
      text(5.5, 6.8, labels = p$eq, cex = 1.2, adj = 1)
    }
    text(5.5, 6.0, labels = paste("P =", format(signif(p$pval, 4), scientific = TRUE)),
         cex = 1.2, adj = 1)
  }
  
  points(dat_sub$added_plant_species, dat_sub$yi,
         bg = cat_colors[dat_sub$Category], col = "black", pch = 21, cex = dat_sub$size)
  
  text(0.4, 7, labels = subplot_labels[i], cex = 1.5, font = 2)
  
 
  if (i == 10) {
    text(5.5, 5.1, labels = "(906/123)", cex = 1.1, adj = 1)
  } else {
    text(5.5, 5.1, labels = p$obs_label, cex = 1.1, adj = 1)
  }
}

dev.off()

