rm=list(ls)
library(metafor)
library(Matrix)
library(matrixcalc)
library(lme4)
library(orchaRd)
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(lqmm)

prefix <- "/Users/yusha/"
data_names  <- read.table(paste0(prefix, "data_names.txt"),  sep = "\t", header = TRUE)[,1]
total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[,1]
n <- c(total_names, data_names)
View(data_names)
res <- NULL
recorded_names <- c()
p <- list()
q <- list()
k <- 1
l <- 1

for (i in n) {
  cat("in process：", i, "\n")
  dat <- read.table(paste0(prefix, i, ".txt"), sep = "\t", header = TRUE)
  dat <- dat[!is.na(dat$yi), ]
  
  
  if (i == "Total_response_of_herbivore") {
    res <- rbind(res, c(nrow(dat), length(unique(dat$S_ID)), rep(NA, 8)))
    recorded_names <- c(recorded_names, i)
    next
  }
  
  fit_success <- FALSE
  I2 <- NA
  dat.fit <- NULL
  
 
  try_v <- try(pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), dat)$reTrms$Zt)) - 1, 0), silent = TRUE)
  if (!inherits(try_v, "try-error")) {
    V <- tcrossprod(sqrt(dat$CONTROL_SD)) / sqrt(outer(dat$CONTROL_R * dat$CONTROL_M^2, dat$CONTROL_R * dat$CONTROL_M^2)) +
      diag((dat$TREATMENT_SD) / (dat$TREATMENT_R * dat$TREATMENT_M^2))
    V_ind <- try_v
    V <- V * V_ind
    
    
    dat.fit <- try(rma.mv(yi, V, random = list(~1 | S_ID,~1 | E_ID), data = dat, method = "ML", test = "t",
                          control = list(optimizer = "BFGS")), silent = TRUE)
    if (!inherits(dat.fit, "try-error")) {
      fit_success <- TRUE
      I2 <- tryCatch({
        100 * dat.fit$sigma2 / (dat.fit$sigma2 + mean(diag(dat.fit$V)))
      }, error = function(e) NA)
    }
  }
  
  
  if (!fit_success) {
    if (length(unique(dat$S_ID)) == 1) {
      dat.fit <- try(rma(yi, vi, data = dat, method = "ML", test = "t"), silent = TRUE)
      if (!inherits(dat.fit, "try-error")) {
        fit_success <- TRUE
        I2 <- tryCatch({
          100 * dat.fit$tau2 / (dat.fit$tau2 + mean(dat$vi))
        }, error = function(e) NA)
      }
    } else {
      dat.fit <- try(rma.mv(yi, vi, random = list(~1 | S_ID,~1 | E_ID), data = dat,
                            method = "ML", control = list(optimizer = "nlminb")), silent = TRUE)
      if (!inherits(dat.fit, "try-error")) {
        fit_success <- TRUE
        I2 <- tryCatch({
          100 * dat.fit$sigma2 / (dat.fit$sigma2 + mean(dat$vi))
        }, error = function(e) NA)
      }
    }
  }
  
 
  if (fit_success) {
    s <- tryCatch(coef(summary(dat.fit)), error = function(e) NULL)
    if (!is.null(s)) {
      res <- rbind(res, c(dat.fit$k,
                          length(unique(dat$S_ID)),
                          s["intrcpt", "estimate"],
                          s["intrcpt", "tval"],
                          s["intrcpt", "pval"],
                          nrow(dat) - 1,
                          s["intrcpt", "ci.lb"],
                          s["intrcpt", "ci.ub"],
                          I2))
    } else {
      res <- rbind(res, c(nrow(dat), length(unique(dat$S_ID)), rep(NA, 8)))
    }
  } else {
    res <- rbind(res, c(nrow(dat), length(unique(dat$S_ID)), rep(NA, 8)))
  }
  
  recorded_names <- c(recorded_names, i)
  

  if (!inherits(dat.fit, "try-error")) {
    plot_obj <- orchard_plot(dat.fit, group = "S_ID", xlab = "Standardised mean difference",
                             transfm = "none", twig.size = 2, trunk.size = 2) +
      theme(axis.title.y = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            axis.title.x = element_text(size = 20),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            plot.tag = element_text(size = 30))
    
    if (i %in% total_names) {
      plot_obj <- plot_obj +
        labs(tag = letters[k]) +
        scale_x_discrete("", labels = capitalize(gsub("Total response of ", "", gsub("_", " ", i))))
      p[[k]] <- plot_obj
      k <- k + 1
    } else {
      plot_obj <- plot_obj +
        labs(tag = letters[l]) +
        scale_x_discrete("", labels = capitalize(gsub(" response", "", gsub("_", " ", i))))
      q[[l]] <- plot_obj
      l <- l + 1
    }
  }
}

library(gridExtra)
library(ggplot2)

# labels a, b, c...
labels <- letters

save_plots_dynamic_custom <- function(plot_list, prefix_name,
                                      fixed_width = NULL,
                                      dynamic_width_logic = NULL,
                                      height = 30, dpi = 300,
                                      max_per_page = 8, max_ncol = 4) {
  num_plots <- length(plot_list)
  label_index <- 1
  
  
  if (!is.null(fixed_width)) {
    width_cm <- fixed_width
  } else if (!is.null(dynamic_width_logic)) {
    width_cm <- eval(dynamic_width_logic)
  } else {
    width_cm <- 80  # default fallback
  }
  
  for (i in seq(1, num_plots, by = max_per_page)) {
    plots_this_page <- plot_list[i:min(i + max_per_page - 1, num_plots)]
    n_plots <- length(plots_this_page)
    
    
    ncol_layout <- min(max_ncol, n_plots)
    nrow_layout <- ceiling(n_plots / ncol_layout)
    
   
    for (j in seq_along(plots_this_page)) {
      plots_this_page[[j]] <- plots_this_page[[j]] + labs(tag = labels[label_index])
      label_index <- label_index + 1
    }
    
  
    ggsave(
      filename = paste0(prefix_name, "_page", ceiling(i / max_per_page), ".png"),
      plot = do.call(grid.arrange, c(plots_this_page, nrow = nrow_layout, ncol = ncol_layout)),
      dpi = dpi,
      width = width_cm, height = height, units = "cm"
    )
  }
}


save_plots_dynamic_custom(p, "SuppFig9")  # 
save_plots_dynamic_custom(q, "SuppFig10")  #



# 输出汇总结果
res_df <- data.frame(Name = recorded_names, res)
colnames(res_df) <- c("Response", "n_obs", "n_groups", "estimate", "tval", "pval", "df", "ci.lb", "ci.ub", "I2")

View(res_df)
xlsx::write.xlsx(res, paste(path,"SuppTab3.xlsx",sep=""), row.names = F)

