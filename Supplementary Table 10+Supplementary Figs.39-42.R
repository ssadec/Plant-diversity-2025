rm(list = ls())
library(metafor)
library(Matrix)      # nearPD
library(ggplot2)
library(gridExtra)
library(orchaRd)

prefix <- "/Users/yusha/"

data_names  <- read.table(paste0(prefix, "data_names.txt"), sep = "\t", header = TRUE)[,1]
total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[,1]

Zone<- c("Temperate","Tropic")
k <- 39
datasets <- c(total_names, data_names)
non_pd_list <- c()

# -------------------------------
# -------------------------------
save_plots_dynamic_custom <- function(
    plot_list, prefix_name,
    per_plot_width = 20,   
    per_plot_height = 15, 
    dpi = 300,
    max_per_page = 8, max_ncol = 4
) {
  num_plots <- length(plot_list)
  label_index <- 1
  
  for (i in seq(1, num_plots, by = max_per_page)) {
    plots_this_page <- plot_list[i:min(i + max_per_page - 1, num_plots)]
    n_plots <- length(plots_this_page)
    
    
    for (j in seq_along(plots_this_page)) {
      plots_this_page[[j]] <- plots_this_page[[j]] +
        labs(tag = letters[(label_index - 1) %% 26 + 1])
      label_index <- label_index + 1
    }
    
    if (n_plots == 1) {
 
      g <- plots_this_page[[1]]
      width_cm  <- per_plot_width
      height_cm <- per_plot_height
    } else {

      ncol_layout <- min(max_ncol, n_plots)
      nrow_layout <- ceiling(n_plots / ncol_layout)
      g <- do.call(gridExtra::arrangeGrob,
                   c(plots_this_page, nrow = nrow_layout, ncol = ncol_layout))
      width_cm  <- per_plot_width * ncol_layout
      height_cm <- per_plot_height * nrow_layout
    }
    
    ggsave(
      filename = paste0(prefix_name, "_page", ceiling(i / max_per_page), ".png"),
      plot = g,
      dpi = dpi,
      width = width_cm,
      height = height_cm,
      units = "cm",
      limitsize = FALSE
    )
  }
}


# -------------------------------
# -------------------------------
res <- data.frame(
  Zone = character(),
  Category = character(),
  Num_obs = numeric(),
  Num_studies = numeric(),
  Effect_size = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  df = numeric(),
  CI_lb = numeric(),
  CI_ub = numeric(),
  I2 = numeric(),
  stringsAsFactors = FALSE
)

# -------------------------------
# -------------------------------
for (i in Zone) {
  l <- 1
  m <- 1
  p <- list()
  q <- list()
  
  for (j in datasets) {

    filepath <- paste0(prefix, j, ".txt")
    if (!file.exists(filepath)) {
      message("File not found: ", filepath, " -> skip")
      res <- rbind(res, data.frame(
        Zone = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    dat <- read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    # ensure column name matching
    if (!("Zone" %in% names(dat))) {
      warning("No Zone column in ", j, " -> skip")
      res <- rbind(res, data.frame(
        Zone = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    dat <- subset(dat, dat$Zone == i)
    
    if (nrow(dat) == 0) {
  
      res <- rbind(res, data.frame(
        Zone = i, Category = j,
        Num_obs = 0, Num_studies = 0,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    

    if ("CONTROL_SD" %in% names(dat)) dat$CONTROL_SD[dat$CONTROL_SD < 0.01] <- 0.01
    if ("TREATMENT_SD" %in% names(dat)) dat$TREATMENT_SD[dat$TREATMENT_SD < 0.01] <- 0.01
    

    V_tmp <- try({
      tcrossprod(sqrt(dat$CONTROL_SD)) /
        sqrt(outer(dat$CONTROL_R * dat$CONTROL_M^2, dat$CONTROL_R * dat$CONTROL_M^2)) +
        diag(dat$TREATMENT_SD / (dat$TREATMENT_R * dat$TREATMENT_M^2))
    }, silent = TRUE)
    if (inherits(V_tmp, "try-error")) {
      message("Error computing V_tmp for ", j)
      res <- rbind(res, data.frame(
        Zone = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
  
    dat.fit <- NULL
    Zt <- try(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), dat)$reTrms$Zt, silent = TRUE)
    if (inherits(Zt, "try-error")) {
      message("Error computing Zt for ", j, "; attempting fixed-effect model...")
 
      dat.fit <- try(rma(yi, vi = diag(V_tmp), data = dat, method = "REML", test = "t"), silent = TRUE)
    } else {

      V_ind <- pmax(crossprod(Zt) - 1, 0)
      V_final <- as.matrix(V_tmp * V_ind)
      

      if (!is.positive.definite(V_final)) {
        message("Non-PD matrix detected for ", j, "; applying nearPD correction...")
        V_fixed <- tryCatch({
          as.matrix(nearPD(V_final)$mat)
        }, error = function(e) NULL)
        if (is.null(V_fixed)) {
          message("Failed to fix covariance matrix for ", j, "; skipping dataset.")
          non_pd_list <- c(non_pd_list, j)
          res <- rbind(res, data.frame(
            Zone = i, Category = j,
            Num_obs = NA, Num_studies = NA,
            Effect_size = NA, t_value = NA, P_value = NA,
            df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
            stringsAsFactors = FALSE
          ))
          next
        }
        V_final <- V_fixed
      }
      
  
      grp_levels <- length(unique(interaction(dat$S_ID, dat$CONTROL_M)))
      if (grp_levels >= nrow(dat)) {
        rand_form <- list(~1 | S_ID,~1 | E_ID)
      } else {
        rand_form <- list(~1 | S_ID / CONTROL_M)
      }
      
      dat.fit <- try(rma.mv(yi, V_final, random = rand_form, data = dat,
                            method = "ML", test = "t",
                            control = list(optimizer = "nlminb", rel.tol = 1e-6)), silent = TRUE)
    
      if (inherits(dat.fit, "try-error")) {
        message("Mixed model fit failed for ", j, "; attempting fixed-effect model...")
        dat.fit <- try(rma(yi, vi = diag(V_tmp), data = dat, method = "REML", test = "t"), silent = TRUE)
      }
    }
    

    if (inherits(dat.fit, "try-error") || is.null(dat.fit)) {
      message("Model fit failed for ", j, "; skipping.")
      res <- rbind(res, data.frame(
        Zone = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # ---------------- Extract statistics: effect, tval, p, CI, I2 ----------------

    beta_mat <- tryCatch(coef(summary(dat.fit)), error = function(e) NULL)
    
    # Estimate
    est <- tryCatch({
      if (!is.null(beta_mat)) {
        beta_mat[1, "estimate"]
      } else if (!is.null(dat.fit$b)) {
        as.numeric(dat.fit$b[1])
      } else NA
    }, error = function(e) NA)
    

    tval <- tryCatch({
      if (!is.null(beta_mat) && "tval" %in% colnames(beta_mat)) beta_mat[1, "tval"]
      else if (!is.null(beta_mat) && "zval" %in% colnames(beta_mat)) beta_mat[1, "zval"]
      else if (!is.null(dat.fit$tval)) dat.fit$tval[1]
      else NA
    }, error = function(e) NA)
    

    pval <- tryCatch({
      if (!is.null(beta_mat) && "pval" %in% colnames(beta_mat)) beta_mat[1, "pval"]
      else if (!is.null(beta_mat) && "p" %in% colnames(beta_mat)) beta_mat[1, "p"]
      else if (!is.null(dat.fit$pval)) dat.fit$pval[1]
      else NA
    }, error = function(e) NA)
    
    # CI
    ci_lb <- tryCatch({
      if (!is.null(beta_mat) && "ci.lb" %in% colnames(beta_mat)) beta_mat[1, "ci.lb"]
      else if (!is.null(dat.fit$ci.lb)) dat.fit$ci.lb[1]
      else NA
    }, error = function(e) NA)
    
    ci_ub <- tryCatch({
      if (!is.null(beta_mat) && "ci.ub" %in% colnames(beta_mat)) beta_mat[1, "ci.ub"]
      else if (!is.null(dat.fit$ci.ub)) dat.fit$ci.ub[1]
      else NA
    }, error = function(e) NA)
    

    I2 <- tryCatch({

      if (!is.null(dat.fit$sigma2) && length(dat.fit$sigma2) >= 1) {

        sigma2 <- if (length(dat.fit$sigma2) == 1) dat.fit$sigma2 else dat.fit$sigma2[1]
        100 * sigma2 / (sigma2 + mean(diag(dat.fit$V)))
      } else if (!is.null(dat.fit$tau2)) {
        
        tau2 <- dat.fit$tau2
      
        denom_mean_vi <- if ("vi" %in% colnames(dat)) mean(dat$vi, na.rm = TRUE) else mean(diag(V_tmp))
        100 * tau2 / (tau2 + denom_mean_vi)
      } else {
        NA
      }
    }, error = function(e) NA)
    

    df_val <- tryCatch({ nrow(dat) - 1 }, error = function(e) NA)
    

    num_obs <- tryCatch({ if (!is.null(dat.fit$k)) dat.fit$k else nrow(dat) }, error = function(e) nrow(dat))
    num_studies <- tryCatch(length(unique(dat$S_ID)), error = function(e) NA)
    

    res <- rbind(res, data.frame(
      Zone = i,
      Category = j,
      Num_obs = num_obs,
      Num_studies = num_studies,
      Effect_size = est,
      t_value = tval,
      P_value = pval,
      df = df_val,
      CI_lb = ci_lb,
      CI_ub = ci_ub,
      I2 = I2,
      stringsAsFactors = FALSE
    ))
    
    # ---------------- 绘图部分 ----------------
    if (j %in% total_names & nrow(dat) >= 1) {
      p[[l]] <- orchard_plot(dat.fit, group = "S_ID", xlab = "Standardised mean difference",
                             transfm = "none", twig.size = 2, trunk.size = 2) +
        theme(axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              plot.tag = element_text(size = 30)) +
        labs(tag = letters[l]) +
        scale_x_discrete("", labels = tools::toTitleCase(gsub("Total response of ", "", gsub("_", " ", j))))
      l <- l + 1
    } else if (nrow(dat) > 0) {
      q[[m]] <- orchard_plot(dat.fit, group = "S_ID", xlab = "Standardised mean difference",
                             transfm = "none", twig.size = 2, trunk.size = 2) +
        theme(axis.title.y = element_text(size = 20),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              axis.title.x = element_text(size = 20),
              axis.text.x = element_text(size = 20),
              axis.text.y = element_text(size = 20),
              plot.tag = element_text(size = 30)) +
        labs(tag = letters[m]) +
        scale_x_discrete("", labels = tools::toTitleCase(gsub(" response", "", gsub("_", " ", j))))
      m <- m + 1
    }
    
  } # end for j datasets
  

  if (length(p) > 0) {
    save_plots_dynamic_custom(
      plot_list = p,
      prefix_name = paste0("SuppFig", k),
      per_plot_width = 20,
      per_plot_height = 15,
      dpi = 300
    )
  }
  k <- k + 1
  

  if (length(q) > 0) {
    save_plots_dynamic_custom(
      plot_list = q,
      prefix_name = paste0("SuppFig", k),
      per_plot_width = 20,
      per_plot_height = 15,
      dpi = 300
    )
  }
  k <- k + 1
}

# -------------------------------

# -------------------------------
colnames(res) <- c("Climate Zone", "Category", "Number of observations",
                   "Number of studies", "Effect size", "t-value",
                   "P-value", "df", "Lower Bound", "Upper Bound", "I2")

if (length(non_pd_list)) {
  cat("Datasets failed to fix covariance:", paste(non_pd_list, collapse = ", "), "\n")
} else {
  cat("All matrices fixed successfully.\n")
}

View(res)
write.csv(res,'/Users/yusha/Stable10.csv')

