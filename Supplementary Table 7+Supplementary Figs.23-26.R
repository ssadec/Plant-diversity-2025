library(metafor)
library(Matrix)      # nearPD
library(ggplot2)
library(gridExtra)
library(orchaRd)

prefix <- "/Users/yusha/"

data_names  <- read.table(paste0(prefix, "data_names.txt"), sep = "\t", header = TRUE)[,1]
total_names <- read.table(paste0(prefix, "total_names.txt"), sep = "\t", header = TRUE)[,1]


NE_c <- c("generalist", "specialist")  
k <- 23
datasets <- c(total_names, data_names)
non_pd_list <- c()

# -------------------------------

# -------------------------------
save_plots_dynamic_custom <- function(
    plot_list, prefix_name,
    per_plot_width = 20, per_plot_height = 15,
    dpi = 300, max_per_page = 8, max_ncol = 4
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
      plot = g, dpi = dpi,
      width = width_cm, height = height_cm,
      units = "cm", limitsize = FALSE
    )
  }
}

# -------------------------------
# -------------------------------
res <- data.frame(
  NE_c = character(),
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
for (i in NE_c) {
  l <- 1
  m <- 1
  p <- list()
  q <- list()
  
  for (j in datasets) {
    filepath <- paste0(prefix, j, ".txt")
    if (!file.exists(filepath)) {
      message("File not found: ", filepath, " -> skip")
      res <- rbind(res, data.frame(
        NE_c = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    dat <- read.table(filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    
  
    if (!("NE_c" %in% names(dat))) {
      warning("No NE_c column in ", j, " -> skip")
      res <- rbind(res, data.frame(
        NE_c = i, Category = j,
        Num_obs = NA, Num_studies = NA,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    dat$NE_c <- tolower(trimws(dat$NE_c))
    
    dat <- subset(dat, dat$NE_c == i)
    if (nrow(dat) == 0) {
      res <- rbind(res, data.frame(
        NE_c = i, Category = j,
        Num_obs = 0, Num_studies = 0,
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    

    dat.fit <- tryCatch(
      rma.mv(yi, vi, random = (~ 1 | S_ID,~1 | E_ID), data = dat),
      error = function(e) NULL
    )
    
    if (is.null(dat.fit)) {
      message("Model failed for ", j)
      res <- rbind(res, data.frame(
        NE_c = i, Category = j,
        Num_obs = nrow(dat), Num_studies = length(unique(dat$S_ID)),
        Effect_size = NA, t_value = NA, P_value = NA,
        df = NA, CI_lb = NA, CI_ub = NA, I2 = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    

    est <- summary(dat.fit)
    
  
    tau2 <- sum(dat.fit$sigma2)
 
    vi_bar <- mean(dat$vi, na.rm = TRUE)
 
    I2_val <- tau2 / (tau2 + vi_bar)
    
    res <- rbind(res, data.frame(
      NE_c = i, Category = j,
      Num_obs = nrow(dat),
      Num_studies = length(unique(dat$S_ID)),
      Effect_size = est$b[1],
      t_value = est$zval,     
      P_value = est$pval,
      df = NA,                
      CI_lb = est$ci.lb,
      CI_ub = est$ci.ub,
      I2 = I2_val,
      stringsAsFactors = FALSE
    ))
    
    # ---------------- pics----------------
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
        scale_x_discrete("", labels = tools::toTitleCase(gsub(" response", "", gsub("_", " ", j))))
      m <- m + 1
    }
    
  } # end for j datasets
  

  if (length(p) > 0) {
    save_plots_dynamic_custom(p, paste0("SuppFig", k), per_plot_width = 20, per_plot_height = 15, dpi = 300)
  }
  k <- k + 1
  
  if (length(q) > 0) {
    save_plots_dynamic_custom(q, paste0("SuppFig", k), per_plot_width = 20, per_plot_height = 15, dpi = 300)
  }
  k <- k + 1
}

# -------------------------------

# -------------------------------
colnames(res) <- c("NEs_category", "Category", "Number of observations",
                   "Number of studies", "Effect size", "t-value",
                   "P-value", "df", "Lower Bound", "Upper Bound", "I2")
View(res)

write.csv(res,'/Users/yusha/Stable7.csv', row.names = FALSE)

