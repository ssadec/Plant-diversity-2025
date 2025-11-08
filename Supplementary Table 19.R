rm(list = ls())

prefix <- "/Users/yusha/"
path   <- "/Users/yusha/"

Category <- c(
  "Total_response_of_predator_performance",
  "Total_response_of_parasitoid_performance",
  "Total_response_of_NEs_performance",
  "Total_response_of_herbivore_performance",
  "Total_response_of_crop_performance"
)

total_names <- read.table(paste0(prefix, "total_names.txt"),
                          sep = "\t", header = TRUE)[, 1]

res_all <- data.frame()

for (i in seq_along(total_names)) {
  cat("\n========== Processing:", Category[i], "==========\n")
  
  # ==== Step 1. Read data ====
  dat <- try(read.table(paste0(prefix, total_names[i], ".txt"),
                        sep = "\t", header = TRUE), silent = TRUE)
  if (inherits(dat, "try-error") || nrow(dat) < 5) {
    cat("⚠️ Data read failed or sample size too small, skipping:", Category[i], "\n")
    next
  }
  
  dat <- dat[!is.na(dat$yi), ]
  if (nrow(dat) < 5) next
  
  # ==== Step 2. Remove extreme effect sizes and variances ====
  # Remove extremely large yi or vi (e.g.,extreme quantiles)
  
  if (nrow(dat) < 5) {
    cat("⚠️ Sample size too small after cleaning, skipping:", Category[i], "\n")
    next
  }
  
  # ==== Step 3. Construct covariance matrix ====
  try_v <- try(
    pmax(as.matrix(crossprod(lFormula(yi ~ 1 + (1 | S_ID / CONTROL_M), dat)$reTrms$Zt)) - 1, 0),
    silent = TRUE
  )
  
  if (!inherits(try_v, "try-error")) {
    V <- tcrossprod(sqrt(dat$CONTROL_SD)) /
      sqrt(outer(dat$CONTROL_R * dat$CONTROL_M^2,
                 dat$CONTROL_R * dat$CONTROL_M^2)) +
      diag((dat$TREATMENT_SD) / (dat$TREATMENT_R * dat$TREATMENT_M^2))
    V <- V * try_v
  } else {
    cat("⚠️ try_v construction failed, using diagonal matrix:", Category[i], "\n")
    V_diag <- (dat$CONTROL_SD / (dat$CONTROL_R * dat$CONTROL_M^2)) +
      (dat$TREATMENT_SD / (dat$TREATMENT_R * dat$TREATMENT_M^2))
    V <- diag(V_diag)
  }
  
  # ==== Step 4. Fit bias-robust fixed effect model ====
  mod_MLFE <- try(
    rma.mv(yi = yi, V = V, random = list(~1 | S_ID), method = "ML", test = "t", dfs = "contain", data = dat),
    silent = TRUE
  )
  if (inherits(mod_MLFE, "try-error")) {
    cat("❌ Model fitting failed:", Category[i], "\n")
    next
  }
  
  # ==== Step 5. Cluster-robust variance estimation (CRVE) ====
  mod_MLFE_RVE <- try(
    robust(mod_MLFE, cluster = dat$S_ID, adjust = TRUE, clubSandwich = TRUE),
    silent = TRUE
  )
  if (inherits(mod_MLFE_RVE, "try-error")) {
    cat("⚠️ CRVE calculation failed:", Category[i], "\n")
    next
  }
  
  # ==== Step 6. Check results ====
  if (is.null(mod_MLFE_RVE$b) || length(mod_MLFE_RVE$b) == 0) {
    cat("⚠️ Model returned no estimates, skipping:", Category[i], "\n")
    next
  }
  
  # ==== Step 7. Save results ====
  if (!is.null(mod_MLFE_RVE$b) && length(mod_MLFE_RVE$b) > 0 && 
      !is.null(mod_MLFE_RVE$se) && length(mod_MLFE_RVE$se) > 0) {
    
    beta_val <- as.numeric(mod_MLFE_RVE$b[1])
    se_val <- as.numeric(mod_MLFE_RVE$se[1])
    
    # Manually calculate t value: t = beta / SE
    if (!is.na(se_val) && se_val != 0) {
      tval_calculated <- beta_val / se_val
    } else {
      tval_calculated <- NA
    }
    
    # Use model p-value if available, otherwise set NA
    if (!is.null(mod_MLFE_RVE$pval) && length(mod_MLFE_RVE$pval) > 0) {
      pval_used <- as.numeric(mod_MLFE_RVE$pval[1])
    } else {
      pval_used <- NA
    }
    
    tmp <- data.frame(
      Category = Category[i],
      k = mod_MLFE_RVE$k,
      n_cluster = length(unique(dat$S_ID)),
      beta = beta_val,
      SE = se_val,
      tval = tval_calculated,
      pval = pval_used,
      CIlb = as.numeric(mod_MLFE_RVE$ci.lb[1]),
      CIub = as.numeric(mod_MLFE_RVE$ci.ub[1])
    )
    
    res_all <- rbind(res_all, tmp)
    cat("✅ Completed:", Category[i], "\n")
    
  } else {
    cat("⚠️ Model returned no valid estimates or SE, skipping:", Category[i], "\n")
  }
}

# ==== Step 8. Output results ====
View(res_all)
write.csv(res_all, paste0(path, "TwoStep_Robust_Meta_Results.csv"), row.names = FALSE)
cat("\n✅ All results saved to: TwoStep_Robust_Meta_Results.csv\n")
