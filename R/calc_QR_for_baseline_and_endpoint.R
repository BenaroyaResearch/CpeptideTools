# Updating the calc_QR function to work for different baseline and endpoint...
# ... times.
calc_QR_for_baseline_and_endpoint <- function(data_to_add_QR_values_to,
                           current_response_time,
                           baseline_date,
                           mean_auc_baseline_col,
                           mean_auc_response_col,
                           age_at_baseline_col) {
  
  coef_label <- paste0("baseline: ", baseline_date, " endpoint: ", current_response_time)
  
  if (!coef_label %in% names(base_QR_models_minimal)) {
    stop(paste0(
      "No model found for '", coef_label, "'.\n",
      "Available models: ", paste(names(base_QR_models_minimal), collapse = "; ")
    ))
  }
  
  model_info <- base_QR_models_minimal[[coef_label]]
  coefs      <- model_info$coef
  
  newdata <- data_to_add_QR_values_to %>%
    dplyr::transmute(
      log_mean_AUC_baseline = log({{ mean_auc_baseline_col }} + 1),
      log_mean_AUC_response = log({{ mean_auc_response_col }} + 1),
      Age_At_Screening      = {{ age_at_baseline_col }}
    )
  
  X           <- model.matrix(~ log_mean_AUC_baseline + Age_At_Screening, data = newdata)
  fitted_vals <- as.numeric(X %*% coefs)
  
  data_to_add_QR_values_to %>%
    dplyr::mutate(
      log_mean_AUC_baseline = log({{ mean_auc_baseline_col }} + 1),
      log_mean_AUC_response = log({{ mean_auc_response_col }} + 1),
      lm_placebo_estimates  = fitted_vals,
      lm_estimated_te       = log_mean_AUC_response - lm_placebo_estimates
    )
}

