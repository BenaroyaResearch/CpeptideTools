#' Calculate quantitative response (QR) values for a given baseline and endpoint
#'
#' This function computes quantitative response (QR) estimates for each subject
#' by comparing their observed log-transformed AUC response to a fitted placebo trajectory.
#' A positive QR value indicates the subject's response exceeded what would be expected
#' under placebo, suggesting a treatment effect. The placebo trajectory is estimated from 
#' a pre-fitted linear model stored in the package's internal \code{base_QR_models_minimal} 
#' object, selected by the combination of baseline visit and response timepoint.
#'
#' @param data_to_add_QR_values_to data frame containing subject-level data, including AUC
#'   values at baseline and the response timepoint, and age at baseline. The QR columns will
#'   be appended to this data frame and returned.
#' @param current_response_time character, the label of the response/endpoint timepoint used
#'   to select the correct pre-fitted model (e.g. \code{"12 months"}).
#' @param baseline_date character, the label of the baseline visit used to select the correct
#'   pre-fitted model (e.g. \code{"randomization"}).
#' @param mean_auc_baseline_col character, the name of the column in \code{data_to_add_QR_values_to}
#'   containing the mean AUC value at baseline.
#' @param mean_auc_response_col character, the name of the column in \code{data_to_add_QR_values_to}
#'   containing the mean AUC value at the response timepoint.
#' @param age_at_baseline_col character, the name of the column in \code{data_to_add_QR_values_to}
#'   containing the subject's age at baseline.
#' @import dplyr
#' @import magrittr
#' @export
#' @return the input data frame with four additional columns appended: \code{log_mean_AUC_baseline}
#'   (log-transformed baseline AUC), \code{log_mean_AUC_response} (log-transformed response AUC),
#'   \code{lm_placebo_estimates} (fitted placebo trajectory from the pre-fitted model), and
#'   \code{QR} (the quantitative response estimate, i.e. the difference between the subject's
#'   observed log-transformed AUC response and the fitted placebo estimate; positive values
#'   indicate a response exceeding placebo expectation).
#' @usage
#' calc_QR_for_baseline_and_endpoint(
#'   data_to_add_QR_values_to,
#'   current_response_time,
#'   baseline_date,
#'   mean_auc_baseline_col,
#'   mean_auc_response_col,
#'   age_at_baseline_col)
#' @export
calc_QR_for_baseline_and_endpoint <- function(data_to_add_QR_values_to,
                                              current_response_time,
                                              baseline_date,
                                              mean_auc_baseline_col,
                                              mean_auc_response_col,
                                              age_at_baseline_col) {
  
  model_lookup_key <- paste0("baseline: ", baseline_date, " endpoint: ", current_response_time)
  
  if (!model_lookup_key %in% names(base_QR_models_minimal)) {
    stop(paste0(
      "No model found for '", model_lookup_key, "'.\n",
      "Available models: ", paste(names(base_QR_models_minimal), collapse = "; ")
    ))
  }
  
  selected_model <- base_QR_models_minimal[[model_lookup_key]]
  coefs          <- selected_model$coef
  
  # Build design matrix and compute fitted placebo values via matrix multiplication
  design_matrix       <- model.matrix(
    ~ log_mean_AUC_baseline + Age_At_Screening,
    data = data_to_add_QR_values_to %>%
      dplyr::transmute(
        log_mean_AUC_baseline = log(.data[[mean_auc_baseline_col]] + 1),
        log_mean_AUC_response = log(.data[[mean_auc_response_col]] + 1),
        Age_At_Screening      = .data[[age_at_baseline_col]]
      )
  )
  fitted_placebo_vals <- as.numeric(design_matrix %*% coefs)
  
  data_to_add_QR_values_to %>%
    dplyr::mutate(
      log_mean_AUC_baseline = log(.data[[mean_auc_baseline_col]] + 1),
      log_mean_AUC_response = log(.data[[mean_auc_response_col]] + 1),
      lm_placebo_estimates  = fitted_placebo_vals,
      QR                    = log_mean_AUC_response - lm_placebo_estimates
    )
}

