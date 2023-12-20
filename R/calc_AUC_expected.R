#' Calculate the expected AUC under the QR model
#'
#' This function calculates the expected AUC under the quantitative response
#' (QR) model for a mixed meal tolerance test (MMTT), following the approach
#' developed by Bundy et al. 2020 (PMID 32704564), with expansions for other
#' timepoints by Ylescupidez et al. 2023 (PMID 37940642). As input,
#' it takes a data frame with age and baseline C-peptide values for each
#' subject, and a timepoint for which to calculate the expected AUC.
#' NOTE: Values for AUC must be in the standardized units used in the QR model,
#' specifically in units of nmol per L per minute from a 120-minute MMTT
#' (not log-transformed). Values for age must be in years.
#' 
#' The equation used is:
#' "Cpep_expected = -0.191 + 0.812 * ln(Cpep_baseline + 1) + 0.00638 * Age"
#' @param input_data data frame containing the C-peptide AUC values from MMTTs and subject ages.
#' @param timepoint numeric, the number of months at which to calculate the expected C-peptide value. Current implementation only supports 12.
#' @param cpeptide_baseline_colname character or numeric, the column containing the C-peptide values. Defaults to "aucMean".
#' @param age_colname character or numeric, the column containing the age values. Defaults to "ageYears".
#' @import checkmate
#' @export
#' @return a vector containing the expected C-peptide values at the specified timepoint, in ln(nmol/L/min + 1).
#' @usage
#' calc_AUC_expected_QR(
#'   input_data, timepoint = 12, cpeptide_baseline_colname = "aucMean", age_colname = "ageYears")
calc_AUC_expected_QR <-
  function(input_data, timepoint = 12, cpeptide_baseline_colname = "auc_mean", age_colname = "age_years") {
    # check input
    assert_data_frame(input_data)
    assert_number(timepoint)
    assert(
      check_string(cpeptide_baseline_colname),
      check_string(age_colname),
      combine = "and")
    assert(
      ifelse(cpeptide_baseline_colname %in% colnames(input_data),
        TRUE,
        "cpeptide_baseline_colname not found in input_data"),
      ifelse(age_colname %in% colnames(input_data),
        TRUE,
        "age_colname not found in input_data"),
      combine = "and")
    assert_numeric(input_data[[cpeptide_baseline_colname]])
    assert_numeric(input_data[[age_colname]])
    
    # calculate expected C-peptide values
    if (timepoint == 12) {
      cpep_expected <-
        -0.191 +
        0.812 * log(input_data[[cpeptide_baseline_colname]] + 1) +
        0.00638 * input_data[[age_colname]]
    } else {
      stop("Only timepoint = 12 is currently supported.")
    }

    return(cpep_expected)
  }
