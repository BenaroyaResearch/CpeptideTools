#' Calculate the quantitative response (QR) of C-peptide AUC
#'
#' This function calculates the quantitative response (QR) of Bundy et al. 2020
#' (PMID 32704564) for a mixed meal tolerance test (MMTT), with expansions for
#' other timepoints by Ylescupidez et al. 2023 (PMID 37940642). As input, it
#' takes a data frame with age, baseline C-peptide AUC, and C-peptide AUC at
#' the specified timepoint (all for each subject), and the timepoint for which
#' to calculate the QR.
#' NOTE: Values for AUC must be in the standardized units used in the QR model,
#' specifically in units of nmol per L per minute from a 120-minute MMTT
#' (not log-transformed). Values for age must be in years.
#' 
#' The equation used is:
#' "QR = Cpep_actual - Cpep_expected = ln(Cpep_timepoint + 1) + 0.191 - 0.812 * ln(Cpep_baseline + 1) - 0.00638 * Age"
#' @param input_data data frame containing the C-peptide AUC values from MMTTs and subject ages.
#' @param timepoint numeric, the number of months at which to calculate the expected C-peptide value. Current implementation only supports 12.
#' @param cpeptide_baseline_colname character or numeric, the column containing the C-peptide values. Defaults to "auc_mean_baseline".
#' @param cpeptide_timepoint_colname character or numeric, the column containing the C-peptide values. Defaults to "auc_mean".
#' @param age_colname character or numeric, the column containing the age values. Defaults to "age_years".
#' @import checkmate
#' @export
#' @return a vector containing the QR values at the specified timepoint, in ln(nmol/L/min + 1).
#' @usage
#' calc_QR(
#'   input_data, timepoint = 12, cpeptide_baseline_colname = "auc_mean_baseline",
#'   cpeptide_timepoint_colname = "auc_mean", age_colname = "age_years")
calc_QR <-
  function(input_data, timepoint = 12,
           cpeptide_baseline_colname = "auc_mean_baseline",
           cpeptide_timepoint_colname = "auc_mean",
           age_colname = "age_years") {
    # check input
    assert_data_frame(input_data)
    assert_number(timepoint)
    assert(
      check_string(cpeptide_baseline_colname),
      check_string(cpeptide_timepoint_colname),
      check_string(age_colname),
      combine = "and")
    assert(
      ifelse(cpeptide_baseline_colname %in% colnames(input_data),
             TRUE,
             "cpeptide_baseline_colname not found in input_data"),
      ifelse(cpeptide_timepoint_colname %in% colnames(input_data),
             TRUE,
             "cpeptide_timepoint_colname not found in input_data"),
      ifelse(age_colname %in% colnames(input_data),
             TRUE,
             "age_colname not found in input_data"),
      combine = "and")
    assert_numeric(input_data[[cpeptide_baseline_colname]])
    assert_numeric(input_data[[cpeptide_timepoint_colname]])
    assert_numeric(input_data[[age_colname]])
    
    # calculate expected C-peptide values
    cpep_expected <- calc_AUC_expected_QR(
      input_data, timepoint = timepoint,
      cpeptide_baseline_colname = cpeptide_baseline_colname,
      age_colname = age_colname)

    # calculate QR from timepoint and expected C-peptide values
    cpep_actual <- log(input_data[[cpeptide_timepoint_colname]] + 1)
    qr <- cpep_actual - cpep_expected
    
    return(qr)
  }
