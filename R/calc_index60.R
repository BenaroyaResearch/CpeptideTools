#' Calculate the T1D Diagnostic Index60 value (Index60), a composite measure of pre-diagnosis T1D status
#'
#' This function calculates the T1D Diagnostic Index60 value (Index60) of
#' Sosenko et al. 2015 (PMID 25519451), a composite measure of pre-diagnosis
#' T1D status from an oral glucose tolerance test (OGTT) that has been found
#' to predict onset of T1D. Published work indicates that Index60 values above
#' 2.00 or 1.00 may predict onset of clinical T1D. The function takes an input
#' data frame with fasting C-peptide levels, 60-minute C-peptide levels from
#' OGTT, and 60-minute glucose levels from OGTT.
#' NOTE: Values for glucose and C-peptide must be in the standardized units used
#' in Index60, specifically ng/mL for C-peptide and mg/dL for glucose.
#' 
#' The equation used is:
#' "Index60 = 0.3695 * log (fasting C-peptide) + 0.0165 * (60-minute glucose) - 0.3644 * (60-minute C-peptide)"
#' @param input_data data frame containing the C-peptide and glucose values
#' @param cpeptide_fasting_colname character or numeric, the column containing the fasting C-peptide values in ng/mL. Defaults to "cpeptide_fasting".
#' @param cpeptide_60min_colname character or numeric, the column containing the OGTT 60-minute C-peptide values in ng/mL. Defaults to "cpeptide_60min".
#' @param glucose_60min_colname character or numeric, the column containing the OGTT 60-minute glucose values in mg/dL. Defaults to "glucose_60min".
#' @import checkmate
#' @export
#' @return a vector containing the Index60 values
#' @usage
#' calc_index60(
#'   input_data,
#'   cpeptide_fasting_colname = "pep0",
#'   cpeptide_60min_colname = "pep60",
#'   glucose_60min_colname = "glu60")
calc_index60 <-
  function(input_data,
           cpeptide_fasting_colname = "pep0",
           cpeptide_60min_colname = "pep60",
           glucose_60min_colname = "glu60") {
    # check input
    assert_data_frame(input_data)
    if (is.numeric(cpeptide_fasting_colname)) {
      cpeptide_fasting_colname <- colnames(input_data)[cpeptide_fasting_colname]
    }
    if (is.numeric(cpeptide_60min_colname)) {
      cpeptide_60min_colname <- colnames(input_data)[cpeptide_60min_colname]
    }
    if (is.numeric(glucose_60min_colname)) {
      glucose_60min_colname <- colnames(input_data)[glucose_60min_colname]
    }
    assert(
      check_string(cpeptide_fasting_colname),
      check_string(cpeptide_60min_colname),
      check_string(glucose_60min_colname),
      combine = "and")
    assert(
      ifelse(cpeptide_fasting_colname %in% colnames(input_data),
             TRUE,
             "cpeptide_fasting_colname not found in input_data"),
      ifelse(cpeptide_60min_colname %in% colnames(input_data),
             TRUE,
             "cpeptide_60min_colname not found in input_data"),
      ifelse(glucose_60min_colname %in% colnames(input_data),
             TRUE,
             "glucose_60min_colname not found in input_data"),
      combine = "and")
    assert_numeric(input_data[[cpeptide_fasting_colname]])
    assert_numeric(input_data[[cpeptide_60min_colname]])
    assert_numeric(input_data[[glucose_60min_colname]])
    
    # calculate Index60 values
    index60 <- 
      0.3695 * log(input_data[[cpeptide_fasting_colname]]) +
      0.0165 * input_data[[glucose_60min_colname]] -
      0.3644 * input_data[[cpeptide_60min_colname]]
    
    return(index60)
  }
