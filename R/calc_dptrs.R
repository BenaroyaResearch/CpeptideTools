#' Calculate the Diabetes Prevention Trial-Type 1 Risk Score (DPTRS), a composite measure of pre-diagnosis T1D status
#'
#' This function calculates the Diabetes Prevention Trial-Type 1 Risk Score of
#' Sosenko et al. 2015 (PMID 26077017) and Sosenko et al. 2008 (PMID  18000175),
#' a composite measure of pre-diagnosis T1D status from an oral glucose
#' tolerance test (OGTT) that has been found to predict onset of T1D. Published
#' work indicates that Index60 values above 7.0 or 9.0 may predict onset of
#' clinical T1D. The function takes an input data frame with fasting C-peptide
#' levels; glucose data from 30, 60, 90, 120 minutes, C-peptide data from 30,
#' 60, 90, 120 minutes; ; and age.
#' NOTE: Values for all input data must must be in the standardized units used
#' in DPTRS, specifically ng/mL for C-peptide, mg/dL for glucose, km/m^2 for
#' BMI, and years for age.
#'
#' The equation used is:
#' "DPTRS = (1.57 * log BMI) − (0.06 * age) + (0.81 * glucose sum from 30 to 120 min/100) − (0.85 * C-peptide sum from 30 to 120 min/10) + (0.48 * log fasting C-peptide)
#' @param input_data data frame containing the C-peptide, glucose values, BMI, and age data
#' @param cpeptide_fasting_colname character or numeric, the column containing the fasting C-peptide values in ng/mL. Defaults to "cpeptide_fasting".
#' @param cpeptide_sum_colname character or numeric, the column containing the sum of OGTT C-peptide values in ng/mL for 30, 60, 90, 120 minutes. Defaults to NULL, and EITHER this or \code{cpeptide_colname_prefix} must be provided, but not both.
#' @param cpeptide_colname_prefix character, the name prefix of the columns containing the OGTT C-peptide values in ng/mL for 30, 60, 90, 120 minutes. Defaults to NULL, and EITHER this or \code{cpeptide_sum_colname} must be provided, but not both.
#' @param glucose_sum_colname character or numeric, the column containing the sum of OGTT C-peptide values in mg/dL for 30, 60, 90, 120 minutes. Defaults to NULL, and EITHER this or \code{glucose_colname_prefix} must be provided, but not both.
#' @param glucose_colname_prefix character, the name prefix of the columns containing the OGTT C-peptide values in mg/dL for 30, 60, 90, 120 minutes. Defaults to NULL, and EITHER this or \code{glucose_sum_colname} must be provided, but not both.
#' @param bmi_colname, character or numeric, the column containing the BMI data. Defaults to "bmi"
#' @param age_colname, character or numeric, the column containing the age data. Defaults to "age"
#' @import checkmate
#' @export
#' @return a vector containing the DPTRS values
#' @usage
#' calc_dptrs(
#'   input_data,
#'   cpeptide_fasting_colname = "pep0",
#'   cpeptide_sum_colname = NULL,
#'   cpeptide_colname_prefix = NULL,
#'   glucose_sum_colname = NULL,
#'   glucose_colname_prefix = NULL,
#'   bmi_colname,
#'   age_colname)
calc_dptrs <-
  function(
    input_data,
    cpeptide_fasting_colname = "pep0",
    cpeptide_sum_colname = NULL,
    cpeptide_colname_prefix = NULL,
    glucose_sum_colname = NULL,
    glucose_colname_prefix = NULL,
    bmi_colname = "bmi",
    age_colname = "age"
  ) {
    # check input
    assert_data_frame(input_data)

    # check C-peptide fasting data
    if (is.numeric(cpeptide_fasting_colname)) {
      cpeptide_fasting_colname <- colnames(input_data)[cpeptide_fasting_colname]
    }
    assert(check_string(cpeptide_fasting_colname))
    assert_names(
      colnames(input_data),
      must.include = cpeptide_fasting_colname
    )
    assert_numeric(input_data[[cpeptide_fasting_colname]])

    cpeptide_fasting <- input_data[[cpeptide_fasting_colname]]

    # check and calculate C-peptide sum
    if (
      (is.null(cpeptide_sum_colname) & is.null(cpeptide_colname_prefix)) ||
        (!is.null(cpeptide_sum_colname) & !is.null(cpeptide_colname_prefix))
    ) {
      stop(
        "Either 'cpeptide_sum_colname' or 'cpeptide_colname_prefix' must be provided, but not both"
      )
    }
    if (!is.null(cpeptide_sum_colname)) {
      assert(
        ifelse(
          cpeptide_sum_colname %in% colnames(input_data),
          TRUE,
          "cpeptide_sum_colname not found in input_data"
        )
      )
      assert_numeric(input_data[[cpeptide_sum_colname]])
      cpeptide_sum <- input_data[[cpeptide_sum_colname]]
    } else if (!is.null(cpeptide_colname_prefix)) {
      cpeptide_timepoints <- c(30, 60, 90, 120)
      cpeptide_colnames <- paste0(cpeptide_colname_prefix, cpeptide_timepoints)
      if (!all(cpeptide_colnames %in% colnames(input_data))) {
        stop(
          "Input data must contain columns with 'cpeptide_colname_prefix' followed by '30', '60', '90', and '120'"
        )
      }
      if (!all(sapply(input_data[, cpeptide_colnames], is.numeric))) {
        stop(
          "Input data must be numeric for all columns with 'cpeptide_colname_prefix' followed by '30', '60', '90', and '120'"
        )
      }
      cpeptide_sum <- rowSums(input_data[, cpeptide_colnames])
    }

    # check and calculate glucose sum
    if (
      (is.null(glucose_sum_colname) & is.null(glucose_colname_prefix)) ||
        (!is.null(glucose_sum_colname) & !is.null(glucose_colname_prefix))
    ) {
      stop(
        "Either 'glucose_sum_colname' or 'glucose_colname_prefix' must be provided, but not both"
      )
    }
    if (!is.null(glucose_sum_colname)) {
      assert(
        ifelse(
          glucose_sum_colname %in% colnames(input_data),
          TRUE,
          "glucose_sum_colname not found in input_data"
        )
      )
      assert_numeric(input_data[[glucose_sum_colname]])
      glucose_sum <- input_data[[glucose_sum_colname]]
    } else if (!is.null(glucose_colname_prefix)) {
      glucose_timepoints <- c(30, 60, 90, 120)
      glucose_colnames <- paste0(glucose_colname_prefix, glucose_timepoints)
      if (!all(glucose_colnames %in% colnames(input_data))) {
        stop(
          "Input data must contain columns with 'glucose_colname_prefix' followed by '30', '60', '90', and '120'"
        )
      }
      if (!all(sapply(input_data[, glucose_colnames], is.numeric))) {
        stop(
          "Input data must be numeric for all columns with 'glucose_colname_prefix' followed by '30', '60', '90', and '120'"
        )
      }
      glucose_sum <- rowSums(input_data[, glucose_colnames])
    }

    # check BMI data
    if (is.numeric(bmi_colname)) {
      bmi_colname <- colnames(input_data)[bmi_colname]
    }
    assert(check_string(bmi_colname))
    assert(
      ifelse(
        bmi_colname %in% colnames(input_data),
        TRUE,
        "bmi_colname not found in input_data"
      )
    )
    assert_numeric(input_data[[bmi_colname]])
    bmi <- input_data[[bmi_colname]]

    # check age data
    if (is.numeric(age_colname)) {
      age_colname <- colnames(input_data)[age_colname]
    }
    assert(check_string(age_colname))
    assert(
      ifelse(
        age_colname %in% colnames(input_data),
        TRUE,
        "age_colname not found in input_data"
      )
    )
    assert_numeric(input_data[[age_colname]])
    age <- input_data[[age_colname]]

    # calculate DPTRS values
    dptrs <-
      1.57 *
      log(bmi) +
      -0.06 * age +
      0.81 * (glucose_sum / 100) +
      -0.85 * (cpeptide_sum / 10) +
      0.48 * log(cpeptide_fasting)

    return(dptrs)
  }
