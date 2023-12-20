#' Convert C-peptide values between standard units
#'
#' This function converts C-peptide values between different standard units. As
#' input, it takes the values (as a vector), the input units, and the desired
#' output units. These can be individual timepoints or AUC or mean AUC from an
#' OGTT or MMTT - the function is ambivalent.
#' @param values numeric vector C-peptide values. If there are issues with values below the detection threshold, those must be resolved before running this function (i.e. values such as "<0.01" are not allowed)
#' @param input_units character, the input units of the C-peptide values. Options include "ng_per_ml", "nmol_per_l", or less commonly "nmol_per_ml" or "ng_per_l". Each of these can be specified with and without "_per_min" appended. If "_per_min" is included, input values are assumed to be means across the time period of the MMTT or OGTT, and the number of minutes must be specified in \code{test_minutes}.
#' @param output_units character, the desired output units of the C-peptide values. Options are the same as for \code{input_units}. If "_per_min" is included, output values will be means across the time period of the MMTT or OGTT, and the number of minutes must be specified in \code{test_minutes}.
#' @param test_minutes numeric, the number of minutes for the OGTT or MMTT. Ignored if "_per_min" is not included as part of \code{input_units} or \code{output_units}. If both input and output units use "_per_min", the number of minutes must be the same (you cannot pass two different values for this argument).
#' @import checkmate
#' @import stringr
#' @export
#' @return a vector containing the C-peptide values, with NA where it could not be calculated
#' @usage
#' convert_Cpeptide_units(
#'   values, input_units, output_units, test_minutes = NULL)
convert_Cpeptide_units <-
  function(values, input_units, output_units, test_minutes = NULL) {
    # set conversion values, with standard units of ng/ml
    quantity <- c("ng" = 1, "nmol" = 1 / 3020.29)
    volume <- c("ml" = 1, "l" = 1000)
    
    # check validity of inputs
    assert_numeric(values) # verify that values are numeric
    assert_string(input_units) # verify that input_units is a string
    assert_string(output_units) # verify that output_units is a string
    # verify that test_minutes is a number if "_per_min" is included in input_units or output_units
    assert(
      ifelse (!str_detect(input_units, "_per_min") & !str_detect(output_units, "_per_min"), TRUE,
              "Argument 'test_minutes' must be a number if '_per_min' is included in 'input_units' or 'output_units'."),
      check_number(test_minutes),
      combine = "or")
    
    # extract units and per_min status from input_units and output_units
    input_units <- str_to_lower(input_units)
    input_units_quantity <- str_extract(input_units, "^(ng|nmol)")
    input_units_volume <- str_extract(input_units, "(?<=_per_)(ml|l)")
    input_units_per_min <- str_detect(input_units, "_per_min")
    output_units <- str_to_lower(output_units)
    output_units_quantity <- str_extract(output_units, "^(ng|nmol)")
    output_units_volume <- str_extract(output_units, "(?<=_per_)(ml|l)")
    output_units_per_min <- str_detect(output_units, "_per_min")

    # verify that input_units and output_units are valid
    assert(
      check_choice(input_units_quantity, c("ng", "nmol")),
      check_choice(input_units_volume, c("ml", "l")),
      check_logical(input_units_per_min, len = 1),
      check_choice(output_units_quantity, c("ng", "nmol")),
      check_choice(output_units_volume, c("ml", "l")),
      check_logical(output_units_per_min, len = 1)
      )

    # convert values
    values <-
      values *
      (quantity[output_units_quantity] * volume[output_units_volume] / (ifelse(output_units_per_min, test_minutes, 1))) /
      (quantity[input_units_quantity] * volume[input_units_volume] / (ifelse(input_units_per_min, test_minutes, 1)))
    
    return(values)
  }
