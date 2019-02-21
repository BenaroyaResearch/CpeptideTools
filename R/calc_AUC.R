#' Calculate the trapezoidal area under the curve for a stimulated C-peptide assay
#'
#' This function calculates the area under the curve for a mixed-meal tolerance test (MMTT), oral glucose tolerance
#' test (OGTT), or similar, using the trapezoidal rule. As input, it takes a data frame of C-peptide or insulin
#' values for each timepoint, with an identifier column. This can be a long data frame with timepoint indicated in
#' a separete column, or a wide data frame with separate columns for each timepoint. Options include 2-hour vs.
#' 4-hour AUC, dropping visits with missing timepoints vs. using the average of the two surrounding timepoints,
#' and whether to use the 0 timepoint or the average of the pre-0 and 0 timepoints as the initial value.
#' @param input_data data frame containing the C-peptide data from MMTTs or stimulated insulin data from OGTTs. Should contain a unique visit/test identifer and the C-peptide level at each timepoint. If there are issues with values below the detection threshold, those must be resolved before running this function (i.e. values such as "<0.01" are not allowed)
#' @param input_format character, the format of the C-peptide data. Defaults to "wide", which has each timepoint in a separate column, with a single row for each visit. Set to "long" if the data frame has a separate row for each timepoint.
#' @param pep_colname character or numeric. Usage depends on the value of \code{input_format}. If \code{input_format} is set to "wide", this should be the prefix of the column identifiers for the timepoint. E.g. if the different timepoints are in columns named "pep0", "pep15", "pep30", then this argument should be set to "pep". If \code{input_format} is set to "long", this argument should be the name or number of the column containing the C-peptide information. Defaults to "pep".
#' @param timepoint_colname character or numeric, the column containing the timepoint values. Used only if \code{input_format} is set to "long".
#' @param initial_timepoint character, how to handle the initial timepoint. Defaults to "zero", which uses the 0-minute timepoint. Alternative option is "average", which uses the average of the pre-0 timepoint and the 0 timepoint.
#' @param treat_missing character, how to handle visits with missing timepoints. Defaults to "drop", which causes any visits missing a timepoint in the AUC window to be dropped. Alternative option is "weighted_average", which causes the AUC to be calculated as if that timepoint were not taken (i.e. the curve is drawn with a straight line between the preceding and following timepoints); this works only for non-consecutive missing timepoints that are not the first or last in the AUC calculation!
#' @param auc_timeframe character, either "2-hour" or "4-hour". Determines whether the AUC is calculated over 2 hours or 4.  If "4-hour" is selected, any 2-hour MMTTs will be dropped. Defaults to "2-hour".
#' @export
#' @return a vector containing the AUC values, with NA where it could not be calculated
#' @usage
#' calc_AUC(
#'   input_data, input_format="wide",
#'   pep_colname="pep", timepoint_colname="timepoint", 
#'   initial_timepoint="zero", treat_missing="drop",
#'   auc_timeframe="2-hour")
calc_AUC <-
  function(input_data, input_format="wide",
           pep_colname="pep", timepoint_colname="timepoint",
           initial_timepoint="zero", treat_missing="drop",
           auc_timeframe="2-hour") {
    
    input_format <- match.arg(input_format, choices = c("wide", "long"))
    treat_missing <- match.arg(treat_missing, choices = c("drop", "weighted_average"))
    
    # convert to wide for simplicity
    if (input_format=="long") {
      input_data <-
        input_data %>%
        tidyr::spread(key=!!timepoint_colname, value=!!pep_colname, sep="") %>%
        setNames(str_replace(names(.), paste0(timepoint_colname, "(?=([_.-]?[0-9]+$))"), pep_colname))
    }
    
    # create a data frame of timepoint columns and values
    timepoints <-
      data.frame(
        timepoint_col=
          colnames(input_data) %>%
          grep(pattern=paste0(pep_colname, "[_.-]?[0-9]+$"), x=., value=TRUE) %>%
          as.character(),
        stringsAsFactors=FALSE) %>%
      dplyr::mutate(
        timepoint_val=
          (timepoint_col %>%
             str_extract("[_.-]?[0-9]+$") %>%
             str_replace("[_.-]", "-") %>%
             as.numeric())) %>%
      dplyr::arrange(timepoint_val)
    
    # if 2-hour MMTT, drop all timepoints after 120
    if (auc_timeframe == "2-hour")
      timepoints <-
      timepoints %>%
      dplyr::filter(timepoint_val <= 120)
    
    # handle baseline value
    if (initial_timepoint == "average") {
      # average 0 and preceding timepoint
      input_data[, timepoints$timepoint_col[timepoints$timepoint_val==0]] <-
        apply(
          input_data[
            , timepoints$timepoint_col[
              c(which(timepoints$timepoint_val==0)-1, which(timepoints$timepoint_val==0))]],
          MARGIN = 1, FUN = mean)
    } else if (initial_timepoint != "zero") {
      stop("Input value for 'initial_timepoint' not recognized. Should be 'zero' or 'average'.")
    }
    
    # remove pre-0 timepoints
    timepoints <-
      timepoints %>%
      dplyr::filter(timepoint_val >= 0)
    
    # account for missing values (fill in if treat_missing is set to "weighted_average")
    if (treat_missing=="weighted_average") {
      for (timepoint.tmp in 2:(nrow(timepoints)-1)) {
        rows.tmp <- which(is.na(input_data[,timepoints$timepoint_col[timepoint.tmp]]))
        if (length(rows.tmp) > 0)
          input_data[rows.tmp, timepoints$timepoint_col[timepoint.tmp]] <-
            apply(
              input_data[rows.tmp, timepoints$timepoint_col[timepoint.tmp+c(-1,1)]],
              MARGIN=1, FUN=weighted.mean,
              w=abs(timepoints$timepoint_val[timepoint.tmp+c(-1,1)] - timepoints$timepoint_val[timepoint.tmp]))
      }
    }
    
    # calculate AUC as the matrix product of the time window averages and the time window lengths
    auc_values <-
      (as.matrix(
        input_data[,timepoints$timepoint_col[-1]] +
          input_data[,timepoints$timepoint_col[-nrow(timepoints)]]) / 2) %*%
      (timepoints$timepoint_val[-1] - timepoints$timepoint_val[-nrow(timepoints)])
      
    as.vector(auc_values)
  }
