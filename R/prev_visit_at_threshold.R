#' Determine which visits for any subject occur after C-peptide AUC has reached the detection limit
#'
#' This function determines which, if any, visits for any subject occur after that subject's C-peptide AUC has
#' reached the detection limit. For AUC values, this means that the C-peptide levels for each timepoint are at the
#' detection limit; the AUC detection limit should be calculated as such. As input, it takes a data frame of
#' C-peptide values, with subject identifiers and a date or visit column on which to sort the visits. Note that
#' visits will be sorted blindly, so character visit names (e.g. Day 1, Week 1, Month 2) should be used with EXTREME
#' CAUTION!
#' @param cpeptide_auc_data data frame containing the C-peptide data. Should contain a unique subject identifer, a column for sorting visits, and the C-peptide AUC data. Each row should be a single visit.
#' @param threshold numeric, the threshold C-peptide AUC value. Values at or below this level will be treated as sub-threshold.
#' @param identifier_column character or numeric, the column containing the subject identifiers. By default it uses "subject".
#' @param sort_column character or numeric, the column on which each subject's visits should be sorted. This should be either numeric, date, or character that will be sorted consistently by dplyr::arrange(). If set to NULL, rows for each subject will be assumed to be in the desired order. By default it uses "visit".
#' @param auc_column character or numeric, the column containing the C-peptide AUC values. By default it uses "auc".
#' @export
#' @return a logical vector, with each element corresponding to a row of the input data frame; contains TRUE if the subject already had C-peptide AUC at or below the threshold at prior visits, and FALSE otherwise.
#' @usage \code{
#' prev_visit_at_threshold(
#'   cpeptide_auc_data, threshold,
#'   identifier_column="subject", sort_column="visit",
#'   auc_column="auc"
#'   )}
prev_visit_at_threshold <-
  function(cpeptide_auc_data, threshold,
           identifier_column="subject", sort_column="visit",
           auc_column="auc") {
    
    # ensure character column names (necessary for dplyr functions)
    if (is.numeric(identifier_column)) identifier_column <- colnames(cpeptide_auc_data)[identifier_column]
    if (is.numeric(sort_column)) sort_column <- colnames(cpeptide_auc_data)[sort_column]
    if (is.numeric(auc_column)) auc_column <- colnames(cpeptide_auc_data)[auc_column]
    
    cpeptide_auc_data$orig_row <- 1:nrow(cpeptide_auc_data)
    
    if (!is.null(sort_column)) cpeptide_auc_data <- dplyr::arrange(cpeptide_auc_data, !! rlang::sym(sort_column))

    cpeptide_auc_data$prev_visit_at_threshold <- FALSE
    
    for (id.tmp in unique(cpeptide_auc_data[[identifier_column]])) {
      rows.tmp <- which(cpeptide_auc_data[[identifier_column]]==id.tmp)
      if (length(rows.tmp) > 1) {
        for (i in 2:length(rows.tmp)) {
          cpeptide_auc_data$prev_visit_at_threshold[rows.tmp[i]] <-
            cpeptide_auc_data$prev_visit_at_threshold[rows.tmp[i-1]] |
            (!is.na(cpeptide_auc_data[rows.tmp[i-1], auc_column, drop=TRUE]) &
               (cpeptide_auc_data[rows.tmp[i-1], auc_column, drop=TRUE] <= threshold))
        }
      }
    }
    
    cpeptide_auc_data <-
      cpeptide_auc_data %>%
      dplyr::arrange(orig_row)

    cpeptide_auc_data$prev_visit_at_threshold
  }