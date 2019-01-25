#' Estimate rate of C-peptide change over time, by subject (with options!)
#'
#' This function uses linear models to estimate the rate of C-peptide AUC change over time. 
#' There are three options for how to fit models, set using the \code{model_type} argument:
#' "independent" fits a separate, independent intercept and slope for each value of
#' \code{identifier_column}; "random_effect" fits an individual-level random effect for both the
#' intercept and the slpope; "grouped_random_effect" fits an individual-level random effect for
#' the intercept, a group-level fixed effect for each unique value of \code{group_column}, and
#' an individual-level random effect around the group means.
#' @param cpeptide_auc_data data frame containing the C-peptide data. Should contain a unique subject identifer, a numeric column for the timing of visits, and C-peptide AUC values. Any non-NA values included in this data frame will be used, so filtering should be done before passing the data to this function.
#' @param model_type character, the type(s) of model to fit. Options are "independent", "random_effect", and "grouped_random_effect". "independent" uses a simple fixed-effects model with slopes and intercepts for each subject. "random_effect" uses a mixed-effects model with subject-level random intercepts and slopes. "grouped_random_effect" uses a mixed-effects model like "random_effect", but with group-level fixed effects for the slopes.
#' @param identifier_column character or numeric, the column containing the subject identifiers. Defaults to "subject".
#' @param time_column character or numeric, the column with the time variable. Should be numeric and based on a baseline visit, as intercept values only make sense if they are based on a common scale. Default is "cpeptide_study_day".
#' @param auc_column character or numeric, the column containing the C-peptide AUC values. Data should be transformed to whatever form the model is to be fit to (typically log-transformed). Defaults to "auc".
#' @param group_column character or numeric, the column containing the subject grouping. Used only if \code{model_type} is set to "grouped_random_effect"
#' @export
#' @return a list containing, for each value of \code{model_type}, an element with the model fit(s) (as $model) and and a data frame with the extracted slopes and intercepts for each subject (as $rates). Slopes are given in the units of \code{auc_column} / unit of \code{time_column}.
#' @usage \code{
#' calc_Cpeptide_rates(
#'   cpeptide_auc_data, model_type=c("independent", "random_effect"),
#'   identifier_column="subject", time_column="cpeptide_study_day",
#'   auc_column="auc",
#'   group_column
#'   )}
calc_Cpeptide_rates <-
  function(cpeptide_auc_data, model_type=c("independent", "random_effect"),
           identifier_column="subject", time_column="cpeptide_study_day",
           auc_column="auc",
           group_column) {
    
    # ensure character column names (necessary for dplyr functions)
    if (is.numeric(identifier_column)) identifier_column <- colnames(cpeptide_auc_data)[identifier_column]
    if (is.numeric(time_column)) time_column <- colnames(cpeptide_auc_data)[time_column]
    if (is.numeric(auc_column)) auc_column <- colnames(cpeptide_auc_data)[auc_column]
    if (is.numeric(group_column)) group_column <- colnames(cpeptide_auc_data)[group_column]
    
    # scale time variable for better model fitting
    time_scaled <-
      scale(cpeptide_auc_data[[time_column]], center=FALSE,
            scale=sd(cpeptide_auc_data[[time_column]], na.rm=TRUE))
    cpeptide_auc_data[[time_column]] <-
      as.vector(time_scaled)
    time.sd <- attr(time_scaled, "scaled:scale")

    
    model_type <- 
      match.arg(
        model_type, c("independent", "random_effect", "grouped_random_effect"),
        several.ok = TRUE)
    
    cpeptide_models_rates <- list()
    
    # need to figure out how to structure the output
    if ("independent" %in% model_type) {
      cpeptide_models_rates[["independent"]] <- list()
      cpeptide_models_rates[["independent"]][["model"]] <-
        lm(
          formula(
            paste(paste0(auc_column, " ~ 0"),
                  identifier_column,
                  paste0(time_column, ":", identifier_column),
                  sep=" + ")),
          data=cpeptide_auc_data)
      cpeptide_models_rates[["independent"]][["rates"]] <-
        data.frame(
          unique(cpeptide_auc_data[[identifier_column]]),
          slope=as.numeric(NA),
          intercept=as.numeric(NA),
          stringsAsFactors = FALSE)
      colnames(cpeptide_models_rates[["independent"]][["rates"]])[1] <-
        identifier_column
      for (i in 1:nrow(cpeptide_models_rates[["independent"]][["rates"]])) {
        if (
          sum(!is.na(coef(cpeptide_models_rates[["independent"]][["model"]])[
            which(str_detect(names(coef(cpeptide_models_rates[["independent"]][["model"]])),
                             cpeptide_models_rates[["independent"]][["rates"]][i, identifier_column]))])) == 2) {
          cpeptide_models_rates[["independent"]][["rates"]][i, "slope"] <-
            coef(cpeptide_models_rates[["independent"]][["model"]])[
              which(str_detect(names(coef(cpeptide_models_rates[["independent"]][["model"]])),
                               cpeptide_models_rates[["independent"]][["rates"]][i, identifier_column]) &
                      str_detect(names(coef(cpeptide_models_rates[["independent"]][["model"]])),
                                 time_column))] / time.sd
          cpeptide_models_rates[["independent"]][["rates"]][i, "intercept"] <-
            coef(cpeptide_models_rates[["independent"]][["model"]])[
              which(str_detect(names(coef(cpeptide_models_rates[["independent"]][["model"]])),
                               cpeptide_models_rates[["independent"]][["rates"]][i, identifier_column]) &
                      !str_detect(names(coef(cpeptide_models_rates[["independent"]][["model"]])),
                                  time_column))]
        }
      }
    }
      
      if ("random_effect" %in% model_type) {
        cpeptide_models_rates[["random_effect"]] <- list()
        cpeptide_models_rates[["random_effect"]][["model"]] <-
          lme4::lmer(
            formula=
              formula(
                paste(paste0(auc_column, " ~ 1"),
                      paste0("(1|", identifier_column, ")"),
                      time_column,
                      paste0("(", time_column, "-1|", identifier_column, ")"),
                      sep=" + ")),
            data=cpeptide_auc_data)
        
        cpeptide_models_rates[["random_effect"]][["rates"]] <-
          data.frame(
            unique(cpeptide_auc_data[[identifier_column]]),
            slope=as.numeric(NA),
            intercept=as.numeric(NA),
            stringsAsFactors = FALSE)
        colnames(cpeptide_models_rates[["random_effect"]][["rates"]])[1] <-
          identifier_column
        
        for (i in 1:nrow(cpeptide_models_rates[["random_effect"]][["rates"]])) {
          if (
            sum(!is.na(coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]])[
              match(
                cpeptide_models_rates[["random_effect"]][["rates"]][i, identifier_column],
                rownames(coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]])),]) == 2) {
            
            # fix from here down
            cpeptide_models_rates[["random_effect"]][["rates"]][i, "slope"] <-
              coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]][
                match(cpeptide_models_rates[["random_effect"]][["rates"]][i, identifier_column],
                      rownames(coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]])),
                time_column] / time.sd
            cpeptide_models_rates[["random_effect"]][["rates"]][i, "intercept"] <-
              coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]][
                match(cpeptide_models_rates[["random_effect"]][["rates"]][i, identifier_column],
                      rownames(coef(cpeptide_models_rates[["random_effect"]][["model"]])[[identifier_column]])),
                "(Intercept)"]
          }
        }
      }
    
    if ("grouped_random_effect" %in% model_type) {
      cpeptide_models_rates[["grouped_random_effect"]] <- list()
      cpeptide_models_rates[["grouped_random_effect"]][["model"]] <-
        lme4::lmer(
          formula=
            formula(
              paste(paste0(auc_column, " ~ 1"),
                    paste0("(1|", identifier_column, ")"),
                    time_column,
                    paste0(time_column, ":", group_column),
                    paste0("(", time_column, "-1|", identifier_column, ")"),
                    sep=" + ")),
          data=cpeptide_auc_data)
      cpeptide_models_rates[["grouped_random_effect"]][["rates"]] <-
        data.frame(
          unique(cpeptide_auc_data[[identifier_column]]),
          stringsAsFactors = FALSE)
      colnames(cpeptide_models_rates[["grouped_random_effect"]][["rates"]])[1] <-
        identifier_column
      cpeptide_models_rates[["grouped_random_effect"]][["rates"]][[group_column]] <-
        cpeptide_auc_data[[group_column]][
          match(cpeptide_models_rates[["grouped_random_effect"]][["rates"]][[identifier_column]],
                cpeptide_auc_data[[identifier_column]])]
      cpeptide_models_rates[["grouped_random_effect"]][["rates"]][["slope"]] <- as.numeric(NA)
      cpeptide_models_rates[["grouped_random_effect"]][["rates"]][["intercept"]] <- as.numeric(NA)
      
      # continue fixing from here 
      for (i in 1:nrow(cpeptide_models_rates[["grouped_random_effect"]][["rates"]])) {
        if (
          sum(!is.na(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])[
            match(
              cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, identifier_column],
              rownames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])),]) == 
          ncol(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])) {
          
          # fix from here down
          cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, "slope"] <-
            ((coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]][
              match(cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, identifier_column],
                    rownames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])),
              time_column]) +
               (sum(
                 coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]][
                   match(cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, identifier_column],
                         rownames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])),
                   str_detect(colnames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]]),
                              paste0(time_column, ":", group_column))] *
                   (str_replace_all(
                     grep(paste0(time_column, ":", group_column),
                          colnames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]]),
                          value=TRUE),
                     paste0(time_column, ":", group_column), "") %in%
                      cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, group_column])))) / time.sd
                
          cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, "intercept"] <-
            coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]][
              match(cpeptide_models_rates[["grouped_random_effect"]][["rates"]][i, identifier_column],
                    rownames(coef(cpeptide_models_rates[["grouped_random_effect"]][["model"]])[[identifier_column]])),
              "(Intercept)"]
        }
      }
    }
    cpeptide_models_rates
  }