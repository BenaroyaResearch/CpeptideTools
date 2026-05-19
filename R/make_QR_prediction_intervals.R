#' Generate quantitative response (QR) prediction interval plots for treated and placebo subjects
#'
#' This function computes prediction intervals around each subject's quantitative response (QR)
#' estimate and generates waterfall plots showing individual-level QR values with uncertainty
#' bounds. QR measures the difference between a subject's observed log-transformed AUC response
#' and the trajectory expected under placebo; prediction intervals around this estimate are used
#' to classify subjects as responders to therapy. Subjects whose entire prediction interval lies
#' above zero are classified as responders, indicating their response exceeds placebo expectation.
#' Plots can be stratified by study and prediction interval level, and optionally saved to PDF
#' or PNG files.
#'
#' @param data_to_add_QR_values_to data frame containing subject-level data, including AUC
#'   values at baseline and the response timepoint, age at baseline, and columns named
#'   \code{ID}, \code{Study_label}, and \code{Active_versus_Placebo}.
#' @param current_response_time character, the label of the response/endpoint timepoint used
#'   to select the correct pre-fitted model (e.g. \code{"12 months"}).
#' @param baseline_date character, the label of the baseline visit used to select the correct
#'   pre-fitted model (e.g. \code{"randomization"}).
#' @param mean_auc_baseline_col character, the name of the column containing mean AUC at baseline.
#' @param mean_auc_response_col character, the name of the column containing mean AUC at the
#'   response timepoint.
#' @param age_at_baseline_col character, the name of the column containing age at baseline.
#' @param QR_prediction_levels character vector, the prediction interval levels to compute and
#'   plot. Must be a subset of \code{c("55", "60", "65", "70", "80", "90", "95")}. A subject
#'   is classified as a responder at a given level if the lower bound of their QR prediction
#'   interval at that level exceeds zero. Defaults to \code{c("65", "80", "95")}.
#' @param groups_to_plot character vector, which treatment groups to include in the responder
#'   count annotation on each plot. Annotation for a group is suppressed if no subjects from
#'   that group are present in the data. Defaults to \code{c("Active", "Placebo")}.
#' @param study_label character or NULL, an optional study label to restrict plotting to a
#'   single study. If NULL (default), plots are generated for all studies in the data.
#' @param ylim numeric vector of length 2 or NULL, y-axis limits for the plots (e.g.
#'   \code{c(-0.5, 1)}). If NULL (default), limits are determined automatically.
#' @param width numeric, width of the output plot in inches. Defaults to 11.
#' @param height numeric, height of the output plot in inches. Defaults to 8.5.
#' @param base_size numeric, base font size for the plot theme. Defaults to 14.
#' @param positive_color character, color used for subjects classified as responders (i.e.
#'   whose entire QR prediction interval lies above zero). Defaults to \code{"darkred"}.
#' @param nonpositive_color character, color used for subjects not classified as responders.
#'   Defaults to \code{"grey70"}.
#' @param pdf_file character or NULL, file path for saving all plots to a single multi-page
#'   PDF. If NULL (default), no PDF is saved.
#' @param png_prefix character or NULL, file path prefix for saving individual plots as PNG
#'   files. If NULL (default), no PNG files are saved. Files are named using the prefix,
#'   study label, and prediction interval level.
#' @param png_units character, units for PNG dimensions, passed to \code{ggplot2::ggsave}.
#'   Defaults to \code{"in"}.
#' @param png_dpi numeric, resolution for PNG output in dots per inch. Defaults to 300.
#' @param show_sd logical, whether to include the placebo model's residual SD in the plot
#'   annotation. The residual SD is used to construct the prediction intervals. Defaults to TRUE.
#' @import dplyr
#' @import ggplot2
#' @import purrr
#' @import stringr
#' @import magrittr
#' @export
#' @return invisibly, a named list with three elements: \code{plots} (a named list of
#'   \code{ggplot} objects, one per study/prediction-level combination), \code{png_files} (a
#'   character vector of any PNG file paths written), and \code{plot_data} (the long-format
#'   data frame used for plotting, with one row per subject per prediction interval level,
#'   including the QR estimate, lower and upper interval bounds, and responder classification).
#' @usage
#' make_QR_prediction_intervals(
#'   data_to_add_QR_values_to,
#'   current_response_time,
#'   baseline_date,
#'   mean_auc_baseline_col,
#'   mean_auc_response_col,
#'   age_at_baseline_col,
#'   QR_prediction_levels = c("65", "80", "95"),
#'   groups_to_plot       = c("Active", "Placebo"),
#'   study_label          = NULL,
#'   ylim                 = NULL,
#'   width                = 11,
#'   height               = 8.5,
#'   base_size            = 14,
#'   positive_color       = "darkred",
#'   nonpositive_color    = "grey70",
#'   pdf_file             = NULL,
#'   png_prefix           = NULL,
#'   png_units            = "in",
#'   png_dpi              = 300,
#'   show_sd              = TRUE)
make_QR_prediction_intervals <- function(data_to_add_QR_values_to,
                                         current_response_time,
                                         baseline_date,
                                         mean_auc_baseline_col,
                                         mean_auc_response_col,
                                         age_at_baseline_col,
                                         QR_prediction_levels = c("65", "80", "95"),
                                         groups_to_plot       = c("Active", "Placebo"),
                                         study_label          = NULL,
                                         ylim                 = NULL,
                                         width                = 11,
                                         height               = 8.5,
                                         base_size            = 14,
                                         positive_color       = "darkred",
                                         nonpositive_color    = "grey70",
                                         pdf_file             = NULL,
                                         png_prefix           = NULL,
                                         png_units            = "in",
                                         png_dpi              = 300,
                                         show_sd              = TRUE) {
  
  # ── Input validation ─────────────────────────────────────────────────────────
  all_qr_levels <- c("95", "90", "80", "70", "65", "60", "55")
  if (!all(QR_prediction_levels %in% all_qr_levels)) {
    stop("QR_prediction_levels must contain only: '95', '90', '80', '70', '65', '60', '55'.")
  }
  if (!is.null(ylim)) {
    if (!is.numeric(ylim) || length(ylim) != 2 || any(!is.finite(ylim))) {
      stop("ylim must be NULL or a numeric vector of length 2, e.g. c(-0.5, 1).")
    }
  }
  
  z_values <- c(
    "95" = 1.960, "90" = 1.645, "80" = 1.281, "70" = 1.036,
    "65" = 0.935, "60" = 0.841, "55" = 0.755
  )
  z_values <- z_values[QR_prediction_levels]
  
  # ── Pull model info from sysdata ─────────────────────────────────────────────
  model_lookup_key <- paste0("baseline: ", baseline_date, " endpoint: ", current_response_time)
  
  if (!model_lookup_key %in% names(base_QR_models_minimal)) {
    stop(paste0(
      "No model found for '", model_lookup_key, "'.\n",
      "Available models: ", paste(names(base_QR_models_minimal), collapse = "; ")
    ))
  }
  
  selected_model <- base_QR_models_minimal[[model_lookup_key]]
  coefs          <- selected_model$coef
  sigma_hat      <- selected_model$sigma_hat
  
  # ── Compute fitted placebo values and SE ─────────────────────────────────────
  design_matrix <- model.matrix(
    ~ log_mean_AUC_baseline + Age_At_Screening,
    data = data_to_add_QR_values_to %>%
      dplyr::transmute(
        log_mean_AUC_baseline = log(.data[[mean_auc_baseline_col]] + 1),
        log_mean_AUC_response = log(.data[[mean_auc_response_col]] + 1),
        Age_At_Screening      = .data[[age_at_baseline_col]]
      )
  )
  fitted_placebo_vals <- as.numeric(design_matrix %*% coefs)
  
  QR_dat <- data_to_add_QR_values_to %>%
    dplyr::mutate(
      log_mean_AUC_baseline = log(.data[[mean_auc_baseline_col]] + 1),
      log_mean_AUC_response = log(.data[[mean_auc_response_col]] + 1),
      lm_placebo_estimates  = fitted_placebo_vals,
      lm_estimated_te       = log_mean_AUC_response - lm_placebo_estimates,
      lm_estimated_te_se    = sigma_hat
    )
  
  # ── Optional study filter ─────────────────────────────────────────────────────
  if (!is.null(study_label)) {
    QR_dat <- QR_dat %>% dplyr::filter(Study_label == study_label)
  }
  if (nrow(QR_dat) == 0) {
    stop("No data remained after filtering. Check study_label or input data.")
  }
  
  # ── Build one row per (patient × QR level) ────────────────────────────────────
  plot_dat <- purrr::map_dfr(names(z_values), function(nm) {
    z_val <- z_values[[nm]]
    QR_dat %>%
      dplyr::mutate(
        QR_level      = paste0("QR ", nm, " prediction interval"),
        QR_estimate   = lm_estimated_te,
        QR_lower      = lm_estimated_te - z_val * lm_estimated_te_se,
        QR_upper      = lm_estimated_te + z_val * lm_estimated_te_se,
        Responder     = QR_lower > 0,
        interval_sign = ifelse(QR_lower > 0, "positive", "nonpositive"),
        point_sign    = ifelse(QR_estimate > 0, "positive", "nonpositive")
      ) %>%
      dplyr::select(
        ID, Study_label, Active_versus_Placebo,
        QR_level, QR_estimate, QR_lower, QR_upper,
        Responder, interval_sign, point_sign
      )
  }) %>%
    dplyr::filter(is.finite(QR_estimate), is.finite(QR_lower), is.finite(QR_upper)) %>%
    dplyr::mutate(
      Active_versus_Placebo = factor(
        dplyr::case_when(
          Active_versus_Placebo == "Active"  ~ "Treated",
          Active_versus_Placebo == "Placebo" ~ "Placebo",
          TRUE ~ as.character(Active_versus_Placebo)
        ),
        levels = c("Placebo", "Treated")
      ),
      interval_sign = factor(interval_sign, levels = c("nonpositive", "positive")),
      point_sign    = factor(point_sign,    levels = c("nonpositive", "positive"))
    )
  
  if (nrow(plot_dat) == 0) {
    stop("No finite QR intervals could be constructed from the data.")
  }
  
  # ── Annotation helper ─────────────────────────────────────────────────────────
  build_annotation_text <- function(study_dat, sigma_hat, groups_to_plot, show_sd) {
    lines <- character(0)
    if (show_sd) {
      lines <- c(lines, paste0("Residual SD: ", round(sigma_hat, 3)))
    }
    if ("Placebo" %in% groups_to_plot && any(study_dat$Active_versus_Placebo == "Placebo", na.rm = TRUE)) {
      n_tot  <- sum(study_dat$Active_versus_Placebo == "Placebo", na.rm = TRUE)
      n_resp <- sum(study_dat$Active_versus_Placebo == "Placebo" & study_dat$Responder, na.rm = TRUE)
      pct    <- round(100 * n_resp / n_tot, 1)
      lines  <- c(lines, paste0("Placebos called as responders: ", n_resp, " / ", n_tot, " (", pct, "%)"))
    }
    if ("Active" %in% groups_to_plot && any(study_dat$Active_versus_Placebo == "Treated", na.rm = TRUE)) {
      n_tot  <- sum(study_dat$Active_versus_Placebo == "Treated", na.rm = TRUE)
      n_resp <- sum(study_dat$Active_versus_Placebo == "Treated" & study_dat$Responder, na.rm = TRUE)
      pct    <- round(100 * n_resp / n_tot, 1)
      lines  <- c(lines, paste0("Treatment group responders: ", n_resp, " / ", n_tot, " (", pct, "%)"))
    }
    paste(lines, collapse = "\n")
  }
  
  # ── Plotting ──────────────────────────────────────────────────────────────────
  qr_levels_ordered <- paste0("QR ", all_qr_levels, " prediction interval")
  qr_levels_ordered <- qr_levels_ordered[all_qr_levels %in% QR_prediction_levels]
  
  png_files <- character(0)
  plot_list <- list()
  
  make_one_plot <- function(study_dat, plot_title) {
    study_dat <- study_dat %>%
      dplyr::arrange(dplyr::desc(QR_estimate)) %>%
      dplyr::mutate(patient_order = dplyr::row_number())
    
    annotation_text <- build_annotation_text(study_dat, sigma_hat, groups_to_plot, show_sd)
    x_annot <- if (nrow(study_dat) > 0) max(study_dat$patient_order) else 1
    y_annot <- if (!is.null(ylim)) ylim[2] else max(study_dat$QR_upper, na.rm = TRUE)
    
    p <- ggplot2::ggplot(study_dat, ggplot2::aes(x = patient_order, y = QR_estimate)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_linerange(
        ggplot2::aes(ymin = QR_lower, ymax = QR_upper, color = interval_sign),
        linewidth = 0.4
      ) +
      ggplot2::geom_point(ggplot2::aes(color = point_sign), size = 1.5) +
      ggplot2::scale_color_manual(
        values = c("nonpositive" = nonpositive_color, "positive" = positive_color),
        guide  = "none"
      ) +
      ggplot2::coord_cartesian(ylim = ylim) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::labs(title = plot_title, x = "Patients", y = "QR value") +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        plot.title   = ggplot2::element_text(hjust = 0.5)
      )
    
    if (nzchar(annotation_text)) {
      p <- p + ggplot2::annotate(
        "text", x = x_annot, y = y_annot, label = annotation_text,
        hjust = 1, vjust = 1.2, size = base_size / 3.2
      )
    }
    p
  }
  
  if (!is.null(pdf_file) && nzchar(pdf_file)) {
    grDevices::pdf(pdf_file, width = width, height = height, onefile = TRUE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  all_studies <- if (!is.null(study_label)) study_label else unique(plot_dat$Study_label)
  
  for (current_study in all_studies) {
    for (current_qr in qr_levels_ordered) {
      
      QR_val    <- stringr::str_split_i(current_qr, " ", 2)
      study_dat <- plot_dat %>%
        dplyr::filter(Study_label == current_study, QR_level == current_qr)
      
      if (nrow(study_dat) == 0) next
      
      plot_title <- paste0(current_study, " PI", QR_val)
      p          <- make_one_plot(study_dat, plot_title)
      list_key   <- paste(current_study, current_qr, sep = " | ")
      plot_list[[list_key]] <- p
      
      print(p)
      
      if (!is.null(png_prefix) && nzchar(png_prefix)) {
        current_png <- paste0(
          png_prefix, "_",
          gsub("[^A-Za-z0-9]+", "_", current_study),
          "_QR_", QR_val, ".png"
        )
        ggplot2::ggsave(current_png, plot = p, width = width, height = height,
                        units = png_units, dpi = png_dpi)
        png_files <- c(png_files, current_png)
      }
    }
  }
  
  invisible(list(
    plots     = plot_list,
    png_files = png_files,
    plot_data = plot_dat
  ))
}