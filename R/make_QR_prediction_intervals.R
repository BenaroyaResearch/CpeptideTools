# Making QR prediction intervals
#' @export
make_QR_prediction_intervals <- function(data_to_add_QR_values_to,
                                         current_response_time,
                                         baseline_date,
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
  coef_label <- paste0("baseline: ", baseline_date, " endpoint: ", current_response_time)
  
  if (!coef_label %in% names(base_QR_models_minimal)) {
    stop(paste0(
      "No model found for '", coef_label, "'.\n",
      "Available models: ", paste(names(base_QR_models_minimal), collapse = "; ")
    ))
  }
  
  model_info <- base_QR_models_minimal[[coef_label]]
  coefs      <- model_info$coef
  sigma_hat  <- model_info$sigma_hat
  
  # ── Compute fitted placebo values and SE ──────────────────────────────────────
  newdata <- data_to_add_QR_values_to %>%
    dplyr::transmute(
      log_mean_AUC_baseline = log_baseline_mean_AUC_cpep,
      Age_At_Screening      = Age_At_Screening
    )
  
  X           <- model.matrix(~ log_mean_AUC_baseline + Age_At_Screening, data = newdata)
  fitted_vals <- as.numeric(X %*% coefs)
  
  QR_dat <- data_to_add_QR_values_to %>%
    dplyr::mutate(
      lm_placebo_estimates = fitted_vals,
      lm_estimated_te      = log_response_mean_AUC_cpep - lm_placebo_estimates,
      lm_estimated_te_se   = sigma_hat
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
    if ("Placebo" %in% groups_to_plot) {
      n_tot  <- sum(study_dat$Active_versus_Placebo == "Placebo", na.rm = TRUE)
      n_resp <- sum(study_dat$Active_versus_Placebo == "Placebo" & study_dat$Responder, na.rm = TRUE)
      pct    <- if (n_tot > 0) round(100 * n_resp / n_tot, 1) else NA
      lines  <- c(lines, paste0("Placebos called as responders: ", n_resp, " / ", n_tot, " (", pct, "%)"))
    }
    if ("Active" %in% groups_to_plot) {
      n_tot  <- sum(study_dat$Active_versus_Placebo == "Treated", na.rm = TRUE)
      n_resp <- sum(study_dat$Active_versus_Placebo == "Treated" & study_dat$Responder, na.rm = TRUE)
      pct    <- if (n_tot > 0) round(100 * n_resp / n_tot, 1) else NA
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