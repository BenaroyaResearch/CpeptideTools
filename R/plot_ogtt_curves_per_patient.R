#' Plot OGTT curves per patient over time
#'
#' This function outputs OGTT curves over time as a multi-page PDF, where each subplot shows
#' the set of curves for a given patient across different draw dates. Each curve represents
#' one draw date, and is colored by the number of days elapsed since that patient's earliest
#' recorded draw date. Patients are randomly shuffled and divided into pages (blocks) for
#' legibility.
#'
#' @param data data frame in long format, where each row represents a specific patient at a
#'   given draw date, OGTT type, and timepoint. Timepoints correspond to the collection times
#'   within a 2-hour MMTT course (e.g. \code{"PEP30"} for C-peptide collected at 30 minutes).
#' @param ogtt_subset character vector, the OGTT data types to include in the plot (e.g.
#'   C-peptide, glucose). Filters rows where \code{ogtt_type_var} matches one of these values.
#'   Also used as the y-axis label.
#' @param ogtt_type_var character, the name of the column identifying the OGTT data type for
#'   each row. Defaults to \code{"ogtt_type"}.
#' @param id_var character, the name of the column containing the unique patient identifier.
#' @param draw_dt_var character, the name of the column containing the draw date. Defaults to
#'   \code{"DRAW_DT"}.
#' @param time_var character, the name of the column containing the numeric OGTT timepoint
#'   within the 2-hour collection period. Defaults to \code{"time_numeric"}.
#' @param value_var character, the name of the column containing the OGTT measurement values.
#'   Defaults to \code{"Ogtt_data"}.
#' @param file_name character, the file path where the output PDF will be saved.
#' @param blocks numeric, the number of pages (chunks) to divide patients across in the output
#'   PDF. Patients are randomly shuffled before chunking. Defaults to 16.
#' @param seed numeric, random seed for reproducibility of the patient shuffling. Defaults to 1.
#' @param pdf_width numeric, width of the output PDF in inches. Defaults to 12.
#' @param pdf_height numeric, height of the output PDF in inches. Defaults to 8.
#' @param hline_value numeric, y-intercept of a reference horizontal line drawn in red on each
#'   subplot. Defaults to 0.15.
#' @import dplyr
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @importFrom stats as.formula
#' @export
#' @return invisibly, the filtered and processed long-format data frame used for plotting,
#'   with one row per patient per draw date per timepoint. As a side effect, writes a
#'   multi-page PDF to \code{file_name}.
#' @usage
#' plot_ogtt_curves_per_patient(
#'   data,
#'   ogtt_subset,
#'   ogtt_type_var = "ogtt_type",
#'   id_var,
#'   draw_dt_var   = "DRAW_DT",
#'   time_var      = "time_numeric",
#'   value_var     = "Ogtt_data",
#'   file_name,
#'   blocks        = 16,
#'   seed          = 1,
#'   pdf_width     = 12,
#'   pdf_height    = 8,
#'   hline_value   = 0.15)
plot_ogtt_curves_per_patient <- function(data,
                                         ogtt_subset,
                                         ogtt_type_var = "ogtt_type",
                                         id_var,
                                         draw_dt_var = "DRAW_DT",
                                         time_var = "time_numeric",
                                         value_var = "Ogtt_data",
                                         file_name,
                                         blocks = 16,
                                         seed = 1,
                                         pdf_width = 12,
                                         pdf_height = 8,
                                         hline_value = 0.15) {
  
  # Check that requested columns exist
  required_vars <- c(id_var, draw_dt_var, time_var, value_var, ogtt_type_var)
  
  missing_vars <- setdiff(required_vars, names(data))
  
  if (length(missing_vars) > 0) {
    stop(
      paste0(
        "ERROR! The following required columns are missing from 'data': ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }
  
  # Keep only the requested OGTT type and create plotting variables
  ogtt_subset_long <- data %>%
    dplyr::filter(.data[[ogtt_type_var]] %in% ogtt_subset) %>%
    dplyr::mutate(
      DRAW_DT_date = as.Date(.data[[draw_dt_var]]),
      curve_id = paste(.data[[id_var]], DRAW_DT_date, sep = "_")
    ) %>%
    dplyr::group_by(.data[[id_var]]) %>%
    dplyr::mutate(
      dates_numeric = as.integer(DRAW_DT_date - min(DRAW_DT_date, na.rm = TRUE))
    ) %>%
    dplyr::ungroup()
  
  # Check that data remain after filtering
  if (nrow(ogtt_subset_long) == 0) {
    stop("ERROR! No rows remain after filtering to the requested ogtt_subset.")
  }
  
  # Remove rows with missing values needed for plotting
  ogtt_subset_long <- ogtt_subset_long %>%
    dplyr::filter(
      !is.na(.data[[id_var]]),
      !is.na(DRAW_DT_date),
      !is.na(.data[[time_var]]),
      !is.na(.data[[value_var]])
    )
  
  if (nrow(ogtt_subset_long) == 0) {
    stop("ERROR! No complete rows remain after removing missing plotting values.")
  }
  
  set.seed(seed)
  
  grDevices::pdf(file_name, width = pdf_width, height = pdf_height)
  
  # Close the PDF device when the function exits (AK: neat!).
  on.exit(grDevices::dev.off(), add = TRUE)
  
  # Shuffle IDs once
  patient_ids <- unique(ogtt_subset_long[[id_var]])
  patient_ids <- sample(patient_ids, length(patient_ids))
  
  # Split into roughly equal chunks
  id_chunks <- split(
    patient_ids,
    cut(seq_along(patient_ids), breaks = blocks, labels = FALSE)
  )
  
  for (i in seq_len(length(id_chunks))) {
    chunk_ids <- id_chunks[[i]]
    
    if (length(chunk_ids) == 0) {
      next
    }
    
    random_set <- ogtt_subset_long %>%
      dplyr::filter(.data[[id_var]] %in% chunk_ids)
    
    current_breaks <- sort(unique(random_set[[time_var]]))
    
    p <- random_set %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data[[time_var]],
          y = .data[[value_var]],
          group = curve_id,
          color = dates_numeric
        )
      ) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(stats::as.formula(paste("~", id_var))) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(breaks = current_breaks) +
      ggplot2::geom_hline(yintercept = hline_value, col = "red") +
      ggplot2::scale_color_viridis_c(name = "Days since first draw") +
      ggplot2::labs(
        x = "OGTT time (min)",
        y = ogtt_subset
      )
    
    print(p)
  }
  
  invisible(ogtt_subset_long)
}