
COLNAMES <-
  c("sample_id", "user_id", "date", "time", "conc_nguL", "A260", "A280",
    "A260_280", "A260_230", "constant", "cursor_pos", "cursor_abs", "raw_340",
    "measurement_type", "serial_num", "config", 220:350)

SAMPLETYPES <- c("50" = "dsDNA", "40" = "RNA", "33" = "ssDNA")

COLSPEC <-
  c(rep("c", 4), rep("d", 9), rep("c", 3), rep("d", 131)) %>%
  paste0(collapse = "")


#' Plot spectrum of Nanodrop data
#' @param ndv_file A path for the ndv file
#' @export
read_ndv <- function(ndv_file) {
  sample_id <- constant <- sample_name <- user_id <-
    measurement_type <- sample_type <- NULL
  if(fs::path_ext(ndv_file) != "ndv")
    stop("The file extention must be `.ndv`.")

  module  <- readr::read_lines(ndv_file, n_max = 1L)
  if(module != "Module:\tNucleic Acid")
    stop("Module must be `Nucleic Acid`")

  ndv_header <-
    readr::read_tsv(ndv_file, n_max = 4, col_names = F, col_types = "cc") %>%
    purrr::set_names(c("name", "value"))

  ndv_df <-
    readr::read_tsv(ndv_file, skip = 4, col_names = T, col_types = COLSPEC) %>%
    purrr::set_names(COLNAMES) %>%
    dplyr::mutate(sample_name = sample_id, sample_id = dplyr::row_number()) %>%
    dplyr::mutate(sample_type = SAMPLETYPES[as.character(constant)]) %>%
    dplyr::select(sample_id, sample_name, user_id:measurement_type,
                  sample_type, dplyr::everything())

  list(file = ndv_file, header = ndv_header, data = ndv_df)
}


#' Plot spectrum of Nanodrop data
#' @param ndv_li list from read_ndv()
#' @import ggplot2
#' @export
plot_ndv_spectrum <- function(ndv_li) {
  sample_id <- sample_type <- wv_len <- value <- NULL

  ndv_df <- ndv_li[["data"]]
  plot_df <-
    ndv_df %>%
    dplyr::mutate(sample_id = as.character(sample_id)) %>%
    dplyr::select(sample_id, sample_type, "220":"350") %>%
    tidyr::pivot_longer(
      cols = -c(sample_id, sample_type),
      names_to = "wv_len",
      names_ptypes = list(wv_len = integer())
    )

  gp <-
    plot_df %>%
    ggplot(aes(wv_len, value, group = sample_id)) +
    geom_line(aes(color = sample_id)) +
    labs(title = ndv_li$file, x = "Wave length (nm)", y = "10 mm Absorbance") +
    scale_x_continuous(expand = expand_scale(mult = 0),
                       breaks = seq(220, 350, 10)) +
    theme_linedraw() +
    theme(legend.position = "right")

  ndv_li$plot_df <- plot_df
  ndv_li$gp <- gp
  ndv_li
}


