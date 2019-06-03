#' Extract PCoA coordinates, as well as percentage variability explained from an 
#' MDS fit
#'
#' @param fit_pcoa an MDS fit from running ape::pcoa
#' @param n_axis number of top axes to extract
#'
#' @return
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
extract_pcoa <- function(fit_pcoa, n_axis) {
  varExplained <- fit_pcoa$values$Relative_eig[1:n_axis]
  tb_axes <- gather_axes(fit_pcoa$vectors, n_axis)
  tb_varExplained <- (1:(n_axis/2)) %>% 
    purrr::map_dfr(function(i.axis) {
      tibble::tibble(axes = paste0(i.axis*2 - 1, "vs.", i.axis*2),
                     x.perc.var = varExplained[i.axis*2 - 1],
                     y.perc.var = varExplained[i.axis*2])
    })
  tb_axes %>% dplyr::left_join(tb_varExplained, by = "axes") %>% return()
}

#' Gather axes of a data frame, for plotting multiple top axes of a PCoA analysis
#' 
#' @param data data frame of PCoA coordinates, each column corresponding to a top PCoA
#' axis
#' @param n_axis number of top axis to gather
#'
#' @return
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
gather_axes <- function(data, n_axis) {
  if(n_axis %% 2 != 0)
    stop("n_axis must be even number!")
  tb_axes <- (1:(n_axis/2)) %>% 
    purrr::map_dfr(function(i.axis) {
      tibble::tibble(x.axis = data[, i.axis*2 - 1, drop = TRUE],
                     y.axis = data[, i.axis*2, drop = TRUE],
                     axes = paste0(i.axis*2 - 1, "vs.", i.axis*2))
    })
  if(!is.null(row.names(data)))
    tb_axes$rowname <- rep(row.names(data), n_axis/2)
  return(tb_axes)
}