# utilities
subset_distance <- function(dist, ind) {
  if(class(dist) != "dist")
    stop("dist must be a distance class!")
  if(!(class(ind) %in% c("character", "logical")))
    stop("ind must be a name/logical vector!")
  
  dist_matrix <- as.matrix(dist)
  if(class(ind) == "logical" & nrow(dist_matrix) != length(ind))
    stop("Dimensions of dist and ind must match!")
  if(class(ind) == "character" & !all(ind %in% rownames(dist_matrix)))
    stop("ind must be a subset of the row/col names of dist!")
  
  return(as.dist(dist_matrix[ind, ind]))
}

# plot
cowplot_title <- function(p, title, rel_heights = c(0.1, 1)) {
  title <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold")
  plot_grid(title, p, ncol =1 , rel_heights = rel_heights)
}
rotate_xaxis <- function(angle) {
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
}
no_label_xaxis <- function() {
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
}
gg_color_hue <- function(values = NULL, n = NULL) {
  if(!is.null(values)) {
    if(anyDuplicated(values)) stop("Values should be unique if provided!")
    if(!is.character(values)) stop("Values should be of character class!")
    n <- length(values)
    hues <- seq(15, 375, length = n + 1)
    colors <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colors) <- values
    return(colors)
  }
  hues <- seq(15, 375, length = n + 1)
  colors <- hcl(h = hues, l = 65, c = 100)[1:n]
}
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

fill_na <- function(variable) {
  if(class(variable) != "character") stop("variable must be character!")
  ifelse(is.na(variable), "NA", variable)
}