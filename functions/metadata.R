#' Fill in a pre-specified filler value for a character vector with missing values
#'
#' @param variable character with missing values
#' @param fill filler value for missingness, default "NA"
#'
#' @return
#' @export
#'
#' @examples
fill_na <- function(variable, fill = "NA") {
  if(class(variable) != "character") stop("variable must be character!")
  ifelse(is.na(variable), fill, variable)
}

#' Generate summary statistic for a variable. If it's categorical, then provide n per category;
#' if it's numeric, then provide mean and standard error.
#'
#' @param x 
#' @param variable_type 
#' @param categories 
#'
#' @return
#' @export
#'
#' @examples
summarise_variable <- function(x, 
                               variable_type,
                               categories = NULL) {
  if(all(is.na(x))) return("")
  
  if(!(variable_type %in% c("character", "numeric")))
    stop("Variable type must be either character or numeric!")
  if(class(x) != variable_type)
    stop("Class of x is different from provided by variable_type!")
  
  if(variable_type == "character") {
    if(is.null(categories))
      stop("Allowed categories must be provided in variable is character!")
    if(!all(x[!is.na(x)] %in% categories))
      stop("x has values not present in categories!")
    
    ns <- sapply(categories, function(cat) sum(x %in% cat))
    ns <- ns[ns != 0]
    return(paste0(names(ns), " " , ns) %>% 
             paste0(collapse = "/"))
  }
  
  if(variable_type == "numeric") {
    mean_variable <- mean(x, na.rm = TRUE)
    sd_variable <- sd(x, na.rm = TRUE)
    return(paste0(round(mean_variable, digits = 2), 
                  " (", 
                  round(sd_variable, digits = 2),
                  ")"))
  }
}