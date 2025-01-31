#' Get all functions that are used in a formula  `expr`.
#'
#' @param expr a formula object
#'
#' @return A character vector with the extracted names.
#' @export
all_funs = function(expr){
  all.names(expr, unique=TRUE)[!(all.names(expr, unique=TRUE) %in% all.vars(expr))]
}
