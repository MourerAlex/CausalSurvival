#' Null-Coalescing Operator
#'
#' Returns `x` if not NULL, otherwise `y`.
#'
#' @param x,y Values to coalesce.
#' @return `x` if not NULL, else `y`.
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
