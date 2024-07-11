#' Hello World Function
#' @param x A character string
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' hello("world")
#' \dontrun{
#' hello("Don")
#' }
hello <- function(x) {
  print(paste0("Hello, ", x, "!"))
}
