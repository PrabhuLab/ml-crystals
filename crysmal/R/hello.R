# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

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
