# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' The amazing hello world function
#' @name hello
#' @description
#' A simple function that prints a welcome
#' @return Prints "Hello, world!"
#' @examples
#' hello()
#' @export

hello <- function() {
  print("Hello, world!")
}
