#' Title
#'
#' @param x 


#' @export
start.DyadSignal <- function (x) {attr(x,"start")}
#' @export
end.DyadSignal <- function (x) {attr(x,"end")}
#' @export
duration <- function(x) {attr(x,"duration")}
#' @export
s1Name <- function(x) {attr(x,"s1Name")}
#' @export
s2Name <- function(x) {attr(x,"s2Name")}

#' @export
session = function(dyadSession){
  attr(dyadSession,"sessionId")
}
#' @export
id = function(dyadSession){
  attr(dyadSession,"dyadId")
} 

#' @export
name <- function(x) {
  attr(x,"name")
}


#' @export
sampRate <-function(x) {
  attr(x,"sampRate")
}