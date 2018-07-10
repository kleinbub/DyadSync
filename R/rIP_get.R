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
session = sessionId = function(dyadSession){
  attr(dyadSession,"sessionId")
}
#' @export
dyad = dyadId = function(dyadSession){
  attr(dyadSession,"dyadId")
} 
#' @export
group.DyadSession = groupId = function(dyadSession){
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

#' @export
color <-function(x){
  attr(x,"col")
}
#' @export
lty <-function(x){
  attr(x,"lty")
}
#' @export
lwd <-function(x){
  # attr(x,"lwd")
}
