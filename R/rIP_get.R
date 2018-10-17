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
sessionId = function(dyadSession){
  attr(dyadSession,"sessionId")
}
#' @rdname sessionId
#' @export
session = sessionId

#' @aliases dyad dyadId
#' @export
dyadId = function(dyadSession){
  attr(dyadSession,"dyadId")
} 
#' @export
dyad = dyadId
#' @export
groupId = function(dyadSession){
  attr(dyadSession,"dyadId")
} 
#' @export
group.DyadSession = groupId

#' @export
name <- function(x) {
  attr(x,"name")
}

#' @export
sampRate <-function(x) {
  if(is.ts(x)) frequency(x)
  else attr(x,"sampRate")
}

#' @export
frequency.DyadSignal = function(x){attr(x,"sampRate")}

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
