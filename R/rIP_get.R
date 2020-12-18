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


## ## GROUP ID
## get
#' @export
groupId = function(x) {
  UseMethod("groupId", x)
}

#' @export
groupId.DyadExperiment = function(x){
  `Group Ids` = sapply(x, attr, "groupId")
  table(`Group Ids`)
}

#' @export
groupId.DyadSession = function(x){
  attr(x,"groupId")
}

#' @export
name <- function(x) {
  attr(x,"name")
}

## set
#' @export
`groupId<-` = function(x, value){
  UseMethod("groupId<-",x)
}
#' @export
`groupId<-.DyadExperiment` <- function(x, value){
  namesx = character(length(x))
  for(i in 1:length(x)){
    groupId(x[[i]]) = value
    namesx[[i]] = name(x[[i]])
  }
  names(x) = namesx
  x
}
#' @export
`groupId<-.DyadSession` <- function(x, value){
  name = paste(value,lead0(sessionId(x)),dyadId(x),sep="_")
  attr(x, "groupId") <- value
  attr(x, "name")    <- name
  x
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
