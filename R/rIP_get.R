#' Title
#'
#' @param x 


#' @export
start.DyadSignal <- function (x) {attr(x,"start")}
#' @export
end.DyadSignal <- function (x) {attr(x,"end")}
#' @export

#' @export
s1Name <- function(x) {attr(x,"s1Name")}
#' @export
s2Name <- function(x) {attr(x,"s2Name")}


## ## SESSION ID
## get
#' @export
sessionId = function(x) {
  UseMethod("sessionId", x)
}

#' @export
sessionId.DyadExperiment = function(x){
  `Session Ids` = sapply(x, attr, "sessionId")
  # table(`Session Ids`)
}

#' @export
sessionId.DyadSession = function(x){
  attr(x,"sessionId")
}

#' @export
sessionId.DyadSignal = function(x){
  attr(x,"sessionId")
}
## ## DYAD ID
## get

#' @export
dyadId = function(x) {
  UseMethod("dyadId", x)
}

#' @export
dyadId.DyadExperiment = function(x){
  `Dyad Ids` = sapply(x, attr, "dyadId")
  # table(`Dyad Ids`)
}

#' @export
dyadId.DyadSession = function(x){
  attr(x,"dyadId")
}

#' @export
dyadId.DyadSignal = function(x){
  attr(x,"dyadId")
}




## ## GROUP ID
## get
#' @export
groupId = function(x) {
  UseMethod("groupId", x)
}

#' @export
groupId.DyadExperiment = function(x){
  `Group Ids` = sapply(x, attr, "groupId")
  # table(`Group Ids`)
}

#' @export
groupId.DyadSession = function(x){
  attr(x,"groupId")
}

#' @export
groupId.DyadSignal = function(x){
  attr(x,"groupId")
}

#' @export
name <- function(x) {
  attr(x,"name")
}

#' @export
UID = function(x){
  if(!all(c("dyadId","sessionId","groupId") %in% names(attributes(x)))){
    stop("Only objects of class DyadSession, or DyadSignal 
         have DyadSync::UID methods")
  }
  paste(attr(x,"groupId"), attr(x,"dyadId"), attr(x,"sessionId"), sep="_" )
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
