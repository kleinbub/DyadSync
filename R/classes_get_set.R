#' Safe Object Attribute Lists
#' @description a wrapper for \code{\link[base]{attributes}} which accesses and sets an object's attributes
#' with the exception of protected ones, currently: \code{c("names", "comment", "dim", "dimnames", "row.names", "tsp")}
#' @param x an object
#' @param value an appropriate named list of attributes, or NULL.
#' @return asd
#' @examples 
#' a <- list("a"=1,"b"=2,"c"=3)
#' attr(a, "foo") <- "bar"
#' #extracts as well names
#' attributes(a) 
#' #extracts only foo
#' classAttr(a) 
#' #overwrites only foo
#' classAttr(a) <- list("bar"="foo") 
#' 
#' @export
#' 
classAttr = function(x){
  attributes(x)[!names(attributes(x)) %in% LOCK_ATTR]
}
#' @rdname classAttr
#' @export
`classAttr<-` = function(x,value){
  if(!is.list(value)||is.null(value)) stop("attributes must be a list or NULL")
  if(!is.null(attr(x,"tsp")))
    LOCK_ATTR = c(LOCK_ATTR,"tsp") #don't overwrite new tsp, but allow inheriting if missing
  value <- value[!names(value)%in%LOCK_ATTR]
  attributes(x) <- c(attributes(x)[LOCK_ATTR],value)
  x
}
LOCK_ATTR = c("names", "dim", "dimnames", "row.names")

#' clone attributes
#'
#' @param from 
#' @param to 
#'
#' @return the object 'to' with the attributes of 'from'. This does not overwrites
#' protected attributes: "names", "dim", "dimnames", "row.names", "tsp"
#' @export

cloneAttr = function(from, to){
  classAttr(to) <- classAttr(from)
  to
}




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
sessionId.default = function(x){
  if("sessionId" %in% names(attributes(x)))
    attr(x,"sessionId")
  else stop("This object doesn't have a sessionId attribute.\n")
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
dyadId.default = function(x){
  if("dyadId" %in% names(attributes(x)))
    attr(x,"dyadId")
  else stop("This object doesn't have a dyadId attribute.\n")
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
groupId.default = function(x){
  if("groupId" %in% names(attributes(x)))
    attr(x,"groupId")
  else stop("This object doesn't have a groupId attribute.\n")
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

#' @export
frequency.DyadSignal = function(x) {
  attr(x,"SR")
}
