scale_nonzero <- function(x) {
  x.nonzero <- x[x != 0]
  if(length(x.nonzero) < 2) return(x)
  x.nonzero <- x.nonzero/median(x.nonzero)
  x[x != 0] <- x.nonzero
  return(x)
}