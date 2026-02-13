## The Pooled Adjacent Violators Algorithm
pava <- function(x, wt = rep(1, length(x))) {
  n = length(x)
  if (n <= 1) {
    return(x)
  }
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets = 1:n
  repeat {
    viol = (as.vector(diff(x)) < 0)
    if (!(any(viol))) {
      break
    }
    i = min((1:(n - 1))[viol])
    lvl1 = lvlsets[i]
    lvl2 = lvlsets[i + 1]
    ilvl = (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] = sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] = lvl1
  }
  return(x)
}