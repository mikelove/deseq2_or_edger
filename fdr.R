library("ggplot2")
library("reshape")
getSensitivity <- function(alpha, alphaOut) {
  t(sapply(1:nreps, function(i) sapply(namesAlgos, function(algo) {
    sigHeldout <- resHeldout[[i]][[algo]] < alphaOut
    mean((resTest[[i]][[algo]] < alpha)[sigHeldout])
  })))
}
getFDR <- function(alpha, alphaOut) {
  t(sapply(1:nreps, function(i) sapply(namesAlgos, function(algo) {
    sigTest <- resTest[[i]][[algo]] < alpha
    if (sum(sigTest) == 0) return(0)
    mean((resHeldout[[i]][[algo]] > alphaOut)[sigTest])
  })))
}
nreps <- length(resTest)
sens <- getSensitivity(.1,.1)
fdr <- getFDR(.1, .1)
boxplot(fdr)
