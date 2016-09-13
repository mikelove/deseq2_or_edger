load("sensFDR.rda")
library("ggplot2")
library("reshape")

getSensitivity <- function(alpha) {
  t(sapply(1:nreps, function(i) sapply(namesAlgos, function(algo) {
    sigHeldout <- resHeldout[[i]][[algo]] < alpha
    if (sum(sigHeldout) == 0) return(0)
    mean((resTest[[i]][[algo]] < alpha)[sigHeldout])
  })))
}
getFDR <- function(alpha) {
  t(sapply(1:nreps, function(i) sapply(namesAlgos, function(algo) {
    sigTest <- resTest[[i]][[algo]] < alpha
    if (sum(sigTest) == 0) return(0)
    mean((resHeldout[[i]][[algo]] > alpha)[sigTest])
  })))
}

nreps <- length(resTest)
sens <- getSensitivity(.1)
fdr <- getFDR(.1)
boxplot(fdr)
boxplot(sens)

## fdr.m <- melt(as.data.frame(fdr))
## sens.m <- melt(as.data.frame(sens))
## names(fdr.m) <- c("algo","FDR")
## data <- cbind(fdr.m, sens=sens.m[,2])
## ggplot(data, aes(FDR, sens, col=algo)) + geom_point() + xlim(0,.5) + ylim(0,1)

getOverlap <- function(a, b, alpha) {
  out <- sapply(1:nreps, function(i) {
    a.padj <- resHeldout[[i]][[a]]
    b.padj <- resHeldout[[i]][[b]]
    over <- sum(a.padj < alpha & b.padj < alpha)
    c(over/sum(a.padj < alpha), over/sum(b.padj < alpha))
  })
  out <- t(out)
  colnames(out) <- c(a,b)
  out
}
       
boxplot(getOverlap("DESeq2","edgeR",.1), ylim=c(0,1))
boxplot(getOverlap("DESeq2","edgeRQL",.1), ylim=c(0,1))
boxplot(getOverlap("DESeq2","limma.voom",.1), ylim=c(0,1))
boxplot(getOverlap("edgeR","edgeRQL",.1), ylim=c(0,1))
boxplot(getOverlap("edgeR","limma.voom",.1), ylim=c(0,1))
boxplot(getOverlap("edgeRQL","limma.voom",.1), ylim=c(0,1))
