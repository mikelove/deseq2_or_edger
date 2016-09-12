library("GenomicRanges")
library("Biobase")
#library("BiocParallel")
library("SummarizedExperiment")
load("bottomly_sumexp.RData")
bottomly <- updateObject(bottomly)
randomSubsets <- read.table("random_subsets.txt",strings=FALSE)
se <- bottomly[,match(randomSubsets[1,],colnames(bottomly))]
colData(se)$run <- colnames(se)
eset <- ExpressionSet(assay(se),
                      AnnotatedDataFrame(as.data.frame(colData(se))))
pData(eset)$condition <- pData(eset)$strain
levels(pData(eset)$condition) <- c("A","B")

library("DESeq")
library("DESeq2")
library("edgeR")
library("limma")
source("runScripts.R")

algos <- list("DESeq2"=runDESeq2,"edgeR"=runEdgeR,"edgeRQL"=runEdgeRQL,"limma-voom"=runVoom)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

resTest <- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()

nreps <- 30

#register(MulticoreParam(workers=8,verbose=TRUE))

nreps <- 3

res <- lapply(1:nreps, function(i) {   
  print(i)
  testSet <- as.character(randomSubsets[i,1:6])
  heldOutSet <- as.character(randomSubsets[i,-(1:6)])
  eTest <- eset[,testSet]
  eHeldout <- eset[,heldOutSet]
  st <- system.time({
    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](eTest))
    resHeldout0 <- lapply(namesAlgos, function(n) algos[[n]](eHeldout))
  })
  print(paste(round(unname(st[3])),"seconds"))
  resTest <- as.data.frame(lapply(resTest0, function(z) z$padj))
  resHeldout <- as.data.frame(lapply(resHeldout0, function(z) z$padj))
  lfcTest <- as.data.frame(lapply(resTest0, function(z) z$beta))
  lfcHeldout <- as.data.frame(lapply(resHeldout0, function(z) z$beta))
  list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout)
})

resTest <- lapply(res, "[[", "resTest")
resHeldout <- lapply(res, "[[", "resHeldout")
lfcTest <- lapply(res, "[[", "lfcTest")
lfcHeldout <- lapply(res, "[[", "lfcHeldout")

save(resTest,resHeldout,lfcTest,lfcHeldout,namesAlgos,file="sensitivityPrecision.RData")

