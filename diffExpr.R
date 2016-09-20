library("Biobase")
library("SummarizedExperiment")

# http://www-huber.embl.de/DESeq2paper/data/bottomly_sumexp.RData
load("bottomly_sumexp.RData")
bottomly <- updateObject(bottomly)

# http://www-huber.embl.de/DESeq2paper/script//bottomly/random_subsets.txt
randomSubsets <- read.table("random_subsets.txt",strings=FALSE)

se <- bottomly[,match(randomSubsets[1,],colnames(bottomly))]
colData(se)$run <- colnames(se)
eset <- ExpressionSet(assay(se),
                      AnnotatedDataFrame(as.data.frame(colData(se))))

###################################################
# load Harold's count table
counts <- readRDS("counts_table_bottomly.rds")
# Harold's count table as eset:
pdata <- as.data.frame(colData(se)[colnames(counts),])
all(rownames(pdata) == colnames(counts))
eset <- ExpressionSet(counts, AnnotatedDataFrame(pdata))
###################################################

pData(eset)$condition <- pData(eset)$strain
levels(pData(eset)$condition) <- c("A","B")
pData(eset)$batch <- factor(pData(eset)$experiment.number)
levels(pData(eset)$batch) <- 1:3

library("DESeq2")
library("edgeR")
library("limma")
source("runScripts.R")

algos <- list("DESeq2"=runDESeq2,"edgeR"=runEdgeR,"edgeRQL"=runEdgeRQL,"limma.voom"=runVoom)
namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos
resTest <- list()
resHeldout <- list()
nreps <- 30

set.seed(1) # not needed AFAIK but for safety (e.g. SAMseq)
res <- lapply(1:nreps, function(i) {   
  print(i)
  testSet <- as.character(randomSubsets[i,1:6])
  heldOutSet <- as.character(randomSubsets[i,-(1:6)])
  eTest <- eset[,testSet]
  eHeldout <- eset[,heldOutSet]
  st <- system.time({
    resTest <- as.data.frame(lapply(namesAlgos, function(n) algos[[n]](eTest)))
    resHeldout <- as.data.frame(lapply(namesAlgos, function(n) algos[[n]](eHeldout)))
  })
  print(paste(round(unname(st[3])),"seconds"))
  list(resTest=resTest,resHeldout=resHeldout)
})

resTest <- lapply(res, "[[", "resTest")
resHeldout <- lapply(res, "[[", "resHeldout")

save(resTest,resHeldout,namesAlgos,file="sensFDR_no_filter.rda")

