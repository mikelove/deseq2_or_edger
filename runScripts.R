runDESeq2 <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  # p-value adjustment in DESeq2 does not occur
  # over rows with all zero count, so this line is not needed:
  ## padj <- p.adjust(pvals,method="BH")
  # DESeq2 filters low count genes using genefilter.
  # set the remaining NA (filtered) adjusted p-values to 1:
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

runEdgeR <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj,
       beta=log2(exp(1)) * edger.fit$coefficients[,"conditionB"])
}

runEdgeRQL <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmQLFit(dgel,design)
  edger.qlf <- glmQLFTest(edger.fit)
  pvals <- edger.qlf$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj,
       beta=log2(exp(1)) * edger.fit$coefficients[,"conditionB"])
}

runVoom <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=tt$logFC)
}
