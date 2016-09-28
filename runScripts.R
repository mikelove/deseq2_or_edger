runDESeq2 <- function(e) {
  padj <- rep(NA, nrow(e))
  # cpm=FALSE means don't do any filtering
  keep <- cpmFilter(e, cpm=FALSE)
  e <- e[keep,]

  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)

  # slightly different bc DESeq2 filters internally
  padj[keep] <- res$padj
  padj[is.na(padj)] <- 1
  padj
}

runEdgeR <- function(e) {
  padj <- rep(NA, nrow(e))
  keep <- cpmFilter(e, cpm=TRUE)
  e <- e[keep,]

  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  pvals <- edger.lrt$table$PValue

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  padj
}

runEdgeRQL <- function(e) {
  padj <- rep(NA, nrow(e))
  keep <- cpmFilter(e, cpm=TRUE)
  e <- e[keep,]
  
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # current recommendation (Sep 2016) according to vignette:
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmQLFit(dgel,design)
  edger.qlf <- glmQLFTest(edger.fit)
  pvals <- edger.qlf$table$PValue

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  padj
}

runVoom <- function(e) {
  padj <- rep(NA, nrow(e))
  keep <- cpmFilter(e, cpm=TRUE)
  e <- e[keep,]
  
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 

  padj[keep] <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  padj
}

cpmFilter <- function(e, cpm=TRUE) {
  if (cpm) {
    # this is the filtering recommendation from Gordon Smyth
    # https://support.bioconductor.org/p/85511/#85514
    L <- min(colSums(exprs(e))/1e6)
    dgel <- DGEList(exprs(e))
    # here I use 3 even for the heldout set (n=7/8),
    # because otherwise the filtering is too strict
    # to the detriment of edgeR and limma-voom performance
    return(rowSums(cpm(dgel) > 10/L) >= 3)
  } else {
    return(rowSums(exprs(e)) > 0)
  }
}
