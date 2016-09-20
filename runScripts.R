runDESeq2 <- function(e) {
  padj <- rep(NA, nrow(e))
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
  keep <- cpmFilter(e, cpm=FALSE)
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
  keep <- cpmFilter(e, cpm=FALSE)
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
  keep <- cpmFilter(e, cpm=FALSE)
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
    dgel <- DGEList(exprs(e))
    # this should actually be min(table(e$condition))
    # instead of 2 according to edgeR docs...
    return(rowSums(cpm(dgel) > 1) >= 2)
  } else {
    return(rowSums(exprs(e)) > 0)
  }
}
