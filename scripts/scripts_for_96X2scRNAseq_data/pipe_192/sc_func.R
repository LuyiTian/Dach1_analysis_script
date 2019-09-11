suppressMessages(library(topGO))
suppressMessages(library(edgeR))

#' perform GO enrichment analysis on top ranked loading genes for the top .
#'
#'
#' @param pca_out class \code{"prcomp"} containing a the matrix of variable loadings.
#' @param top_gene the number of top ranked loading genes used for GO enrichment analysis
#' @param k the number of principle component to analyse
#' @param p.val the the p-value threshold of GO enrichment test
#' @return a list of tables that gives the significant GO terms
#' 
PCA_GO_enrichment = function(pca_out, 
                             top_gene = 100,
                             k = 15,
                             p.val = 0.01){
  result = list()
  for (i in 1:k){
    all_gene = rep(0,nrow(pca_out$rotation))
    names(all_gene) = rownames(pca_out$rotation)
    all_gene[names(all_gene) %in% names(pca_out$rotation[,i][order(abs(pca_out$rotation[,i]),decreasing = TRUE)][1:50])] = 1
    sampleGOdata <- new("topGOdata",
                        description = "test", ontology = "BP",
                        allGenes = all_gene, geneSelectionFun = function(a_g){a_g == 1},
                        nodeSize = 10,
                        annot=annFUN.org, mapping="org.Mm.eg.db", ID = "SYMBOL")
    
    res = runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    tabl = GenTable(sampleGOdata, classic = res,
                    orderBy = "weight", ranksOf = "classic", topNodes = 10)
    tabl$classic = as.numeric(tabl$classic)
    if (sum(tabl$classic < p.val) == 0)
    {
      result[[i]] = "NA"
    }
    else
    {
      result[[i]] = tabl[tabl$classic<p.val,]
    }
  }
  result
}

Diff_GO_enrichment = function(diff,
                             gene_expr,
                             k = ncol(diff@eigenvectors),
                             p.val = 0.001){
  result = list()
  result1 = list()
  for (i in 1:k)
  {
    design = model.matrix(~ dif@eigenvectors[,i])
    DGE_obj = DGEList(counts=gene_expr)
    DGE_obj = estimateDisp(DGE_obj, design=design, robust=TRUE)
    
    fit = glmQLFit(DGE_obj,design=design)
    lrt = glmQLFTest(fit, coef=2)
    
    tmp_tag = topTags(lrt, n=Inf, p.value=p.val)
    result1[[i]] = tmp_tag
    
    all_gene = rep(0,nrow(gene_expr))
    names(all_gene) = rownames(gene_expr)
    if (nrow(tmp_tag@.Data[[1]])>100){
      all_gene[names(all_gene) %in% rownames(tmp_tag@.Data[[1]])[1:100]] = 1
    }
    else
    {
      all_gene[names(all_gene) %in% rownames(tmp_tag@.Data[[1]])] = 1
    }
    sampleGOdata <- new("topGOdata",
                        description = "test", ontology = "BP",
                        allGenes = all_gene, geneSelectionFun = function(a_g){a_g == 1},
                        nodeSize = 10,
                        annot=annFUN.org, mapping="org.Mm.eg.db", ID = "SYMBOL")
    
    res = runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    tabl = GenTable(sampleGOdata, classic = res,
                    orderBy = "weight", ranksOf = "classic", topNodes = 10)
    tabl$classic = as.numeric(tabl$classic)
    if (sum(tabl$classic < p.val) == 0)
    {
      result[[i]] = "NA"
    }
    else
    {
      result[[i]] = tabl[tabl$classic<p.val,]
    }
  }
  combied = list(GO_enrich=result, DE_genes=result1)
}