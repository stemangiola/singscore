#' @include singscore.R rankAndScoring.R
NULL

#'@title  Rank genes by the gene expression intensities
#'@description The \code{rankGenes} function is a generic function that can deal
#'  with mutilple types of inputs. Given a matrix of gene expression that has
#'  samples in columns, genes in rows, and values being gene expression
#'  intensity,\code{rankGenes} ranks gene expression intensities in each sample.
#'
#'  It can also work with S4 objects that have gene expression matrix as a
#'  component (i.e ExpressionSet, DGEList,SummarizedExperiment). It calls the
#'  \code{rank} function in the base package which ranks the gene expression
#'  matrix by its absolute expression level. If the input is S4 object of
#'  \code{DGEList, ExpressionSet, or SummarizedExperiment}, it will extract the
#'  gene expression matrix from the object and rank the genes. The default
#'  'tiesMethod' is set to 'min'.
#'
#'@param expreMatrix matrix, data.frame, ExpressionSet, DGEList or
#'  SummarizedExperiment storing gene expression measurements
#'@param tiesMethod character, indicating what method to use when dealing with
#'  ties
#'@param stableGenes character, containing a list of stable genes to be used to
#'  rank genes using expression of stable genes. This is required when using the
#'  stable genes dependent version of singscore (see details in `simpleScore`).
#'  Stable genes for solid cancers (carcinomas) and blood transcriptomes can be
#'  obtained using the `getStableGenes` function
#'@param workers integer, number of workers for parallel processing when using
#'  stable genes with DelayedMatrix objects (default: 1)
#'
#'@name rankGenes
#'
#'@return The ranked gene expression matrix that has samples in columns and
#'  genes in rows. Unit normalised ranks are returned if data is ranked using
#'  stable genes
#'@examples
#' rankGenes(toy_expr_se) # toy_expr_se is a gene expression dataset
#'
#' # ExpressionSet object
#' emat <- SummarizedExperiment::assay(toy_expr_se)
#' e <- Biobase::ExpressionSet(assayData = as.matrix(emat))
#' rankGenes(e)
#'
#' #scoring using the stable version of singscore
#' rankGenes(e, stableGenes = c('2', '20', '25'))
#'
#' # With progress bar for DelayedMatrix objects
#' \dontrun{
#' #for real cancer or blood datasets, use getStableGenes()
#' rankGenes(cancer_expr, stableGenes = getStableGenes(5))
#' rankGenes(blood_expr, stableGenes = getStableGenes(5, type = 'blood'))
#' 
#' # Enable progress bar for large DelayedMatrix objects
#' rankGenes(delayed_matrix, stableGenes = stable_genes, 
#'           workers = 4)
#' }
#'@seealso \code{\link{getStableGenes}}, \code{\link{simpleScore}},
#'  \code{\link{rank}}, ["ExpressionSet"][ExpressionSet-class],
#'  ["SummarizedExperiment"][SummarizedExperiment-class],
#'  ["DGEList"][DGEList-class]
#'
#'
#'@export
setGeneric("rankGenes",
           function(expreMatrix, 
                    tiesMethod = 'min',
                    stableGenes = NULL,
                    workers = 1) standardGeneric("rankGenes"))

#' @rdname rankGenes
setMethod("rankGenes",
          signature('matrix','ANY', 'ANY'), 
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL, 
                   workers = 1){
            stopifnot(tiesMethod %in% c("max", "average","min"))
            if (is.null(stableGenes)) {
              rankMat = rankExpr(expreMatrix, tiesMethod)
            } else {
              stopifnot(is.character(stableGenes))
              stopifnot(length(stableGenes) > 0)
              rankMat = rankExprStable_delayed(expreMatrix, tiesMethod, stableGenes, 
                                               workers)
            }
            return(rankMat)
})

#' @rdname rankGenes
setMethod("rankGenes", 
          signature('data.frame', 'ANY', 'ANY'), 
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL,
                   workers = 1){
            rankGenes(as.matrix(expreMatrix), tiesMethod, stableGenes, workers)
})

#' @rdname rankGenes
setMethod("rankGenes", 
          signature('DGEList','ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL,
                   workers = 1){
            rankGenes(expreMatrix$counts, tiesMethod, stableGenes, workers)
})
#' @rdname rankGenes
setMethod("rankGenes", 
          signature('ExpressionSet', 'ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL,
                   workers = 1){
            rankGenes(Biobase::exprs(expreMatrix), tiesMethod, stableGenes, workers)
})

#' @rdname rankGenes
setMethod("rankGenes",
          signature('SummarizedExperiment', 'ANY', 'ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL,
                   workers = 1){
            rankGenes(SummarizedExperiment::assay(expreMatrix), tiesMethod, stableGenes, workers)
          })

#' @rdname rankGenes
setMethod("rankGenes",
          signature('DelayedMatrix','ANY','ANY'),
          function(expreMatrix, tiesMethod = 'min', stableGenes = NULL,
                   workers = 1){
            if (is.null(stableGenes)) {
              return(rankExpr(expreMatrix, tiesMethod))
            } else {
              return(rankExprStable_delayed(expreMatrix, tiesMethod, stableGenes, 
                                           workers))
            }
          })
