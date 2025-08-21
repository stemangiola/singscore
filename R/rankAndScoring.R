#' @include singscore.R
#' @importFrom DelayedMatrixStats colRanks
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedArray blockApply
#' @importFrom DelayedArray colAutoGrid
#' @importFrom DelayedArray rowAutoGrid
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom HDF5Array HDF5Array
#' @importFrom BiocParallel MulticoreParam
NULL

rankExpr <- function(exprsM, tiesMethod = "min") {
  rname= rownames(exprsM)
  cname = colnames(exprsM)

  # For DelayedMatrix, use our custom block-wise ranking
  if (is(exprsM, "DelayedArray")) {
    message("Using block-wise ranking for DelayedMatrix")
    rankedData <- rankExprDelayed(exprsM, tiesMethod, workers = 1)
    return(rankedData)
  }
  
  # For regular matrices, use colRanks
  rankedData = colRanks(
    exprsM,
    ties.method = tiesMethod,
    preserveShape = TRUE
  )
  
  rownames(rankedData) = rname
  colnames(rankedData) = cname
  
  #indicator of the type of ranks
  attr(rankedData, 'stable') = FALSE
  return (rankedData)
}

#' DelayedArray-compatible ranking function with parallel processing
#'
#' This function ranks gene expression data and returns a DelayedArray that remains
#' on-disk, preventing memory issues with large datasets. It supports parallel processing
#' and uses block-wise operations to handle large DelayedMatrix objects efficiently.
#'
#' @param exprsM A DelayedArray or matrix of expression values
#' @param tiesMethod Method for handling ties
#' @param workers Number of parallel workers (default: 1)
#' @param output_file Optional HDF5 file path for output storage
#'
#' @return A DelayedArray with ranked gene expression data
#' @export
rankExprDelayed <- function(exprsM, tiesMethod = "min", workers = 1, output_file = NULL) {
  rname = rownames(exprsM)
  cname = colnames(exprsM)

  # Convert to DelayedArray if not already
  if (!is(exprsM, "DelayedArray")) {
    exprsM <- DelayedArray(exprsM)
  }

  # Use block-wise ranking to avoid colRanks bugs
  message("Using block-wise ranking for DelayedMatrix")
  
  # Get grid for column-wise processing with smaller blocks for large datasets
  if (ncol(exprsM) > 10000) {
    # For very large datasets, use smaller blocks
    col_grid <- DelayedArray::colAutoGrid(exprsM, ncol = min(1000, ncol(exprsM) %/% 10))
    message("Using smaller blocks for large dataset (", ncol(exprsM), " columns)")
  } else {
    col_grid <- DelayedArray::colAutoGrid(exprsM)
  }
  
  # Configure parallel processing with progress bar
  # Reduce workers for very large datasets to prevent memory issues
  # if (ncol(exprsM) > 50000) {
  #   actual_workers <- min(workers, 2)
  #   message("Reducing workers to ", actual_workers, " for very large dataset")
  # } else {
  #   actual_workers <- workers
  # }
  
  if (workers > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
    message("Using ", workers, " workers for parallel processing with progress bar")
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = TRUE)
    message("Using single worker with progress bar")
  }
  
  # Process each column block
  ranked_blocks <- DelayedArray::blockApply(
    exprsM,
    FUN = function(block) {
      # Use colRanks on each block
      colRanks(
        block,
        ties.method = tiesMethod,
        preserveShape = TRUE
      )
    },
    grid = col_grid,
    BPPARAM = bpparam
  )
  
  # Combine blocks into a single DelayedArray
  # Use cbind instead of arbind to avoid vector size issues
  rankedData <- do.call(cbind, ranked_blocks)
  
  # Ensure it's a DelayedArray
  if (!is(rankedData, "DelayedArray")) {
    rankedData <- DelayedArray(rankedData)
  }
  
  # Save to HDF5 if output file specified
  if (!is.null(output_file)) {
    rankedData <- writeHDF5Array(rankedData, filepath = output_file, name = "ranked_data")
  }
  
  #indicator of the type of ranks
  attr(rankedData, 'stable') = FALSE
  return(rankedData)
}

#' Block-wise stable gene ranking helper
#'
#' Computes stable gene-based ranks for a block of expression data. For each column (sample),
#' each gene is ranked according to its position relative to the sorted values of the stable genes.
#' Used internally by \code{rankExprStable_delayed} to enable block processing of large matrices.
#'
#' @param block A numeric matrix (genes x samples) block of expression values.
#' @param st_idx Integer vector of row indices corresponding to stable genes.
#'
#' @return A numeric matrix of the same dimensions as \code{block} with stable gene-based ranks.
#'
#' @importFrom stats findInterval
#' @keywords internal
#' @noRd
.rank_block_by_stable_genes <- function(block, st_idx) {
  
  nb <- nrow(block)
  kb <- ncol(block)
  out <- matrix(NA_real_, nrow = nb, ncol = kb)
  s_idx <- st_idx
  s_idx <- s_idx[s_idx >= 1 & s_idx <= nb]
  for (j in seq_len(kb)) {
    x <- block[, j]
    # compute ranks: number of stable genes with value less than x, plus 1
    rs <- rowSums(outer(x, x[s_idx], '>')) + 1
    out[, j] <- rs
  }
  out
}

rankExprStable <- function(exprsM, tiesMethod = "min", stgenes) {
  stgenes = intersect(stgenes, rownames(exprsM))
  stopifnot(length(stgenes) > 0)
  
  rname = rownames(exprsM)
  cname = colnames(exprsM)
  
  #compute ranks
  rankedData = apply(exprsM, 2, function(x) {
    rowSums(outer(x, x[stgenes], '>')) + 1
  })
  
  #normlise ranks
  rankedData = rankedData / (length(stgenes) + 1)
  
  rownames(rankedData) = rname
  colnames(rankedData) = cname
  
  #indicator of the type of ranks
  attr(rankedData, 'stable') = TRUE
  return(rankedData)
}

#' @importFrom DelayedArray blockApply
#' @importFrom DelayedArray colAutoGrid
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel MulticoreParam
rankExprStable_delayed <- function(exprsM, tiesMethod = "min", stgenes, 
                                   workers = 1) {
 
  stgenes = intersect(stgenes, rownames(exprsM))
  stopifnot(length(stgenes) > 0)

  rname = rownames(exprsM)
  cname = colnames(exprsM)

  # Work with matrix-like via DelayedArray, block over columns
  dx <- DelayedArray(exprsM)
  st_idx <- match(stgenes, rname)
  
  # Set up BiocParallel param with progress bar enabled
  param <- MulticoreParam(workers, progressbar = TRUE)
  
  parts <- blockApply(
    dx,
    FUN = .rank_block_by_stable_genes,
    grid = colAutoGrid(dx),
    st_idx = st_idx,
    BPPARAM = param
  )
  rankedData <- do.call(cbind, parts)
  rankedData = rankedData / (length(stgenes) + 1)
  rownames(rankedData) = rname
  colnames(rankedData) = cname

  #indicator of the type of ranks
  attr(rankedData, 'stable') = TRUE
  return(rankedData)
}

#define the main helper functions
#helper: check if all genes in the gene set are in the data
checkGenes <- function (geneset, background) {
  #check if there are some missing genes in the geneset
  missingGenes = setdiff(geneIds(geneset), background)
  if (length(missingGenes) > 0) {
    warningMsg = paste(length(missingGenes), "genes missing:", sep = ' ')
    warningMsg = paste(warningMsg, paste(missingGenes, collapse = ', '), sep = ' ')
    warning(warningMsg)
  }
  geneIds(geneset) = setdiff(geneIds(geneset), missingGenes)

  return(geneset)
}

checkGenesMulti <- function(geneset_colc, background) {
  geneset_colc = endoapply(geneset_colc, function(geneset) {
    geneIds(geneset) = intersect(geneIds(geneset), background)
    return(geneset)
  })

  #remove genesets with no genes
  geneset_colc = geneset_colc[sapply(geneset_colc, function(x) length(geneIds(x))) > 0]

  return(geneset_colc)
}

#helper: compute the theoretical boundaries of scores when direction is KNOWN.
# Should be used after filtering out missing genes
calcBounds <- function(gsSize, bgSize, stableSc = FALSE) {
  if (stableSc) {
    lowBound = 0
    upBound = 1
  } else {
    lowBound = (gsSize + 1) / 2
    upBound = (2 * bgSize - gsSize + 1) / 2
  }

  return(list('lowBound' = lowBound, 'upBound' = upBound))
}

calcBoundsMulti <- function(gsSizes, bgSize, stableSc = FALSE) {
  bounds = lapply(gsSizes, calcBounds, bgSize, stableSc)

  return(bounds)
}

#helper: compute the theoretical boundaries of scores when direction is UNKNOWN.
# Should be used after filtering out missing genes
calcBoundsUnknownDir <- function(gsSize, bgSize, stableSc = FALSE) {
  gsSize = floor(gsSize / 2) #number of unique ranks in the geneset
  bgSize = ceiling(bgSize / 2) #number of unique ranks in the matrix

  return(calcBounds(gsSize, bgSize, stableSc))
}

calcBoundsUnknownDirMulti <- function(gsSizes, bgSize, stableSc = FALSE) {
  gsSizes = floor(gsSizes / 2) #number of unique ranks in the geneset
  bgSize = ceiling(bgSize / 2) #number of unique ranks in the matrix

  return(calcBoundsMulti(gsSizes, bgSize, stableSc))
}

#helper: empty data frame of results
emptyScoreDf <- function(onegs = TRUE) {
  if (onegs) {
    df = data.frame("TotalScore" = numeric(),
                    "TotalDispersion" = numeric())
  } else{
    df = data.frame(
      "TotalScore" = numeric(),
      "TotalDispersion" = numeric(),
      "UpScore" = numeric(),
      "UpDispersion" = numeric(),
      "DownScore" = numeric(),
      "DownDispersion" = numeric()
    )
  }

  return(df)
}

#helper: compute scores for one gene set given a rank matrix and bounds
# Assume genes are filtered
calcScores <- function(ranks, geneset, dispersionFun, bounds, scoffset = 0) {
  #compute raw score
  gsranks = ranks[geneIds(geneset), , drop = FALSE]
  
  # Ensure gsranks is a proper matrix for colMeans
  if (is(gsranks, "DelayedArray")) {
    gsranks <- as.matrix(gsranks)
  }
  
  gsscore = colMeans(gsranks)
  gsscore = (gsscore - bounds$lowBound) / (bounds$upBound - bounds$lowBound)

  #calculate dispersion
  dispersion = apply(gsranks, 2, dispersionFun)

  return(list('TotalScore' = gsscore + scoffset, 'TotalDispersion' = dispersion))
}

#helper: compute scores for two gene sets given a rank matrix and bounds
# Assume genes are filtered
calcScoresUpDn <- function(ranks, upset, downset, dispersionFun, upbounds, downbounds, scoffset = 0) {
  #compute independent scores
  upscores = calcScores(ranks, upset, dispersionFun, upbounds, scoffset)
  dnscores = calcScores(ranks, downset, dispersionFun, downbounds, 0)
  dnscores$TotalScore = 1 - dnscores$TotalScore + scoffset

  #compute total scores
  totalscores = list(
    'TotalScore' = upscores$TotalScore + dnscores$TotalScore,
    'TotalDispersion' = (upscores$TotalDispersion + dnscores$TotalDispersion) / 2,
    'UpScore' = upscores$TotalScore,
    'UpDispersion' = upscores$TotalDispersion,
    'DownScore' = dnscores$TotalScore,
    'DownDispersion' = dnscores$TotalDispersion
  )

  return(totalscores)
}

#singscore function for a single geneset/signature
singleSingscore <-
  function (rankData,
            upSet,
            downSet = NULL,
            subSamples = NULL,
            centerScore = TRUE,
            dispersionFun = mad,
            knownDirection = TRUE) {
    #0. determine type of singscore being used
    stableSc = attr(rankData, 'stable')

    #1. subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = FALSE]
    }

    #2. center scores?
    center_offset = ifelse(centerScore, -0.5, 0)

    #3.1 filter genes - return empty data.frame if no scores
    upSet = checkGenes(upSet, rownames(rankData))
    if (length(geneIds(upSet)) == 0)
      return(emptyScoreDf(onegs = is.null(downSet)))

    #known direction
    if (!knownDirection) {
      #2.3 do not center scores
      warning('\'centerScore\' is disabled for this setting')
      center_offset = 0

      #modify rank matrix
      rankData = apply(rankData, 2, function(x) {
        x = abs(x - ceiling(median(x)))
        return(x)
      })

      #4.3 calculate bounds
      upset_bounds = calcBoundsUnknownDir(length(geneIds(upSet)), nrow(rankData), stableSc)
    } else {
      #4.1 calculate bounds
      upset_bounds = calcBounds(length(geneIds(upSet)), nrow(rankData), stableSc)
    }

    #for two gene sets
    if (!is.null(downSet)) {
      #3.2 filter genes - return empty data.frame if no scores
      downSet = checkGenes(downSet, rownames(rankData))
      if (length(geneIds(downSet)) == 0)
        return(emptyScoreDf(onegs = FALSE))

      #4.2 calculate bounds
      downset_bounds = calcBounds(length(geneIds(downSet)), nrow(rankData), stableSc)

      #5.2 compute scores
      scores = calcScoresUpDn(
        rankData,
        upSet,
        downSet,
        dispersionFun,
        upset_bounds,
        downset_bounds,
        center_offset
      )
    } else{
      #5.1 & 5.3 compute scores
      scores = calcScores(rankData,
                          upSet,
                          dispersionFun,
                          upset_bounds,
                          center_offset)
    }

    #6 process the results
    scores = data.frame(scores)
    rownames(scores) = colnames(rankData)

    return(scores)
  }

#singscore function for a multiple geneset/signature
multiSingscore <-
  function (rankData,
            upSetColc,
            downSetColc = NULL,
            subSamples = NULL,
            centerScore = TRUE,
            dispersionFun = mad,
            knownDirection = TRUE) {
    #0. determine type of singscore being used
    stableSc = attr(rankData, 'stable')
    
    #1. subset the data for samples whose calculation is to be performed
    if (!is.null(subSamples)) {
      rankData <- rankData[, subSamples, drop = FALSE]
    }

    #2. center scores?
    center_offset = ifelse(centerScore, -0.5, 0)

    #3.1 filter genes - empty genesets will be removed
    upNames = names(upSetColc)
    upSetColc = checkGenesMulti(upSetColc, rownames(rankData))
    if (length(upSetColc) == 0)
      return(NULL)

    #known direction
    upSizes = sapply(upSetColc, function(x) length(geneIds(x)))
    if (!knownDirection) {
      #2.3 do not center scores
      warning('\'centerScore\' is disabled for this setting')
      center_offset = 0

      #modify rank matrix
      rankData = apply(rankData, 2, function(x) {
        x = abs(x - ceiling(median(x)))
        return(x)
      })

      #4.3 calculate bounds
      upset_bounds = calcBoundsUnknownDirMulti(upSizes, nrow(rankData), stableSc)
    } else {
      #4.1 calculate bounds
      upset_bounds = calcBoundsMulti(upSizes, nrow(rankData), stableSc)
    }

    #for two gene sets
    if (!is.null(downSetColc)) {
      #3.2 filter genes - empty genesets will be removed
      downNames = names(downSetColc)
      stopifnot(all(upNames == downNames))
      
      downSetColc = checkGenesMulti(downSetColc, rownames(rankData))
      if (length(downSetColc) == 0)
        return(NULL)

      #drop those not matching with upSets
      commonNames = intersect(names(upSetColc), names(downSetColc))
      downSetColc = downSetColc[names(downSetColc) %in% commonNames]

      #4.2 calculate bounds
      downSizes = sapply(downSetColc, function(x) length(geneIds(x)))
      downset_bounds = calcBoundsMulti(downSizes, nrow(rankData), stableSc)

      #check that names match
      upselect = names(upSetColc) %in% commonNames
      upSetColc = upSetColc[upselect]
      upset_bounds = upset_bounds[upselect]

      #5.2 compute scores
      scores = mapply(
        calcScoresUpDn,
        list(rankData),
        upSetColc,
        downSetColc,
        list(dispersionFun),
        upset_bounds,
        downset_bounds,
        center_offset
      )
    } else{
      #5.1 & 5.3 compute scores
      scores = mapply(
        calcScores,
        list(rankData),
        upSetColc,
        list(dispersionFun),
        upset_bounds,
        center_offset
      )
    }

    #6 process the results
    dispersions = scores['TotalDispersion', ]
    dispersions = mapply(c, dispersions)
    scores = scores['TotalScore', ]
    scores = mapply(c, scores)
    
    #create the correct return structure
    if (is.matrix(scores) & is.matrix(dispersions)) {
      scores = t(scores)
      dispersions = t(dispersions)
    } else {
      scores = matrix(scores, ncol = 1)
      dispersions = matrix(dispersions, ncol = 1)
    }
    rownames(scores) = rownames(dispersions) = names(upSetColc)
    colnames(scores) = colnames(dispersions) = colnames(rankData)
    
    #create NA entries for empty gene sets
    if (nrow(scores) < length(upNames)) {
      #create NA matrix
      naNames = setdiff(upNames, rownames(scores))
      namat = matrix(NA, nrow = length(naNames), ncol = ncol(scores))
      rownames(namat) = naNames
      
      #combine to score and dispersion matrices
      scores = rbind(scores, namat)[upNames, ]
      dispersions = rbind(dispersions, namat)[upNames, ]
    }

    return(list('Scores' = scores, 'Dispersions' = dispersions))
  }
