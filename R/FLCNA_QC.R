#' @title FLCNA_QC
#'
#' @description Perform QC step on single cells and bins.
#'
#'
#' @param Y_raw raw read count matrix.
#' @param ref_raw raw GRanges object with corresponding GC content and mappability for quality control.
#' @param cov_thresh scalar variable specifying the lower bound of read count summation of each cell. Default is \code{0}.
#' @param minCountQC the minimum read coverage required for normalization and EM fitting. Defalut is \code{20}.
#' @param mapp_thresh scalar variable specifying mappability of each genomic bin. Default is \code{0.9}.
#' @param gc_thresh vector specifying the lower and upper bound of GC content threshold for quality control. Default is \code{20-80}.
#'
#' @return A list with components after quality control.
#' \item{Y}{Read depth matrix after quality control.}
#' \item{ref}{A GRanges object specifying whole genomic bin positions after quality control.}
#'
#'
#' @import stats
#' @export
FLCNA_QC <- function(Y_raw, 
                     ref_raw,
                     cov_thresh = 0, 
                     minCountQC = 20, 
                     mapp_thresh = 0.9,
                     gc_thresh = c(20, 80)) {
  

  if (length(ref_raw) != nrow(Y_raw)) {
    stop("Invalid inputs: length of ref and # of rows
            in read count matrix must be the same")
  }
  
  mapp <- ref_raw$mapp
  gc <- ref_raw$gc
  
  ##sample based QC
  sampfilter1 <- (apply(Y_raw, 2, sum) <= cov_thresh)
  message("Removed ", sum(sampfilter1),
          " samples due to failed library preparation.")
  sampfilter2 <- (apply(Y_raw, 2, mean) <= minCountQC)
  message("Removed ", sum(sampfilter2),
          " samples due to failure to meet min coverage requirement.")    
  if (sum(sampfilter1 | sampfilter2) != 0) {
    Y <- Y_raw[, !(sampfilter1 | sampfilter2)]
  } else {
    Y <- Y_raw
  }
  
  ##Bin based QC
  binfilter1 <- (gc < gc_thresh[1] | gc > gc_thresh[2])
  message("Excluded ", sum(binfilter1),
          " bins due to extreme GC content.")
  binfilter2 <- (mapp < mapp_thresh)
  message("Excluded ", sum(binfilter2),
          " bins due to low mappability.")
  if (sum(binfilter1 | binfilter2) != 0) {
    ref <- ref_raw[!(binfilter1 | binfilter2)]
    Y <- Y[!(binfilter1 | binfilter2), ]
  } else {
    ref <- ref_raw
    Y <- Y
  }

  message("There are ", ncol(Y), " samples and ",
          nrow(Y), " bins after QC step. ")
  list(Y = Y, ref = ref)
}