mcgetCoverage <- function(BSseq, regions = NULL, type = c("Cov", "M"),
                    what = c("perBase", "perRegionAverage", "perRegionTotal")) {
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    what <- match.arg(what)
    if(is.null(regions)) {
        if(what == "perBase")
            return(getBSseq(BSseq, type = type))
        if(what == "perRegionTotal")
            return(colSums(getBSseq(BSseq, type = type)))
        if(what == "perRegionAverage")
            return(colMeans(getBSseq(BSseq, type = type)))
    }
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(is(regions, "GenomicRanges"))
    grBSseq <- granges(BSseq)
    mm <- as.matrix(findOverlaps(regions, grBSseq))
    mmsplit <- split(mm[,2], mm[,1])
    if(what == "perBase") {
        if(type == "Cov") {
            out <- mclapply(mmsplit, function(xx) {
                getBSseq(BSseq, "Cov")[xx,,drop = FALSE]
            })
        }
        if(type == "M") {
            out <- mclapply(mmsplit, function(xx) {
                getBSseq(BSseq, "M")[xx,,drop = FALSE]
            })
        }
        outList <- vector("list", length(regions))
        outList[as.integer(names(mmsplit))] <- out
        return(outList)
    }
    if(what == "perRegionAverage") {
        if(type == "Cov") {
            out <- mclapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        if(type == "M") {
            out <- mclapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "M")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
    }
    if(what == "perRegionTotal") {
        if(type == "Cov") {
            out <- mclapply(mmsplit, function(xx) {
                colSums(getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        if(type == "M") {
            out <- mclapply(mmsplit, function(xx) {
                colSums(getBSseq(BSseq, "M")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
    }
    out <- do.call(rbind, out)
    outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
    colnames(outMatrix) <- sampleNames(BSseq)
    outMatrix[as.integer(rownames(out)),] <- out
    return(outMatrix)
}