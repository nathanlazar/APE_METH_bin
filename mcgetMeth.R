mcgetMeth <- function(BSseq, regions = NULL, type = c("smooth", "raw"),
                      what = c("perBase", "perRegion"), confint = FALSE, alpha = 0.95) {
    p.conf <- function(p, n, alpha) {
        z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
        upper <- (p + z^2/(2*n) + z*sqrt(  (p*(1-p) + z^2/(4*n)) / n)) /
            (1+z^2/n)
        lower <- (p + z^2/(2*n) - z*sqrt(  (p*(1-p) + z^2/(4*n)) / n)) /
            (1+z^2/n)
        return(list(meth = p, lower = lower, upper = upper))
    }
    stopifnot(is(BSseq, "BSseq"))
    type <- match.arg(type)
    if(type == "smooth" & !hasBeenSmoothed(BSseq))
        stop("'type=smooth' requires the object to have been smoothed.")
    what <- match.arg(what)
    z <- abs(qnorm((1 - alpha)/2, mean = 0, sd = 1))
    if(is.null(regions) && type == "smooth") {
        meth <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef"))
        if(confint) {
            upper <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef") +
                                      z * getBSseq(BSseq, type = "se.coef"))
            lower <- getBSseq(BSseq, type = "trans")(getBSseq(BSseq, type = "coef") -
                                      z * getBSseq(BSseq, type = "se.coef"))
            return(list(meth = meth, lower = lower, upper = upper))
        } else {
            return(meth)
        }
    }
    if(is.null(regions) && type == "raw") {
        meth <- getBSseq(BSseq, type = "M") / getBSseq(BSseq, type = "Cov")
        if(confint) {
            return(p.conf(meth, n = getBSseq(BSseq, type = "Cov"), alpha = alpha))
        } else {
            return(meth)
        }
    }
    ## At this point, regions have been specified
    if(class(regions) == "data.frame")
        regions <- data.frame2GRanges(regions)
    stopifnot(is(regions, "GenomicRanges"))
    if(confint) stop("'confint = TRUE' is not supported by 'getMeth' when regions is given")
    grBSseq <- granges(BSseq)
    mm <- as.matrix(findOverlaps(regions, grBSseq))
    mmsplit <- split(mm[,2], mm[,1])
    if(what == "perBase") {
        if(type == "smooth") {
            out <- mclapply(mmsplit, function(xx) {
                getBSseq(BSseq, "trans")(getBSseq(BSseq, "coef")[xx,,drop = FALSE])
            })
        }
        if(type == "raw") {
            out <- mclapply(mmsplit, function(xx) {
                getBSseq(BSseq, "M")[xx,,drop = FALSE] / getBSseq(BSseq, "Cov")[xx,,drop = FALSE]
            })
        }
        outList <- vector("list", length(regions))
        outList[as.integer(names(mmsplit))] <- out
        return(outList)
    }
    if(what == "perRegion") {
        if(type == "smooth") {
            out <- mclapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "trans")(getBSseq(BSseq, "coef")[xx,,drop = FALSE]), na.rm = TRUE)
            })
        }
        if(type == "raw") {
            out <- mclapply(mmsplit, function(xx) {
                colMeans(getBSseq(BSseq, "M")[xx,,drop = FALSE] / getBSseq(BSseq, "Cov")[xx,,drop = FALSE], na.rm = TRUE)
            })
        }
        out <- do.call(rbind, out)
        outMatrix <- matrix(NA, ncol = ncol(BSseq), nrow = length(regions))
        colnames(outMatrix) <- sampleNames(BSseq)
        outMatrix[as.integer(rownames(out)),] <- out
        return(outMatrix)
    }
}
