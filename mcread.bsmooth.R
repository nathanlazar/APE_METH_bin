mcread.bsmoothDirRaw <- function(dir, seqnames = NULL, keepCycle = FALSE, keepFilt = FALSE,
                                   header = TRUE, verbose = TRUE) {

    dir <- normalizePath(dir)
    inpattern <- "\\.cpg\\.tsv(|\\.gz)$"
    if(length(dir) != 1 || !file.info(dir)$isdir)
        stop("argument 'dir' needs to be a single directory")
    allChrFiles <- list.files(dir, pattern = inpattern, full.names = TRUE)
    if(!is.null(seqnames))
        allChrFiles <- allChrFiles[sub(inpattern, "", basename(allChrFiles)) %in% seqnames]
    if(length(allChrFiles) == 0) {
        warning(sprintf("dir '%s' is empty or has no output from 'seqnames'", dir))
        return(NULL)
    }
#####################################################################################################
# Cut out to speed up (assumes all files have the same header
#####################################################################################################
#
#    if(header) {
#        columnHeaders <- sapply(allChrFiles, function(thisfile) {
#            if(grepl("\\.gz$", thisfile))
#                con <- gzfile(thisfile)
#            else
#                con <- file(thisfile, open = "r")
#            out <- readLines(con, n = 1)
#            close(con)
#            out
#        })
#        columnHeaders <- strsplit(columnHeaders, "\t")
#        if(!all(sapply(columnHeaders, function(xx) all.equal(columnHeaders[[1]], xx))))
#            stop(sprintf("input files in dir '%s' does not have the same headers", dir))
#        columnHeaders <- columnHeaders[[1]]
#    } else
#        columnHeaders <- c("ref", "off", "strand", "Mstr", "Mcy", "Ustr", "Ucy",
#                           "filt_cycle", "filt_readlen", "filt_allele", "filt_mapq", "filt_baseq")
#####################################################################################################

    columnHeaders <- c("ref", "off", "strand", "Mstr", "Mcy", "Ustr", "Ucy",
                       "filt_cycle", "filt_readlen", "filt_allele", "filt_mapq", "filt_baseq")

    what0 <- replicate(length(columnHeaders), character(0))
    names(what0) <- columnHeaders
    int <- c(which(columnHeaders %in% c("off", "Mcy", "Ucy")), grep("^filt", columnHeaders))
    what0[int] <- replicate(length(int), integer(0))

    if(!keepCycle)
        what0[c("Mcy", "Ucy")] <- replicate(2, NULL)
    if(!keepFilt) 
        what0[grep("^filt", names(what0))] <- replicate(length(grep("^filt", names(what0))), NULL)

    outList <- mclapply(allChrFiles, function(thisfile) {
        if(verbose)
            cat(sprintf("[read.bsmoothDirRaw] Reading '%s'\n", thisfile))
        if(grepl("\\.gz$", thisfile))
            con <- gzfile(thisfile)
        else
            con <- file(thisfile)

        out <- scan(con, skip = header, what = what0, sep = "\t",
                    quote = "", na.strings = "NA", quiet = TRUE)
        close(con)
        out
    })

    listNames <- names(outList[[1]])
    names(listNames) <- listNames
    out <- lapply(listNames, function(name) {
        do.call(c, lapply(outList, function(xx) xx[[name]]))
    })
    rm(outList)
    seqn <- out[["ref"]]
    seqn[!grepl('chr', seqn)] <- paste0("chr", seqn)
    gr <- GRanges(seqnames = seqn,
                  ranges = IRanges(start = out[["off"]], width = 1))
    out[["ref"]] <- out[["off"]] <- NULL
    names(out)[names(out) == "strand"] <- "bstrand"
    out <- out[!sapply(out, is.null)]
    df <- DataFrame(out)
    elementMetadata(gr) <- df
    gr
}

mcread.bsmooth <- function(dirs, sampleNames = NULL, seqnames = NULL, returnRaw = FALSE,
                             qualityCutoff = 20, rmZeroCov = FALSE, verbose = TRUE) {
# Function adapted from BSseq package to parallelize the process of reading an 
# evidence drive using HTCondor 

    dirs <- normalizePath(dirs, mustWork = TRUE)
    if(!(all(file.info(dirs)$isdir)))
        stop("argument 'dirs' has to be directories")
    if(anyDuplicated(dirs))
        stop("duplicate entries in 'dirs'")
    if(is.null(sampleNames) && !anyDuplicated(basename(dirs)))
        sampleNames <- basename(dirs)
    if(is.null(sampleNames))
        sampleNames <- dirs
    if(!is.null(sampleNames) && (length(sampleNames) != length(dirs) || anyDuplicated(sampleNames)))
        stop("argument 'sampleNames' (if not NULL) has to have the same length as argument 'dirs', without duplicate entries")
    idxes <- seq_along(dirs)
    names(idxes) <- sampleNames
    allOut <- lapply(idxes, function(ii) {
        if(verbose) cat(sprintf("[read.bsmooth] Reading dir '%s' ... ", dirs[ii]))
        ptime1 <- proc.time()
        if(returnRaw) {
            out <- mcread.bsmoothDirRaw(dir = dirs[ii], seqnames = seqnames, keepCycle = TRUE,
                                      keepFilt = TRUE, header = TRUE, verbose = FALSE)
        } else {
            raw <- mcread.bsmoothDirRaw(dir = dirs[ii], seqnames = seqnames, keepCycle = FALSE,
                                      keepFilt = FALSE, header = TRUE, verbose = FALSE)
            out <- sampleRawToBSseq(raw, qualityCutoff = qualityCutoff, rmZeroCov, sampleName = sampleNames[ii])
        }
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) cat(sprintf("done in %.1f secs\n", stime)) 
        out
    })
    if(!returnRaw) {
        if(verbose) cat(sprintf("[read.bsmooth] Joining samples ... "))
        ptime1 <- proc.time()
        allOut <- combineList(allOut)
        ptime2 <- proc.time()
        stime <- (ptime2 - ptime1)[3]
        if(verbose) cat(sprintf("done in %.1f secs\n", stime))
    }
    allOut
}


sampleRawToBSseq <- function(gr, qualityCutoff = 20, sampleName = NULL, rmZeroCov = FALSE) {
    numberQualsGreaterThan <- function(cvec) {
        onestring <- paste(cvec, collapse = "")
        greater <- (as.integer(charToRaw(onestring)) - 33L >= qualityCutoff)
        out <- tapply(greater, rep(1:length(cvec), times = nchar(cvec)), sum)
        out
    }
    strToCov <- function(vec) {
        Cov <- rep(0, length(vec))
        wh <- which(! vec %in% c("", "0"))
        if(length(wh) > 0)
            Cov[wh] <- numberQualsGreaterThan(vec[wh])
        Cov
    }
    M <- matrix(strToCov(elementMetadata(gr)[, "Mstr"]), ncol = 1)
    Cov <- M + strToCov(elementMetadata(gr)[, "Ustr"])
    elementMetadata(gr) <- NULL
    BSseq(gr = gr, M = M, Cov = Cov, sampleNames = sampleName, rmZeroCov = rmZeroCov)
}
