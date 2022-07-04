EMnormalize_custom = function (object, G = 2:6, priorScale = 5, peakThresh = 0.5, mergeVal = 0.1, Title = NA, verbose = TRUE) {
    if (!rCGH:::.validrCGHObject(object)) 
        return(NULL)
    op <- options()
    options(warn = -1)
    ploidy <- as.numeric(getInfo(object, "ploidy"))
    segTable <- getSegTable(object)
    if (nrow(segTable) == 0) {
        stop("Please run the segmentation step before centralizing.")
    }
    simulLR <- rCGH:::.simulateLRfromST(segTable)
    EM <- Mclust(simulLR, G = G, prior = priorControl(scale = priorScale))
    nG <- EM$G
    m <- EM$parameters$mean
    p <- EM$parameters$pro
    s <- EM$parameters$variance$sigmasq
    if (length(s) == 1) 
        s <- rep(s, length(m))
    ord <- order(m)
    m <- m[ord]
    p <- p[ord]
    s <- s[ord]
    if (mergeVal > 0) {
        if (verbose) 
            message("Merging peaks closer than ", mergeVal, 
                " ...")
        mergedPars <- rCGH:::.mergePeaks(nG, simulLR, m, s, p, mergeVal, 
            verbose)
        m <- mergedPars$m
        s <- mergedPars$s
        p <- mergedPars$p
        nG <- length(m)
    }
    peaks <- sapply(seq_len(nG), function(ii) {
        d <- dnorm(simulLR, m[ii], sqrt(s[ii]))
        max(d * p[ii])
    })
    bestPeak <- which(peaks >= max(peaks) * peakThresh)[1]
    if (verbose) {
        message("Gaussian mixture estimation:")
        message("n.peaks =  ", nG)
        message("\nGroup parameters:")
        for (grp in seq_len(nG)) {
            msg <- sprintf("Grp %s:\nprop: %s,\tmean: %s,\tSd: %s,\tpeak height: %s", 
                grp, round(p[grp], 3), round(m[grp], 3), round(sqrt(s[grp]), 
                  3), round(peaks[grp], 3))
            message(msg)
        }
        message()
    }
    correct <- m[bestPeak]
    if (verbose) 
        message("Correction value:  ", round(correct, 3))
    if (is.na(Title)) {
        Title <- sprintf("%s\nCorrection value = %s", getInfo(object, 
            "sampleName"), round(correct, 5))
    }
    segTable$seg.mean <- segTable$seg.mean - correct
    segTable$seg.med <- segTable$seg.med - correct
    segTable <- rCGH:::.estimateCopy(segTable, ploidy)
    cnSet <- getCNset(object)
    cnSet$Log2Ratio <- cnSet$Log2Ratio - correct
    cnSet$estimCopy <- rCGH:::.probeCopyValue(segTable)
    object@segTable <- segTable
    object@cnSet <- cnSet
    object@param$EMcentralized <- TRUE
    object@param$nPeak <- nG
    object@param$peakProp <- as.numeric(p)
    object@param$peakMeans <- as.numeric(m)
    object@param$peakSigmaSq <- as.numeric(s)
    object@param$centralPeak <- as.numeric(bestPeak)
    object@param$correctionValue <- as.numeric(correct)
    if (verbose) 
        message("Use plotDensity() to visualize the LRR densities.")
    options(op)
    return(object)
}





















############################################################### methods


segmentCGH_custom = function (object, Smooth = TRUE, UndoSD = NULL, minLen = 10, nCores = NULL, verbose = TRUE) 
{
    if (!rCGH:::.validrCGHObject(object)) 
        return(NULL)
    ploidy <- as.numeric(getInfo(object, "ploidy"))
    cnSet <- getCNset(object)
    cnSet <- cnSet[order(cnSet$ChrNum, cnSet$ChrStart), ]
    params <- getParam(object)
    params$minSegLen <- minLen
    if (Smooth) {
        mad <- rCGH:::.getMAD(object)
        params$ksmooth <- floor(150 * mad) * 2 + 1
    }
    if (is.null(UndoSD)) {
        mad <- rCGH:::.getMAD(object)
        alpha <- 0.5
        if (inherits(object, "rCGH-Illumina")) {
            alpha <- 0.95
        }
        if (inherits(object, "rCGH-oncoScan")) {
            alpha <- 0.3
        }
        params$UndoSD <- alpha * mad^(1/2)
    } else {
        params$UndoSD <- UndoSD
    }
    L2R <- cnSet$Log2Ratio
    Chr <- cnSet$ChrNum
    Pos <- cnSet$ChrStart
    sampleName <- getInfo(object, "sampleName")
    if (is.na(sampleName)) {
        sampleName <- "sample_x"
    }
    nCores <- rCGH:::.setCores(nCores)
    if (verbose) {
        usd <- format(params$UndoSD, digits = 3)
        message("Computing LRR segmentation using UndoSD: ", 
                usd)
    }
    segTable <- rCGH:::.computeSegmentation(L2R, Chr, Pos, sampleName, params, nCores)
    if (!is.null(minLen) && minLen < 0) {
        message("'minLen', the minimal segment length can't be < 0")
        minLen <- NULL
    }
    if (!is.null(minLen)) {
        if (verbose) 
            message("Merging segments shorter than ", minLen, "Kb.")
            print("marker 1")
        segTable <- rCGH:::.smoothSeg(segTable, minLen)
    }
    print("marker 2")
    segTable <- rCGH:::.computeMedSegm(segTable, L2R)
    print("marker 3")
    # segTable$chrom = as.numeric(segTable$chrom)
    # segTable$loc.start = as.numeric(segTable$loc.start)
    # segTable$loc.end = as.numeric(segTable$loc.end)
    # segTable$num.mark = as.numeric(segTable$num.mark)
    # segTable$seg.mean = as.numeric(segTable$seg.mean)
    # segTable$seg.med = as.numeric(segTable$seg.med)
    # segTable$probes.Sd = as.numeric(segTable$probes.Sd)
    for (col in 1:length(segTable)) {
        currCol = segTable[[col]]
        # print(c("currCol: ", currCol))
        print(table(names(currCol)))
    }
    segTable = segTable[1:length(segTable[[1]])-1,]
    segTable <- rCGH:::.mergeLevels(segTable)
    print("marker 4")
    segTable <- rCGH:::.estimateCopy(segTable, ploidy)
    print("marker 5")
    probeValues <- rCGH:::.probeSegValue(segTable)
    if (verbose) 
        message("Number of segments: ", nrow(segTable))
    params$nSegment <- nrow(segTable)
    object@param <- params
    object@segTable <- segTable
    object@cnSet <- cbind.data.frame(cnSet, Segm = probeValues)
    return(object)
}




## conclusion of this: I don't know how CN are estimated during segmentation phase. It is made using simple math, Running it with our samples should tell whether it is interesting or not.

.estimateCopy = function (st, ploidy, expect = log2(seq(1, 60)/2)) {
    ## st = seg table
    ## expect = expected distribution
    P <- lapply(1:nrow(st), function(ii) {
        mi <- st$seg.med[ii] #seg.med = segment median
        if (mi > 5) 
            return(c(rep(0, length(expect) - 1), 1))
        si <- st$probes.Sd[ii]
        p <- sapply(expect, function(e) {
            dnorm(mi, e, si)
        })
        return(p)
    })
    #P contains a distribution: either a dnorm or a vector : c(0,0,0,0,0...0,0,0,1) (30 values total). p too
    ratio <- sapply(P, function(p) { ## here we apply on P a function that returns either 0 or 2**(a value of expect) at each run
        rCGH:::.estimateRatio(p, expect)
    })
    copies <- ifelse(ratio == 0, 0, ifelse(ratio > 0, 2 * ratio, 
        -1/ratio))
    st$estimCopy <- copies + (ploidy - 2)
    return(st)
}


.estimateRatio = function (p, expect) 
{
    if (max(p, na.rm = TRUE) < 0.001) 
        return(0)
    return(2^expect[which.max(p)]) ## return the value of expect that is the same index as p's max value
}





plotLOH_custom = function (object, Title = NULL) 
{
    if (!rCGH:::.validrCGHObject(object)) 
        return(NULL)
    hg18 <- hg18
    hg19 <- hg19
    hg38 <- hg38
    HG <- switch(getInfo(object, "genome"), hg18 = hg18, 
        hg19 = hg19, hg38 = hg38)
    cnSet <- getCNset(object)
    if (!"modelAllDif" %in% colnames(cnSet)) {
        message("No data available for plotting LOH")
        return(NULL)
    }
    snpSet <- cnSet[grep("^S|^rs", cnSet$ProbeName), ]
    if (is.null(Title)) {
        Title = paste(getInfo(object, "sampleName"), "-", 
            getInfo(object, "analysisDate"))
    }
    ss <- split(snpSet, snpSet$ChrNum)
    gLocs <- lapply(ss, function(tmp) {
        n <- nrow(tmp)
        chr <- unique(tmp$ChrNum)
        return(tmp$ChrStart + HG$cumlen[chr])
    })
    gLocs <- do.call(c, gLocs)
    marks <- sapply(2:nrow(HG), function(ii) (HG$cumlen[ii - 
        1] + HG$cumlen[ii])/2)
    idx <- sample(1:nrow(snpSet), min(nrow(snpSet), 50000))
    X <- data.frame(loc = gLocs[idx], AD = snpSet$Allele.Difference[idx], 
        adjAD = snpSet$modelAllDif[idx])
    cumCentr <- 1/2 * HG$length + HG$cumlen
    gPlot <- ggplot(data = X, aes_string(x = "loc", y = "AD")) + 
        geom_point(pch = 19, cex = 0.1, color = rgb(0.4, 1, 1, 
            0.1)) + geom_point(aes_string(y = "adjAD"), 
        pch = 19, cex = 0.2, color = rgb(0, 0, 0, 0.75)) + geom_hline(yintercept = c(-1, 
        1), color = "blue", linetype = 2) + geom_hline(yintercept = seq(-1.5, 
        1.5, by = 1), color = "lightblue", linetype = 2) + 
        geom_vline(xintercept = HG$cumlen[1:23], color = "red", 
            linetype = 2, size = 0.25) + ggtitle(Title) + xlab("Genomic position (bp)") + 
        ylab("Allelic Difference") + theme_bw() + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.margin = unit(c(0, 
            4, 4, 0), "mm"), plot.title = element_text(lineheight = 0.8, 
            size = rel(2), face = "bold"), axis.title.x = element_text(size = rel(1.8), 
            angle = 0), axis.text.x = element_text(size = rel(1.5)), 
        axis.title.y = element_text(size = rel(1.8), angle = 90), 
        axis.text.y = element_text(size = rel(1.5)))
    if (inherits(object, "rCGH-Illumina")) {
        gPlot <- gPlot + coord_cartesian(ylim = range(0, 1.25)) + 
            scale_y_continuous(breaks = seq(0, 1, by = 0.25)) + 
            annotate("text", x = c(-1e+08, cumCentr[1:23]), 
                y = rep(1.2, 24), label = c("Chr", seq(1, 
                  23)), size = 4, colour = "grey30")
    }
    else {
        gPlot <- gPlot + coord_cartesian(ylim = range(-1.99, 
            1.99)) + scale_y_continuous(breaks = seq(-1.5, 1.5, 
            by = 0.5)) + annotate("text", x = c(-1e+08, 
            cumCentr[1:23]), y = rep(1.75, 24), label = c("Chr", 
            seq(1, 23)), size = 4, colour = "grey30")
    }
    return(gPlot)
}




adjustSignal_custom = function (object, Scale = TRUE, Cy = TRUE, GC = TRUE, Ref = "cy3", 
          suppOutliers = TRUE, nCores = NULL, verbose = TRUE) 
{
    if (!rCGH:::.validrCGHObject(object)) 
        return(NULL)
    cnSet <- getCNset(object)
    if (!inherits(object, "rCGH-Agilent")) {
        Cy <- FALSE
        GC <- FALSE
    }
    if (Cy) {
        if (verbose) {
            message("Recall you are using ", Ref, " as reference.")
            message("Cy effect adjustment...")
        }
        cnSet <- rCGH:::.CyAdjust(cnSet, Ref)
    }
    if (GC) {
        if (verbose) 
            message("GC% adjustment...")
        cnSet <- rCGH:::.GCadjust(cnSet)
    }
    object@param$CyAdjusted = Cy
    object@param$GCAdjusted = GC
    object@param$dLRs <- rCGH:::.dlrs(cnSet$Log2Ratio)
    object@param$MAD <- rCGH:::.MAD(cnSet$Log2Ratio)
    if (verbose) {
        message("Log2Ratios QCs:")
        message("\tdLRs: ", round(object@param$dLRs, 3))
        message("\tMAD: ", round(object@param$MAD, 3))
        message()
    }
    if (Scale) {
        if (verbose) 
            message("Scaling...")
        cnSet$Log2Ratio <- scale(cnSet$Log2Ratio, center = FALSE)
        cnSet$Log2Ratio <- cnSet$Log2Ratio * 1.2
    }
    nCores <- rCGH:::.setCores(nCores, verbose)
    if (suppOutliers) {
        if (verbose) 
            message("Signal filtering...")
        L2R <- cnSet$Log2Ratio
        Chr <- cnSet$ChrNum
        S <- NA
        cnSet$Log2Ratio <- rCGH:::.modelSignal(L2R, Chr, G = 1:5, method = "lr", 
                                        alpha = 2000, S, nCores, verbose)
    }
    print(c("Allele.Difference %in% colnames(cnSet) = ", "Allele.Difference" %in% colnames(cnSet)))
    if ("Allele.Difference" %in% colnames(cnSet)) {
        if (!all(is.na(cnSet$Allele.Difference))) {
            if (verbose) 
                message("Modeling allelic Difference...")
            signal <- cnSet$Allele.Difference
            chr <- cnSet$ChrNum
            S <- ifelse(inherits(object, "rCGH-oncoScan"), 
                        0.05, 0.04)
            print("- modelling allele diff -")
            modelAllDif <- rCGH:::.modelSignal(signal, chr, G = 2:5, 
                                        method = "loh", alpha = 2500, S, nCores, 
                                        verbose)
            print(c("modelAllDif: ", modelAllDif))                            
            cnSet$modelAllDif <- modelAllDif
        }
    }
    object@cnSet <- cnSet
    return(object)
}












#################################################### for error "Value of SET_STRING_ELT() must be a 'CHARSXP' not a 'NULL'" when running rCGH::segmentCGH().   ## 07/06/2022

## aCGH::mergeLevels
mergeLevels = function (vecObs, vecPred, pv.thres = 1e-04, ansari.sign = 0.05, 
    thresMin = 0.05, thresMax = 0.5, verbose = 1, scale = TRUE) 
{
    if (thresMin > thresMax) {
        cat("Error, thresMax should be equal to or larger than thresMin\n")
        return()
    }
    thresAbs = thresMin
    sq <- numeric()
    j = 0
    ansari = numeric()
    lv = numeric()
    flag = 0
    if (thresMin == thresMax) {
        flag = 2
    }
    else {
        l.step <- signif((thresMax - thresMin)/10, 1)
        s.step <- signif((thresMax - thresMin)/200, 1)
    }
    while (1) {
        if (verbose >= 1) {
            cat("\nCurrent thresAbs: ", thresAbs, "\n")
        }
        j = j + 1
        sq[j] <- thresAbs
        vecPredNow = vecPred
        mnNow = unique(vecPred)
        mnNow = mnNow[!is.na(mnNow)]
        cont = 0
        while (cont == 0 & length(mnNow) > 1) {
            mnNow = sort(mnNow)
            n <- length(mnNow)
            if (verbose >= 2) {
                cat("\r", n, ":", length(unique(vecPred)), "\t")
            }
            if (scale) {
                d <- (2 * 2^mnNow)[-n] - (2 * 2^mnNow)[-1]
            }
            else {
                d <- (mnNow)[-n] - (mnNow)[-1]
            }
            dst <- cbind(abs(d)[order(abs(d))], (2:n)[order(abs(d))], 
                (1:(n - 1))[order(abs(d))])
            for (i in 1:nrow(dst)) {
                cont = 1
                out = combine.func(diff = dst[i, 1], vecObs, 
                  vecPredNow, mnNow, mn1 = mnNow[dst[i, 2]], 
                  mn2 = mnNow[dst[i, 3]], pv.thres = pv.thres, 
                  thresAbs = if (scale) {
                    2 * 2^thresAbs - 2
                  }
                  else {
                    thresAbs
                  })
                if (out$pv > pv.thres) {
                  cont = 0
                  vecPredNow = out$vecPredNow
                  mnNow = out$mnNow
                  break
                }
            }
        }
        ansari[j] = ansari.test(sort(vecObs - vecPredNow), sort(vecObs - 
            vecPred))$p.value
        if (is.na(ansari[j])) {
            ansari[j] = 0
        }
        lv[j] = length(mnNow)
        if (flag == 2) {
            break
        }
        if (ansari[j] < ansari.sign) {
            flag = 1
        }
        if (flag) {
            if (ansari[j] > ansari.sign | thresAbs == thresMin) {
                break
            }
            else {
                thresAbs = signif(thresAbs - s.step, 3)
                if (thresAbs <= thresMin) {
                  thresAbs = thresMin
                }
            }
        }
        else {
            thresAbs = thresAbs + l.step
        }
        if (thresAbs >= thresMax) {
            thresAbs = thresMax
            flag = 2
        }
    }
    return(list(vecMerged = vecPredNow, mnNow = mnNow, sq = sq, 
        ansari = ansari))
}




byGeneTable_custom = function (segTable, symbol = NULL, genome = c("hg19", "hg18", "hg38"), columns = NA, verbose = TRUE) 
{
    genome <- match.arg(genome)
    .ByGene_custom(segTable, symbol, genome, columns, verbose)
}



.ByGene_custom = function (segTable, symbol, genome, columns, verbose) 
{
    hg18 <- hg18
    hg19 <- hg19
    hg38 <- hg38
    HG <- switch(genome, hg18 = hg18, hg19 = hg19, hg38 = hg38)
    geneDB <- rCGH:::.createGeneDB(genome)
    print(c("geneDB: ", geneDB))
    if (!"seg.med" %in% colnames(segTable)) 
        segTable$seg.med <- segTable$seg.mean
    if (!is.null(symbol)) 
        return(rCGH:::.getSegFromGene(segTable, symbol, HG, geneDB, 
            columns))
    if (verbose) 
        message("Creating byGene table...")
    bygene <- lapply(seq_len(nrow(segTable)), function(ii) {
        chr <- segTable$chrom[ii]
        Start <- segTable$loc.start[ii]
        End <- segTable$loc.end[ii]
        l <- abs(End - Start)/1000
        lrr <- segTable$seg.med[ii]
        nm <- segTable$num.mark[ii]
        copy <- segTable$estimCopy[ii]
        g <- rCGH:::.getGenesFromSeg(chr, Start, End, geneDB, columns)
        if (is.null(g)) 
            return(NULL)
        cbind.data.frame(g, Log2Ratio = lrr, num.mark = nm, segNum = ii, 
            `segLength(kb)` = round(l, 2), estimCopy = ifelse(is.null(copy), 
                NA, copy))
    })
    bygene <- do.call(rbind, bygene)
    cmValues <- rCGH:::.cmValues(segTable, HG)
    bygene$relativeLog <- rCGH:::.relativeLog(bygene, cmValues, HG)
    rCGH:::.addGenomeLoc(bygene, HG)
}
