## functions to plot values for segments of allele A and B
plotSeg = function(seg_df, allChrs, drawPolygons=F, drawDistances=F){
    lengthOfChrs = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432, 154913754, 57772954)
    print(lengthOfChrs)
    ## get nb cols
    nbcols = length(seg_df)
    ## get current chromosome
    currChr = which(allChrs==seg_df["chr"])
    # print(c("currChr: ", currChr))
    # print(c("class(currChr): ", class(currChr)))
    ## get endpos of last segment of previous chromosome
    if((seg_df["chr"]=="X" ) || (seg_df["chr"]=="Y")){
        pos0CurrChr = sum(lengthOfChrs[1:22])
    }
    else {
        print(c("currChr in ASCAT.R: ", currChr))
        if (currChr!=1){
            pos0CurrChr = sum(lengthOfChrs[1:currChr-1])
        } else {
            pos0CurrChr=0
        }
    }
    ## drawing a segment on plot for each segment of the genome
    segStartPos = as.numeric(seg_df[3])
    segStartPos = segStartPos + pos0CurrChr
    segEndPos = as.numeric(seg_df[4])
    # CN values are on last 2 columns
    A_allele_vals = nbcols-1
    deviation = 0.01
    logRAlleleA = as.numeric(seg_df[[nbcols-1]])
    logRAlleleB = as.numeric(seg_df[[nbcols]])
    logRAlleleA_dev = logRAlleleA + deviation
    logRAlleleB_dev = logRAlleleB - deviation
    # print(c("pos0CurrChr, segStartPos: ", pos0CurrChr+segStartPos))
    segments(pos0CurrChr+segStartPos, logRAlleleA_dev, pos0CurrChr+segEndPos, logRAlleleA_dev, col="dark blue", lwd=2)
    segments(pos0CurrChr+segStartPos, logRAlleleB_dev, pos0CurrChr+segEndPos, logRAlleleB_dev, col="dark red", lwd=2)
    ## drawing polygone joining this segment to closest non-null integer
    if(drawPolygons){
        polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos), 
                y=c(logRAlleleB_dev,logRAlleleB_dev,round(logRAlleleB),round(logRAlleleB)), col = "dark red", density = 10,angle=135)
        polygon(x=c(pos0CurrChr+segStartPos, pos0CurrChr+segEndPos, pos0CurrChr+segEndPos, pos0CurrChr+segStartPos),
                y=c(logRAlleleA_dev,logRAlleleA_dev,round(logRAlleleA),round(logRAlleleA)), col = "dark blue", density = 10, angle=45)
    }
    if(drawDistances) { ## distances between each segment and closest non-negative integer
        pass
    }
}

generateGrid = function(graph_title) {
    #create empty plot to add things in
    plot(1, ylim=c(0,2), xlim=c(0,3*10**9),col="white", xaxt="n", yaxt="n", xlab="nombre de bases", ylab="nombre de copies", main=graph_title)
    #  add horizontal grid
    for (CN in c(-3:8)) {
        abline(h=CN, col="dark grey")
    }
    # X-axis
    axis(1, at = c(0, 5*10**8, 1*10**9, 1.5*10**9, 2*10**9, 2.5*10**9, 3*10**9))
    # Y-axis
    axis(2, at = c(-3:8))
}

## functions for GI calculation
getNbChrs = function(segmentsTable) {
    # print(c("segmentsTable: ", segmentsTable))
    chrs = as.vector(segmentsTable$chr)
    # print(c("chrs: ", chrs))
    nbChr = length(unique(chrs))
    # print(c("nbChr: ", nbChr))
    return(nbChr)
}

calcGI_ASCAT = function(segmentsTable) {
    # print(c("segmentsTable: ", segmentsTable))
    nbChr = getNbChrs(segmentsTable)
    # print(c("nbChr: ", nbChr))
    nbAlter = dim(segmentsTable)[1]
    # print(c("nbAlter: ", nbAlter))
    GI = calcGI(nbAlter, nbChr)
    return(list(GI,nbAlter,nbChr))
}

removeSexChrData = function(rawData) {
    ## remove sex chromosomes data
    ### from ch
    rawData_ch = rawData[["data"]][["ch"]]
    rawData[["data"]][["ch"]] = rawData_ch[names(rawData_ch) %in% c("chrX", "chrY") == FALSE] 
    ### from chr
    lastChrProbes = rawData[["data"]][["ch"]][[length(rawData[["data"]][["ch"]])]]
    lastProbe = lastChrProbes[length(lastChrProbes)]
    chr = rawData[["data"]][["chr"]]
    for(i in length(chr):1) {
        # print(c("i: ", i))
        DNAregion = chr[i]
        if(any(as.vector(DNAregion[[1]])<=lastProbe)) {
            regionLen = length(DNAregion[[1]])
            # print(c("regionLen: ", regionLen))
            # print(c("last probe of this region: ", (DNAregion[[1]][regionLen])))
            if(DNAregion[[1]][regionLen]==lastProbe){
                # print("Last probe of this region is the last probe of chromosome 22")
                chr = chr[1:i]
            } else {
                # print("last probe of chr 22 is not last probe of this region")
                chr = chr[1:i-1]
            }
            break
        }
            
    }
    rawData[["data"]][["chr"]] = chr
    ### from chrs
    rawData_chrs = rawData[["data"]][["chrs"]]
    rawData_chrs = rawData_chrs[rawData_chrs!="chrX"]; rawData_chrs = rawData_chrs[rawData_chrs!="chrY"]
    rawData[["data"]][["chrs"]] = rawData_chrs
    ### from Tumor_*
    rawData[["data"]][["Tumor_LogR.ori"]] = rawData[["data"]][["Tumor_LogR.ori"]] %>% dplyr::slice(1:lastProbe)
    rawData[["data"]][["Tumor_LogR"]] = rawData[["data"]][["Tumor_LogR"]] %>% dplyr::slice(1:lastProbe)
    rawData[["data"]][["Tumor_BAF"]] = rawData[["data"]][["Tumor_BAF"]] %>% dplyr::slice(1:lastProbe)
    rawData[["data"]][["Tumor_AD"]] = rawData[["data"]][["Tumor_AD"]] %>% dplyr::slice(1:lastProbe)
    ### from SNPpos
    rawData[["data"]][["SNPpos"]] = rawData[["data"]][["SNPpos"]] %>% dplyr::slice(1:lastProbe)
    ### from additional data
    rawData[["data"]][["additional"]] = rawData[["data"]][["additional"]] %>% dplyr::slice(1:lastProbe)

    return(rawData)
}

cleanASCATSegData = function(callData, trimData = F) {
    ## get seg data from call result
    # segTable_raw = callData$segments_raw
    segTable = callData$segments
    
    ## check seg data
    # graph_title = paste0(sampleName, " called data")
    # generateGrid(graph_title)
    # allChrs = unique(as.vector(segTable[2]))
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    
    if(trimData) {
        ## removing segments shorter than 300 Kbp
        segTable = dplyr::filter(segTable, endpos-startpos>300000)
    }
    # graph_title = paste0(sampleName, " after removing segments shorter than 300 Kbp")
    # generateGrid(graph_title)
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    
    ## removing segments with CN=1 for each allele
    segTableClean = dplyr::filter(segTable, nMajor!=1 | nMinor!=1)
    nbOrigin = dim(segTable)[1]
    nbClean = dim(segTableClean)[1]
    nbRemoved = nbOrigin-nbClean
    print(paste0(nbRemoved, " segments removed out of ", nbOrigin))
    # graph_title = paste0(sampleName, " after removing segments of copy number=2")
    # generateGrid(graph_title)
    # apply(segTable, 1, plotSeg, allChrs, lengthOfChrs)
    return(segTableClean)   
}
    









































    


custom_ascat.aspcf = function (ASCATobj, selectsamples = 1:length(ASCATobj$samples), 
    ascat.gg = NULL, penalty = 70, out.dir = ".", out.prefix = "") 
{
    gg = NULL
    if (!is.null(ascat.gg)) {
        gg = ascat.gg$germlinegenotypes
    }
    else {
        gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 
            0.7
    }
    ghs = ASCAT:::predictGermlineHomozygousStretches(ASCATobj$chr, gg)
    segmentlengths = unique(c(penalty, 35, 50, 70, 100, 140))
    segmentlengths = segmentlengths[segmentlengths >= penalty]
    Tumor_LogR_segmented = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], 
        ncol = dim(ASCATobj$Tumor_LogR)[2])
    rownames(Tumor_LogR_segmented) = rownames(ASCATobj$Tumor_LogR)
    colnames(Tumor_LogR_segmented) = colnames(ASCATobj$Tumor_LogR)
    Tumor_BAF_segmented = list()
    print("mark")
    for (sample in selectsamples) {
        print.noquote(paste("Sample ", ASCATobj$samples[sample], 
            " (", sample, "/", length(ASCATobj$samples), 
            ")", sep = ""))
        logrfilename = file.path(out.dir, paste(out.prefix, filename = distancepng, width = 1000, 
            ".LogR.PCFed.txt", sep = ""))
        baffilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], 
            ".BAF.PCFed.txt", sep = ""))
        logRPCFed = numeric(0)
        bafPCFed = numeric(0)
        if (!is.null(ASCATobj$X_nonPAR) && ASCATobj$gender[sample] == 
            "XY") {
            nonPAR_index = which(ASCATobj$SNPpos$Chr == "X" & 
                ASCATobj$SNPpos$Position >= ASCATobj$X_nonPAR[1] & 
                ASCATobj$SNPpos$Position <= ASCATobj$X_nonPAR[2] & 
                !is.na(gg[, sample]))
            autosomes_info = table(gg[which(ASCATobj$SNPpos$Chr %in% 
                setdiff(ASCATobj$chrs, ASCATobj$sexchromosomes)), 
                sample])
            if (length(nonPAR_index) > 5) {
                gg[nonPAR_index, sample] = T
                if (!is.null(ASCATobj$Germline_BAF)) {
                    DIST = 1 - sapply(ASCATobj$Germline_BAF[nonPAR_index, sample], function(x) {
                        if (x > 0.5) 
                            return(x)
                        else return(1 - x)
                        })
                    gg[nonPAR_index[which(rank(DIST, ties.method = "random") <= round(length(DIST) * (autosomes_info["FALSE"]/sum(autosomes_info))))], sample] = F
                    rm(DIST)
                }
                else {
                    gg[sample(nonPAR_index, round(length(nonPAR_index) * (autosomes_info["FALSE"]/sum(autosomes_info)))), sample] = F
                }
            }
            rm(nonPAR_index, autosomes_info)
        }
        for (segmentlength in segmentlengths) {
            logRPCFed = numeric(0)
            bafPCFed = numeric(0)
            tbsam = ASCATobj$Tumor_BAF[, sample]
            names(tbsam) = rownames(ASCATobj$Tumor_BAF)
            homosam = gg[, sample]
            for (chrke in 1:length(ASCATobj$chr)) {
                lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]], sample]
                lrwins = vector(mode = "numeric", length = length(lr))
                lrwins[is.na(lr)] = NA
                lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)], 2.5, 25)$ywin
                baf = tbsam[ASCATobj$chr[[chrke]]]
                homo = homosam[ASCATobj$chr[[chrke]]]
                Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
                bafsel = baf[Select_het]
                bafselwinsmirrored = madWins(ifelse(bafsel > 0.5, bafsel, 1 - bafsel), 2.5, 25)$ywin
                bafselwins = ifelse(bafsel > 0.5, bafselwinsmirrored, 1 - bafselwinsmirrored)
                indices = which(Select_het)
                logRaveraged = NULL
                if (length(indices) != 0) {
                    averageIndices = c(1, (indices[1:(length(indices) - 1)] + indices[2:length(indices)])/2, length(lr) + 0.01)
                    startindices = ceiling(averageIndices[1:(length(averageIndices) - 1)])
                    endindices = floor(averageIndices[2:length(averageIndices)] - 0.01)
                    if (length(indices) == 1) {
                    startindices = 1
                    endindices = length(lr)
                    }
                    nrIndices = endindices - startindices + 1
                    logRaveraged = vector(mode = "numeric", length = length(indices))
                    for (i in 1:length(indices)) {
                    if (is.na(endindices[i])) {
                        endindices[i] = startindices[i]
                    }
                    logRaveraged[i] = mean(lrwins[startindices[i]:endindices[i]], na.rm = T)
                    }
                }
                if (length(logRaveraged) > 0) {
                    logRASPCF = NULL
                    bafASPCF = NULL
                    if (length(logRaveraged) < 6) {
                    logRASPCF = rep(mean(logRaveraged), length(logRaveraged))
                    bafASPCF = rep(mean(bafselwins), length(logRaveraged))
                    }
                    else {
                    PCFed = fastAspcf(logRaveraged, bafselwins, 6, segmentlength) #(logR, allB, kmin, gamma)
                    logRASPCF = PCFed$yhat1
                    bafASPCF = PCFed$yhat2
                    }
                    names(bafASPCF) = names(indices)
                    logRc = numeric(0)
                    for (probe in 1:length(logRASPCF)) {
                        if (probe == 1) {
                            logRc = rep(logRASPCF[probe], indices[probe])
                        }
                        if (probe == length(logRASPCF)) {
                            logRc = c(logRc, rep(logRASPCF[probe], 
                            length(lr) - indices[probe]))
                        }
                        else if (logRASPCF[probe] == logRASPCF[probe + 
                            1]) {
                            logRc = c(logRc, rep(logRASPCF[probe], 
                            indices[probe + 1] - indices[probe]))
                        }
                        else {
                            d = numeric(0)
                            totall = indices[probe + 1] - indices[probe]
                            for (bp in 0:(totall - 1)) {
                            dis = sum(abs(lr[(1:bp) + indices[probe]] - 
                                logRASPCF[probe]), na.rm = T)
                            if (bp != totall) {
                                dis = sum(dis, sum(abs(lr[((bp + 1):totall) + 
                                indices[probe]] - logRASPCF[probe + 
                                1]), na.rm = T), na.rm = T)
                            }
                            d = c(d, dis)
                            }
                            breakpoint = which.min(d) - 1
                            logRc = c(logRc, rep(logRASPCF[probe], breakpoint), rep(logRASPCF[probe + 1], totall - breakpoint))
                        }
                    }
                    logRd = numeric(0)
                    seg = rle(logRc)$lengths
                    startprobe = 1
                    endprobe = 0
                    for (i in 1:length(seg)) {
                        endprobe = endprobe + seg[i]
                        level = mean(lr[startprobe:endprobe], na.rm = T)
                        logRd = c(logRd, rep(level, seg[i]))
                        startprobe = startprobe + seg[i]
                    }
                    logRPCFed = c(logRPCFed, logRd)
                    bafPCFed = c(bafPCFed, bafASPCF)
                }
                else {
                    level = mean(lr, na.rm = T)
                    reps = length(lr)
                    logRPCFed = c(logRPCFed, rep(level, reps))
                }
                homsegs = ghs[[sample]][ghs[[sample]][, 1] == chrke, ]
                startchr = min(ASCATobj$chr[[chrke]])
                endchr = max(ASCATobj$chr[[chrke]])
                if (length(homsegs) == 3) {
                    homsegs = t(as.matrix(homsegs))
                }
                if (!is.null(homsegs) && !is.na(homsegs) && dim(homsegs)[1] != 0) {
                    for (i in 1:dim(homsegs)[1]) {
                        startpos = max(homsegs[i, 2], startchr)
                        endpos = min(homsegs[i, 3], endchr)
                        startpos2 = max(homsegs[i, 2] - 100, startchr)
                        endpos2 = min(homsegs[i, 3] + 100, endchr)
                        startpos3 = max(homsegs[i, 2] - 5, startchr)
                        endpos3 = min(homsegs[i, 3] + 5, endchr)
                        towins = ASCATobj$Tumor_LogR[startpos2:endpos2, 
                        sample]
                        winsed = madWins(towins[!is.na(towins)], 
                        2.5, 25)$ywin
                        pcfed = vector(mode = "numeric", length = length(towins))
                        pcfed[!is.na(towins)] = exactPcf(winsed, 6, floor(segmentlength/4))
                        pcfed2 = pcfed[(startpos3 - startpos2 + 1):(endpos3 - 
                        startpos2 + 1)]
                        dif = abs(pcfed2 - logRPCFed[startpos3:endpos3])
                        if (!is.na(dif) && sum(dif > 0.3) > 5) {
                        logRPCFed[startpos3:endpos3] = ifelse(dif > 
                            0.3, pcfed2, logRPCFed[startpos3:endpos3])
                        }
                    }
                }
            }
            logRPCFed = fillNA(logRPCFed, zeroIsNA = TRUE)
            seg = rle(logRPCFed)$lengths
            logRPCFed = numeric(0)
            startprobe = 1
            endprobe = 0
            prevlevel = 0
            for (i in 1:length(seg)) {
                endprobe = endprobe + seg[i]
                level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe, sample], na.rm = T)
                if (is.nan(level)) {
                    level = prevlevel
                }
                else {
                    prevlevel = level
                }
                logRPCFed = c(logRPCFed, rep(level, seg[i]))
                startprobe = startprobe + seg[i]
            }
            names(logRPCFed) = rownames(ASCATobj$Tumor_LogR)
            if (length(unique(logRPCFed)) < 800) {
                break
            }
        }
        print(c("logRPCFed:" , logRPCFed))
        print(c("logrfilename:" , logrfilename))
        write.table(logRPCFed, logrfilename, sep = "\t", 
            col.names = F)
        write.table(bafPCFed, baffilename, sep = "\t", 
            col.names = F)
        bafPCFed = as.matrix(bafPCFed)
        Tumor_LogR_segmented[, sample] = logRPCFed
        Tumor_BAF_segmented[[sample]] = 1 - bafPCFed
    }
    ASCATobj$Tumor_LogR_segmented = Tumor_LogR_segmented
    ASCATobj$Tumor_BAF_segmented = Tumor_BAF_segmented
    ASCATobj$failedarrays = ascat.gg$failedarrays
    return(ASCATobj)
}



## je veux voir logrfilename donc je lance ASCAT. Le but est de voir si plusieurs logRPCFed sont créés en fichier texte et je peux peut-être les comparer.





fastAspcf = function (logR, allB, kmin, gamma) {
    N <- length(logR)
    w <- 1000
    d <- 100
    startw = -d
    stopw = w - d
    nseg = 0
    var2 = 0
    var3 = 0
    breakpts = 0
    larger = TRUE
    repeat {
        from <- max(c(1, startw))
        to <- min(c(stopw, N))
        logRpart <- logR[from:to]
        allBpart <- allB[from:to]
        allBflip <- allBpart
        allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]
        sd1 <- getMad(logRpart)
        sd2 <- getMad(allBflip)
        sd3 <- getMad(allBpart)
        sd.valid <- c(!is.na(sd1), !is.na(sd2), sd1 != 0, sd2 != 
            0)
        if (all(sd.valid)) {
            part.res <- aspcfpart(logRpart = logRpart, allBflip = allBflip, 
                a = startw, b = stopw, d = d, sd1 = sd1, sd2 = sd2, 
                N = N, kmin = kmin, gamma = gamma)
            breakptspart <- part.res$breakpts
            larger = breakptspart > breakpts[length(breakpts)]
            breakpts <- c(breakpts, breakptspart[larger])
            var2 <- var2 + sd2^2
            var3 <- var3 + sd3^2
            nseg = nseg + 1
        }
        if (stopw < N + d) {
            startw <- min(stopw - 2 * d + 1, N - 2 * d)
            stopw <- startw + w
        }
        else {
            break
        }
    }
    breakpts <- unique(c(breakpts, N))
    if (nseg == 0) {
        nseg = 1
    }
    sd2 <- sqrt(var2/nseg)
    sd3 <- sqrt(var3/nseg)
    frst <- breakpts[1:length(breakpts) - 1] + 1
    last <- breakpts[2:length(breakpts)]
    nseg <- length(frst)
    yhat1 <- rep(NA, N)
    yhat2 <- rep(NA, N)
    for (i in 1:nseg) {
        yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), 
            last[i] - frst[i] + 1)
        yi2 <- allB[frst[i]:last[i]]
        if (length(yi2) == 0) {
            mu <- 0
        }
        else {
            mu <- mean(abs(yi2 - 0.5))
        }
        if (sqrt(sd2^2 + mu^2) < 2 * sd2) {
            mu <- 0
        }
        yhat2[frst[i]:last[i]] <- rep(mu + 0.5, last[i] - frst[i] + 
            1)
    }
    return(list(yhat1 = yhat1, yhat2 = yhat2))
}
















custom_ascat.runAscat = function (ASCATobj, gamma = 0.55, pdfPlot = F, y_limit = 5, circos = NA, 
    min_ploidy = 1.5, max_ploidy = 5.5, rho_manual = NA, psi_manual = NA, 
    img.dir = ".", img.prefix = "") 
{
    goodarrays = NULL
    N_samples = dim(ASCATobj$Tumor_LogR)[2]
    print(c("N_samples  : ", N_samples))
    res = vector("list", N_samples)
    stopifnot(length(rho_manual) == length(psi_manual))
    if (length(rho_manual) == 1 && is.na(rho_manual) && N_samples > 
        1) {
        rho_manual = rep(NA, N_samples)
        psi_manual = rep(NA, N_samples)
    }
    else {
        stopifnot(length(rho_manual) == N_samples)
    }
    for (arraynr in 1:N_samples) {
        print.noquote(paste("Sample ", ASCATobj$samples[arraynr], 
            " (", arraynr, "/", length(ASCATobj$samples), 
            ")", sep = ""))
        lrr = ASCATobj$Tumor_LogR[, arraynr]
        names(lrr) = rownames(ASCATobj$Tumor_LogR)
        baf = ASCATobj$Tumor_BAF[, arraynr]
        names(baf) = rownames(ASCATobj$Tumor_BAF)
        lrrsegm = ASCATobj$Tumor_LogR_segmented[, arraynr]
        names(lrrsegm) = rownames(ASCATobj$Tumor_LogR_segmented)
        bafsegm = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , 
            drop = FALSE]
        names(bafsegm) = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
        failedqualitycheck = F
        if (ASCATobj$samples[arraynr] %in% ASCATobj$failedarrays) {
            failedqualitycheck = T
        }
        ending = ifelse(pdfPlot, "pdf", "png")
        circosName = NA
        if (!is.na(circos)) {
            circosName = paste(circos, "_", ASCATobj$samples[arraynr], 
                sep = "")
        }
        ### printing paths given to save plots
        print("sunrise, ascatprofile, rawprofile:")
        print(file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".sunrise.png", sep = "")))
        print(file.path(img.dir, 
                paste(img.prefix, ASCATobj$samples[arraynr], ".ASCATprofile.", ending, sep = "")))
        print( file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], ".rawprofile.", ending, sep = "")))
        print(c("lrr: ", lrr))
        res[[arraynr]] = runASCAT(lrr, baf, lrrsegm, bafsegm, 
            ASCATobj$gender[arraynr], ASCATobj$SNPpos, ASCATobj$ch, 
            ASCATobj$chrs, ASCATobj$sexchromosomes, failedqualitycheck, 
            file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], 
                ".sunrise.png", sep = "")), file.path(img.dir, 
                paste(img.prefix, ASCATobj$samples[arraynr], 
                  ".ASCATprofile.", ending, sep = "")), 
            file.path(img.dir, paste(img.prefix, ASCATobj$samples[arraynr], 
                ".rawprofile.", ending, sep = "")), 
            NA, gamma, rho_manual[arraynr], psi_manual[arraynr], 
            pdfPlot, y_limit, circosName, min_ploidy, max_ploidy, 
            ASCATobj$X_nonPAR)
        if (!is.na(res[[arraynr]]$rho)) {
            goodarrays[length(goodarrays) + 1] = arraynr
        }
    }
    if (length(goodarrays) > 0) {
        n1 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
        n2 = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = length(goodarrays))
        rownames(n1) = rownames(ASCATobj$Tumor_LogR)
        rownames(n2) = rownames(ASCATobj$Tumor_LogR)
        colnames(n1) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
        colnames(n2) = colnames(ASCATobj$Tumor_LogR)[goodarrays]
        for (i in 1:length(goodarrays)) {
            n1[, i] = res[[goodarrays[i]]]$nA
            n2[, i] = res[[goodarrays[i]]]$nB
        }
        distance_matrix = vector("list", length(goodarrays))
        names(distance_matrix) <- colnames(ASCATobj$Tumor_LogR)[goodarrays]
        for (i in 1:length(goodarrays)) {
            distance_matrix[[i]] = res[[goodarrays[i]]]$distance_matrix
        }
        tp = vector(length = length(goodarrays))
        psi = vector(length = length(goodarrays))
        ploidy = vector(length = length(goodarrays))
        goodnessOfFit = vector(length = length(goodarrays))
        naarrays = NULL
        for (i in 1:length(goodarrays)) {
            tp[i] = res[[goodarrays[i]]]$rho
            psi[i] = res[[goodarrays[i]]]$psi
            ploidy[i] = mean(res[[goodarrays[i]]]$nA + res[[goodarrays[i]]]$nB, 
                na.rm = T)
            goodnessOfFit[i] = res[[goodarrays[i]]]$goodnessOfFit
            if (res[[goodarrays[i]]]$nonaberrant) {
                naarrays = c(naarrays, ASCATobj$samples[goodarrays[i]])
            }
        }
        fa = colnames(ASCATobj$Tumor_LogR)[-goodarrays]
        names(tp) = colnames(n1)
        names(ploidy) = colnames(n1)
        names(psi) = colnames(n1)
        names(goodnessOfFit) = colnames(n1)
        seg = NULL
        for (i in 1:length(goodarrays)) {
            segje = res[[goodarrays[i]]]$seg
            seg = rbind(seg, cbind(ASCATobj$samples[goodarrays[i]], 
                as.vector(ASCATobj$SNPpos[segje[, 1], 1]), ASCATobj$SNPpos[segje[, 
                  1], 2], ASCATobj$SNPpos[segje[, 2], 2], segje[, 
                  3], segje[, 4]))
        }
        colnames(seg) = c("sample", "chr", "startpos", 
            "endpos", "nMajor", "nMinor")
        seg = data.frame(seg, stringsAsFactors = F)
        seg[, 3] = as.numeric(seg[, 3])
        seg[, 4] = as.numeric(seg[, 4])
        seg[, 5] = as.numeric(seg[, 5])
        seg[, 6] = as.numeric(seg[, 6])
        seg_raw = NULL
        for (i in 1:length(goodarrays)) {
            segje = res[[goodarrays[i]]]$seg_raw
            seg_raw = rbind(seg_raw, cbind(ASCATobj$samples[goodarrays[i]], 
                as.vector(ASCATobj$SNPpos[segje[, 1], 1]), ASCATobj$SNPpos[segje[, 
                  1], 2], ASCATobj$SNPpos[segje[, 2], 2], segje[, 
                  3], segje[, 4:ncol(segje)]))
        }
        colnames(seg_raw) = c("sample", "chr", "startpos", 
            "endpos", "nMajor", "nMinor", "nAraw", 
            "nBraw")
        seg_raw = data.frame(seg_raw, stringsAsFactors = F)
        seg_raw[, 3] = as.numeric(seg_raw[, 3])
        seg_raw[, 4] = as.numeric(seg_raw[, 4])
        seg_raw[, 5] = as.numeric(seg_raw[, 5])
        seg_raw[, 6] = as.numeric(seg_raw[, 6])
        seg_raw[, 7] = as.numeric(seg_raw[, 7])
        seg_raw[, 8] = as.numeric(seg_raw[, 8])
    }
    else {
        n1 = NULL
        n2 = NULL
        tp = NULL
        ploidy = NULL
        psi = NULL
        goodnessOfFit = NULL
        fa = colnames(ASCATobj$Tumor_LogR)
        naarrays = NULL
        seg = NULL
        seg_raw = NULL
        distance_matrix = NULL
    }
    return(list(nA = n1, nB = n2, purity = tp, aberrantcellfraction = tp, 
        ploidy = ploidy, psi = psi, goodnessOfFit = goodnessOfFit, 
        failedarrays = fa, nonaberrantarrays = naarrays, segments = seg, 
        segments_raw = seg_raw, distance_matrix = distance_matrix))
}













































custom_ascat.metrics = function (ASCAT_input_object, ASCAT_output_object) 
{
    METRICS = do.call(rbind, lapply(1:length(ASCAT_input_object$samples), 
        function(nSAMPLE) {
            SAMPLE = ASCAT_input_object$samples[nSAMPLE]
            sex = ASCAT_input_object$gender[nSAMPLE]
            tumour_mapd = round(median(abs(diff(na.omit(ASCAT_input_object$Tumor_LogR[,SAMPLE])))), 4)
            if (!is.null(ASCAT_input_object$Germline_LogR) && 
                any(SAMPLE %in% colnames(ASCAT_input_object$Germline_LogR))) {
                normal_mapd = round(median(abs(diff(na.omit(ASCAT_input_object$Germline_LogR[, 
                  SAMPLE])))), 4)
            }
            else {
                normal_mapd = NA
            }
            if ("GC_correction_before" %in% names(ASCAT_input_object)) {
                GC_correction_before = ASCAT_input_object$GC_correction_before[SAMPLE]
            }
            else {
                GC_correction_before = NA
            }
            if ("GC_correction_after" %in% names(ASCAT_input_object)) {
                GC_correction_after = ASCAT_input_object$GC_correction_after[SAMPLE]
            }
            else {
                GC_correction_after = NA
            }
            if ("RT_correction_before" %in% names(ASCAT_input_object)) {
                RT_correction_before = ASCAT_input_object$RT_correction_before[SAMPLE]
            }
            else {
                RT_correction_before = NA
            }
            if ("RT_correction_after" %in% names(ASCAT_input_object)) {
                RT_correction_after = ASCAT_input_object$RT_correction_after[SAMPLE]
            }
            else {
                RT_correction_after = NA
            }
            if (!is.null(ASCAT_input_object$Tumor_LogR_segmented) && 
                !is.null(ASCAT_input_object$Tumor_BAF_segmented[[nSAMPLE]])) {
                n_het_SNP = length(ASCAT_input_object$Tumor_BAF_segmented[[nSAMPLE]])
                n_segs_logR = length(rle(ASCAT_input_object$Tumor_LogR_segmented[, 
                  SAMPLE])$values)
                n_segs_BAF = length(rle(ASCAT_input_object$Tumor_BAF_segmented[[nSAMPLE]][, 
                  1])$values)
                n_segs_logRBAF_diff = abs(n_segs_logR - n_segs_BAF)
                segm_baf = ASCAT_input_object$Tumor_BAF[rownames(ASCAT_input_object$Tumor_BAF_segmented[[nSAMPLE]]), 
                  SAMPLE]
                frac_homo = round(length(which(segm_baf < 0.1 | 
                  segm_baf > 0.9))/length(segm_baf), 4)
                rm(segm_baf)
            }
            else {
                n_het_SNP = NA
                n_segs_logR = NA
                n_segs_BAF = NA
                n_segs_logRBAF_diff = NA
                frac_homo = NA
            }
            if (!is.null(ASCAT_output_object$segments) && SAMPLE %in% 
                ASCAT_output_object$segments$sample) {
                purity = round(as.numeric(ASCAT_output_object$purity[SAMPLE]), 
                  4)
                ploidy = round(as.numeric(ASCAT_output_object$ploidy[SAMPLE]), 
                  4)
                goodness_of_fit = round(ASCAT_output_object$goodnessOfFit[SAMPLE], 
                  4)
                profile = ASCAT_output_object$segments[ASCAT_output_object$segments$sample == 
                  SAMPLE, ]
                profile$size = profile$endpos - profile$startpos + 
                  1
                n_segs = nrow(profile)
                segs_size = sum(profile$size)
                n_segs_1kSNP = round(n_segs/(length(ASCAT_input_object$Tumor_BAF_segmented[[nSAMPLE]])/1000), 
                  4)
                INDEX_HD = which(profile$nMajor == 0 & profile$nMinor == 
                  0)
                if (length(INDEX_HD) > 0) {
                  homdel_segs = length(INDEX_HD)
                  homdel_largest = max(profile$size[INDEX_HD])
                  homdel_size = sum(profile$size[INDEX_HD])
                  homdel_fraction = round(homdel_size/sum(profile$size), 
                    4)
                }
                else {
                  homdel_segs = homdel_largest = homdel_size = homdel_fraction = 0
                }
                rm(INDEX_HD)
                profile = profile[which(profile$chr %in% setdiff(ASCAT_input_object$chrs, 
                  ASCAT_input_object$sexchromosomes)), ]
                LOH = round(sum(profile$size[which(profile$nMinor == 
                  0)])/sum(profile$size), 4)
                mode_minA = modeAllele(profile, "nMinor")
                mode_majA = modeAllele(profile, "nMajor")
                if (mode_majA == 0 || !(mode_majA %in% 1:5)) {
                  WGD = NA
                  GI = NA
                }
                else {
                  if (mode_majA == 1) {
                    WGD = 0
                    GI = computeGIscore(WGD, profile)
                  }
                  else if (mode_majA == 2) {
                    WGD = 1
                    GI = computeGIscore(WGD, profile)
                  }
                  else if (mode_majA %in% 3:5) {
                    WGD = "1+"
                    GI = computeGIscore(1, profile)
                  }
                }
                rm(profile)
            }
            else {
                purity = NA
                ploidy = NA
                goodness_of_fit = NA
                n_segs = NA
                segs_size = NA
                n_segs_1kSNP = NA
                homdel_segs = NA
                homdel_largest = NA
                homdel_size = NA
                homdel_fraction = NA
                LOH = NA
                mode_minA = NA
                mode_majA = NA
                WGD = NA
                GI = NA
            }
            OUT = data.frame(sex = sex, tumour_mapd = tumour_mapd, 
                normal_mapd = normal_mapd, GC_correction_before = GC_correction_before, 
                GC_correction_after = GC_correction_after, RT_correction_before = RT_correction_before, 
                RT_correction_after = RT_correction_after, n_het_SNP = n_het_SNP, 
                n_segs_logR = n_segs_logR, n_segs_BAF = n_segs_BAF, 
                n_segs_logRBAF_diff = n_segs_logRBAF_diff, frac_homo = frac_homo, 
                purity = purity, ploidy = ploidy, goodness_of_fit = goodness_of_fit, 
                n_segs = n_segs, segs_size = segs_size, n_segs_1kSNP = n_segs_1kSNP, 
                homdel_segs = homdel_segs, homdel_largest = homdel_largest, 
                homdel_size = homdel_size, homdel_fraction = homdel_fraction, 
                LOH = LOH, mode_minA = mode_minA, mode_majA = mode_majA, 
                WGD = WGD, GI = GI, stringsAsFactors = F)
            rownames(OUT) = SAMPLE
            return(OUT)
        }))
    return(METRICS)
}


































custom_ascat.plotAdjustedAscatProfile = function (ASCAT_output_object, REF, y_limit = 5, plot_unrounded = F, 
    png_prefix = "") 
{
    if (plot_unrounded) {
        SEGMENTS = ASCAT_output_object$segments_raw[, c(1:4, 
            7:8)]
        colnames(SEGMENTS)[5:6] = c("nMajor", "nMinor")
        SEGMENTS$nMajor = SEGMENTS$nMajor + SEGMENTS$nMinor
        colourA = "#943CC3"
        colourB = "#60AF36"
    }
    else {
        SEGMENTS = ASCAT_output_object$segments
        SEGMENTS$nMajor = SEGMENTS$nMajor - 0.1
        SEGMENTS$nMinor = SEGMENTS$nMinor + 0.1
        colourA = "#E03546"
        colourB = "#3557E0"
    }
    SEGMENTS$nMajor = ifelse(SEGMENTS$nMajor > y_limit, y_limit + 
        0.1, SEGMENTS$nMajor)
    SEGMENTS$nMinor = ifelse(SEGMENTS$nMinor > y_limit, y_limit + 
        0.1, SEGMENTS$nMinor)
    if (REF == "hg19") {
        REF = data.frame(chrom = c(1:22, "X"), start = rep(1, 
            23), end = c(249250621, 243199373, 198022430, 191154276, 
            180915260, 171115067, 159138663, 146364022, 141213431, 
            135534747, 135006516, 133851895, 115169878, 107349540, 
            102531392, 90354753, 81195210, 78077248, 59128983, 
            63025520, 48129895, 51304566, 155270560))
    }
    else if (REF == "hg38") {
        REF = data.frame(chrom = c(1:22, "X"), start = rep(1, 
            23), end = c(248956422, 242193529, 198295559, 190214555, 
            181538259, 170805979, 159345973, 145138636, 138394717, 
            133797422, 135086622, 133275309, 114364328, 107043718, 
            101991189, 90338345, 83257441, 80373285, 58617616, 
            64444167, 46709983, 50818468, 156040895))
    }
    else {
        stopifnot(is.data.frame(REF))
        stopifnot(identical(colnames(REF), c("chrom", "start", 
            "end")))
    }
    SEGMENTS$chr = gsub("^chr", "", SEGMENTS$chr)
    print(c("REF$chrom: ", REF$chrom))
    stopifnot(all(ASCAT_output_object$segments$chr %in% REF$chrom))
    REF$size = REF$end - REF$start + 1
    REF$middle = 0
    for (i in 1:nrow(REF)) {
        if (i == 1) {
            REF$middle[i] = REF$size[i]/2
        }
        else {
            REF$middle[i] = sum(as.numeric(REF$size[1:(i - 1)])) + 
                REF$size[i]/2
        }
    }
    rm(i)
    REF$cumul = cumsum(as.numeric(REF$size))
    REF$add = cumsum(as.numeric(c(0, REF$size[1:(nrow(REF) - 
        1)])))
    SEGMENTS$startpos_adjusted = SEGMENTS$startpos
    SEGMENTS$endpos_adjusted = SEGMENTS$endpos
    for (CHR in unique(REF$chrom)) {
        INDEX = which(SEGMENTS$chr == CHR)
        if (length(INDEX) > 0) {
            SEGMENTS$startpos_adjusted[INDEX] = SEGMENTS$startpos_adjusted[INDEX] + 
                REF$add[which(REF$chrom == CHR)]
            SEGMENTS$endpos_adjusted[INDEX] = SEGMENTS$endpos_adjusted[INDEX] + 
                REF$add[which(REF$chrom == CHR)]
        }
        rm(INDEX)
    }
    rm(CHR)
    for (SAMPLE in sort(unique(SEGMENTS$sample))) {
        SEGS = SEGMENTS[which(SEGMENTS$sample == SAMPLE), ]
        if (nrow(SEGS) == 0) 
            warning(paste0("No segments for sample: ", 
                SAMPLE))
        maintitle = paste("Ploidy: ", sprintf("%1.2f", 
            ASCAT_output_object$ploidy[SAMPLE]), ", purity: ", 
            sprintf("%2.0f", ASCAT_output_object$purity[SAMPLE] * 
                100), "%, goodness of fit: ", sprintf("%2.1f", 
                ASCAT_output_object$goodnessOfFit[SAMPLE]), "%", 
            ifelse(isTRUE(ASCAT_output_object$nonaberrantarrays[SAMPLE]), 
                ", non-aberrant", ""), sep = "")
        print(paste0("saving plot to ", png_prefix, SAMPLE, ".adjusted", ifelse(plot_unrounded, "rawprofile", "ASCATprofile"), ".png"))
        png(filename = paste0(png_prefix, SAMPLE, ".adjusted", 
            ifelse(plot_unrounded, "rawprofile", "ASCATprofile"), 
            ".png"), width = 2000, height = (y_limit * 
            100), res = 200)
        par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, 
            cex.axis = 2.5)
        ticks = seq(0, y_limit, 1)
        plot(c(1, REF$cumul[nrow(REF)]), c(0, y_limit), type = "n", 
            xaxt = "n", yaxt = "n", main = maintitle, 
            xlab = "", ylab = "")
        axis(side = 2, at = ticks)
        abline(h = ticks, col = "lightgrey", lty = 1)
        rect(SEGS$startpos_adjusted, (SEGS$nMajor - 0.07), SEGS$endpos_adjusted, 
            (SEGS$nMajor + 0.07), col = ifelse(SEGS$nMajor >= 
                y_limit, adjustcolor(colourA, red.f = 0.75, green.f = 0.75, 
                blue.f = 0.75), colourA), border = ifelse(SEGS$nMajor >= 
                y_limit, adjustcolor(colourA, red.f = 0.75, green.f = 0.75, 
                blue.f = 0.75), colourA))
        rect(SEGS$startpos_adjusted, (SEGS$nMinor - 0.07), SEGS$endpos_adjusted, 
            (SEGS$nMinor + 0.07), col = ifelse(SEGS$nMinor >= 
                y_limit, adjustcolor(colourB, red.f = 0.75, green.f = 0.75, 
                blue.f = 0.75), colourB), border = ifelse(SEGS$nMinor >= 
                y_limit, adjustcolor(colourB, red.f = 0.75, green.f = 0.75, 
                blue.f = 0.75), colourB))
        abline(v = c(1, REF$cumul), lty = 1, col = "lightgrey")
        text(REF$middle, y_limit, REF$chrom, pos = 1, cex = 2)
        dev.off()
        rm(SEGS, ticks, maintitle)
    }
    rm(SAMPLE)
}

