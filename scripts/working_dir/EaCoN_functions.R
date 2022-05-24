source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_mini_functions.R")
source("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir/EaCoN_renorm_functions.R")
# from apt.oncoscan.2.4.0 package
tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

# Handles GZ, BZ2 or ZIP -compressed CEL files
custom_compressed_handler <- function(CELz = NULL) {
  `%do%` <- foreach::"%do%"
  CELz2 <- foreach(CEL = CELz, .combine = "c") %do% {
    tmsg(paste0("Decompressing ", CEL, " ..."))
    if (tolower(tools::file_ext(CEL)) == "bz2") {
      uncomp_file <- tempfile(fileext = ".CEL")
      R.utils::bunzip2(filename = CEL, destname = uncomp_file, FUN = bzfile, remove = FALSE)
      CEL <- uncomp_file
    } else if (tolower(tools::file_ext(CEL)) == "gz") {
      uncomp_file <- tempfile(fileext = ".CEL")
      R.utils::gunzip(filename = CEL, destname = uncomp_file, FUN = gzfile, remove = FALSE)
      CEL <- uncomp_file
    } else if (tolower(tools::file_ext(CEL)) == "zip") {
      zlist <- utils::unzip(CEL, list = TRUE)
      if (length(grep(zlist$Name, pattern = "\\.CEL", ignore.case = TRUE)) != 1) stop(tmsg(paste0(CEL, "archive file does not contain a single and unique CEL file !")), call. = FALSE)
      zname <- zlist$Name[1]
      utils::unzip(zipfile = CEL, files = zname, exdir = tempdir(), overwrite = TRUE)
      CEL <- file.path(tempdir(), zname)
    } else if (tolower(tools::file_ext(CEL)) != "cel") stop(tmsg(paste0("File ", CEL, " is not recognized as raw nor compressed (gz, bz2, zip) CEL file !")), call. = FALSE)
    return(CEL)
  }
  return(CELz2)
}


######################################### EaCoN pipeline functions

custom_OS.Process= function (ATChannelCel = NULL, GCChannelCel = NULL, samplename = NULL, 
          dual.norm = TRUE, l2r.level = "weighted", gc.renorm = TRUE, 
          gc.rda = NULL, wave.renorm = TRUE, wave.rda = NULL, mingap = 1e+06, 
          out.dir = getwd(), oschp.keep = FALSE, force.OS = NULL, apt.version = "2.4.0", 
          apt.build = "na33.r2", genome.pkg = "BSgenome.Hsapiens.UCSC.hg19", 
          return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE, oschp_file = NULL) 
{
    if (is.null(ATChannelCel)) 
        stop(tmsg("An ATChannel CEL file is required !"), 
             call. = FALSE)
    if (is.null(GCChannelCel)) 
        stop(tmsg("A GCChannel CEL file is required !"), 
             call. = FALSE)
    if (is.null(samplename)) 
        stop(tmsg("A samplename is required !"), call. = FALSE)
    if (gc.renorm) {
        if (!is.null(gc.rda)) {
            if (!file.exists(gc.rda)) 
                stop(tmsg(paste0("Could not find gc.rda file ", 
                                 gc.rda)), call. = FALSE)
        }
    }
    if (wave.renorm) {
        if (!is.null(wave.rda)) {
            if (!file.exists(wave.rda)) 
                stop(tmsg(paste0("Could not find wave.rda file ", 
                                 wave.rda)), call. = FALSE)
        }
    }
    if (!file.exists(ATChannelCel)) 
        stop(tmsg(paste0("Could not find ATChannelCel file ", 
                         ATChannelCel, " !")), call. = FALSE)
    if (!file.exists(GCChannelCel)) 
        stop(paste0("Could not find GCChannelCel file ", 
                    GCChannelCel, " !"), call. = FALSE)
    if (ATChannelCel == GCChannelCel) 
        stop(tmsg("ATChannelCel and GCChannelCel files are identical !"), 
             call. = FALSE)
    if (!genome.pkg %in% BSgenome::installed.genomes()) {
        if (genome.pkg %in% BSgenome::available.genomes()) {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), 
                 call. = FALSE)
        }
        else {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), 
                 call. = FALSE)
        }
    }
    if (dir.exists(samplename)) {
        if (!force) 
            stop(tmsg(paste0("A [", samplename, "] dir already exists !")), 
                 call. = FALSE)
        else unlink(samplename, recursive = TRUE, force = FALSE)
    }
    l2r.lev.conv <- list(normal = "Log2Ratio", weighted = "WeightedLog2Ratio")
    if (!(l2r.level %in% names(l2r.lev.conv))) 
        stop(tmsg("Option 'l2r.level' should be 'normal' or 'weighted' !"), 
             call. = FALSE)
    CEL.OS <- custom_compressed_handler(c(ATChannelCel, GCChannelCel))
    ATChannelCel <- CEL.OS[1]
    GCChannelCel <- CEL.OS[2]
    sup.array <- c("OncoScan", "OncoScan_CNV")
    arraytype.celA = affxparser::readCelHeader(filename = ATChannelCel)$chiptype
    arraytype.celC = affxparser::readCelHeader(filename = GCChannelCel)$chiptype
    if (!arraytype.celA %in% sup.array) 
        stop(tmsg(paste0("Identified array type for ATChannelCel file '", 
                         arraytype.celA, "' is not supported by this function !")), 
             call. = FALSE)
    if (!arraytype.celC %in% sup.array) 
        stop(tmsg(paste0("Identified array type for GCChannelCel file '", 
                         arraytype.celC, "' is not supported by this function !")), 
             call. = FALSE)
    if (arraytype.celA != arraytype.celC) 
        stop(tmsg(paste0("ATChannelCel and GCChannelCel files are not of the same design : ", 
                         arraytype.celA, " != ", arraytype.celC, " !")), 
             call. = FALSE)
    valid.apt.versions <- c("2.4.0")
    if (!(apt.version %in% valid.apt.versions)) 
        tmsg(paste0("APT version ", apt.version, " is not supported. Program may fail !"))
    valid.builds <- c("na33.r2", "na33.r4", "na36.r1")
    if (!(tolower(apt.build) %in% valid.builds)) 
        tmsg(paste0("Build ", apt.build, " is not supported. Program may fail !"))
    apt.onco.pkg.name <- paste0("apt.oncoscan.", apt.version)
    if (!(apt.onco.pkg.name %in% installed.packages())) 
        stop(tmsg(paste0("Package ", apt.onco.pkg.name, 
                         " not found !")), call. = FALSE)
    suppressPackageStartupMessages(require(package = apt.onco.pkg.name, 
                                           character.only = TRUE))
    if (dir.exists(samplename) && force) 
        unlink(samplename, recursive = TRUE, force = FALSE)
    ## process CEL files on linux or with APT itself since function apt.oncoscan.process() is from R package apt.oncoscan.2.4.0 which is 16-bit on windows (meaning it doesn't work on 32 or 34-bit machines)
    # oscf <- apt.oncoscan.process(ATChannelCel = ATChannelCel, 
    #                              GCChannelCel = GCChannelCel, samplename = samplename, 
    #                              dual.norm = dual.norm, out.dir = out.dir, temp.files.keep = FALSE, 
    #                              force.OS = force.OS, apt.build = apt.build)
    oscf = oschp_file
    my.oschp <- oschp.load(file = oscf)
    if (length(names(my.oschp$Meta)) == 3) 
        names(my.oschp$Meta) <- c("ATChannelCel", "GCChannelCel", 
                                  "analysis")
    sex.chr <- c("chrX", "chrY")
    if (!("affymetrix-chipsummary-snp-qc" %in% names(my.oschp$Meta$analysis))) 
        my.oschp$Meta$analysis[["affymetrix-chipsummary-snp-qc"]] <- NA
    tmsg(paste0("Loading ", genome.pkg, " ..."))
    suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
    BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
    genome2 <- metadata(BSg.obj)$genome
    cs <- chromobjector(BSg.obj)
    genome <- getmeta("affymetrix-algorithm-param-genome-version", 
                      my.oschp$Meta$analysis)
    if (genome != genome2) 
        stop(tmsg(paste0("Genome build name given with BSgenome package '", 
                         genome.pkg, "', (", genome2, ") is different from the genome build specified by provided APT build version '", 
                         apt.build, "' (", genome, ") !")), call. = FALSE)
    arraytype <- getmeta("affymetrix-array-type", my.oschp$Meta$analysis)
    manufacturer <- getmeta("program-company", my.oschp$Meta$analysis)
    species <- getmeta("affymetrix-algorithm-param-genome-species", 
                       my.oschp$Meta$analysis)
    gender.conv <- list(female = "XX", male = "XY", 
                        `NA` = "NA")
    pgender <- gender.conv[[(getmeta("affymetrix-chipsummary-Y-gender-call", 
                                     my.oschp$Meta$analysis))]]
    if (!(arraytype %in% sup.array)) 
        stop(tmsg(paste0("Unsupported array : '", arraytype, 
                         "' !")), call. = FALSE)
    meta.b <- list(samplename = samplename, source = "microarray", 
                   source.file = list(ATChannelCel = ATChannelCel, GCChannelCel = GCChannelCel), 
                   type = arraytype, manufacturer = manufacturer, species = species, 
                   genome = genome, genome.pkg = genome.pkg, predicted.gender = pgender)
    # print(c("length(my.oschp$ProbeSets$CopyNumber$ProbeSetName): ", length(my.oschp$ProbeSets$CopyNumber$ProbeSetName)))
    # print(c("length(my.oschp$AlgorithmData$MarkerABSignal$PobeSetName): ", length(my.oschp$AlgorithmData$MarkerABSignal$PobeSetName)))
    # print(c("length(my.oschp)", length(my.oschp)))
    # print(c("my.oschp[1:5]", my.oschp[1:5]))
    # return(my.oschp)
    # print(c("head( my.oschp$AlgorithmData$MarkerABSignal): ",head(my.oschp$AlgorithmData$MarkerABSignal)))
    pbn_217251 = as.vector(my.oschp$Genotyping$Calls$ProbeSetName)
    CNdf = my.oschp$ProbeSets$CopyNumber
    pbn_217611 = as.vector(CNdf$ProbeSetName)
    missingProbes = setdiff(pbn_217611, pbn_217251)
    
    # to prevent the dataframes from having different numbers of rows. 360 rows are removed because they have NULL values of log ratio.
    my.oschp$AlgorithmData$MarkerABSignal = dplyr::filter(my.oschp$AlgorithmData$MarkerABSignal, !(ProbeSetName %in% missingProbes))
    my.oschp$ProbeSets$AllelicData = dplyr::filter(my.oschp$ProbeSets$AllelicData, !(ProbeSetName %in% missingProbes))
    my.oschp$ProbeSets$CopyNumber = dplyr::filter(my.oschp$ProbeSets$CopyNumber, !(ProbeSetName %in% missingProbes))
    ### to check that NA probes were removed
    # print(c("my.oschp: ", my.oschp[1:length(my.oschp)-1]))
    # print(c("dim(my.oschp$AlgorithmData$MarkerABSignal): ", dim(my.oschp$AlgorithmData$MarkerABSignal)))
    # print(c("dim(my.oschp$Genotyping$Calls): ", dim(my.oschp$Genotyping$Calls)))
    # print(c("dim(my.oschp$ProbeSets$AllelicData): ", dim(my.oschp$ProbeSets$AllelicData)))
    # print(c("dim(my.oschp$ProbeSets$CopyNumber): ", dim(my.oschp$ProbeSets$CopyNumber)))
    
    ao.df <- data.frame(ProbeSetName = my.oschp$ProbeSets$CopyNumber$ProbeSetName, 
                        chrs = as.vector(my.oschp$ProbeSets$CopyNumber$Chromosome), 
                        pos = as.vector(my.oschp$ProbeSets$CopyNumber$Position), 
                        L2R.ori = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), 
                        L2R = as.vector(my.oschp$ProbeSets$CopyNumber[[l2r.lev.conv[[l2r.level]]]]), 
                        BAF = as.vector(my.oschp$ProbeSets$AllelicData$BAF), 
                        AD = as.vector(my.oschp$ProbeSets$AllelicData$AllelicDifference), 
                        CallF = as.vector(my.oschp$Genotyping$Calls$ForcedCall), 
                        ASignal = my.oschp$Genotyping$Calls$ASignal, BSignal = my.oschp$Genotyping$Calls$BSignal, 
                        stringsAsFactors = FALSE)
    affy.chrom <- my.oschp$Chromosomes$Summary
    ak <- affy.chrom$Display
    names(ak) <- affy.chrom$Chromosome
    ao.df$chrA <- as.vector(ak[as.character(ao.df$chrs)])
    ao.df$chr <- paste0("chr", ao.df$chrA)
    ao.df$chrN <- unlist(cs$chrom2chr[ao.df$chr])
    ao.df <- ao.df[order(ao.df$chrN, ao.df$pos, ao.df$ProbeSetName), 
    ]
    ao.df <- ao.df[!(is.na(ao.df$L2R) & is.na(ao.df$BAF)), ]
    rcmat <- round(cbind(ao.df$BAF * (ao.df$ASignal + ao.df$BSignal), 
                         (1 - ao.df$BAF) * (ao.df$ASignal + ao.df$BSignal)))
    ao.df$LOR <- log(rcmat[, 1] + 1/6) - log(rcmat[, 2] + 1/6)
    ao.df$LORvar <- 1/(rcmat[, 1] + 1/6) + 1/(rcmat[, 2] + 1/6)
    rm(rcmat)
    gc()
    smo <- round(nrow(ao.df)/550)
    if (smo%%2 == 0) 
        smo <- smo + 1
    if (wave.renorm) {
        tmsg("Wave re-normalization ...")
        ren.res <- renorm.go(input.data = ao.df, renorm.rda = wave.rda, 
                             track.type = "Wave", smo = smo, arraytype = arraytype, 
                             genome = genome)
        ao.df <- ren.res$data
        fitted.l2r <- ren.res$renorm$l2r$l2r
        if (is.null(ren.res$renorm$pos)) {
            meta.b <- setmeta("wave.renorm", "None", 
                              meta.b)
            tmsg(" No positive fit.")
        }
        else {
            sex.idx <- ao.df$chr %in% sex.chr
            auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
            auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
            if (any(sex.idx)) {
                for (k in sex.chr) {
                    k.idx <- ao.df$chr == k
                    if (any(k.idx)) {
                        k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], 
                                                na.rm = TRUE) - auto.ori.med
                        k.rn.diffmed <- median(fitted.l2r[k.idx], 
                                               na.rm = TRUE) - auto.rn.med
                        fitted.l2r[k.idx] <- fitted.l2r[k.idx] - 
                            k.rn.diffmed + k.ori.diffmed
                    }
                }
            }
            meta.b <- setmeta("wave.renorm", paste0(ren.res$mrenorm$pos, 
                                                    collapse = ","), meta.b)
        }
        rm(ren.res)
        ao.df[["L2R.Wave"]] <- fitted.l2r - median(fitted.l2r, 
                                                   na.rm = TRUE)
        ao.df$L2R <- ao.df[["L2R.Wave"]]
    }
    else {
        meta.b <- setmeta("wave.renorm", "FALSE", 
                          meta.b)
    }
    if (gc.renorm) {
        tmsg("GC renormalization ...")
        ren.res <- renorm.go(input.data = ao.df, renorm.rda = gc.rda, 
                             track.type = "GC", smo = smo, arraytype = arraytype, 
                             genome = genome)
        ao.df <- ren.res$data
        fitted.l2r <- ren.res$renorm$l2r$l2r
        if (is.null(ren.res$renorm$pos)) {
            meta.b <- setmeta("gc.renorm", "None", 
                              meta.b)
            message(tmsg(" No positive fit."))
        }
        else {
            sex.idx <- ao.df$chr %in% sex.chr
            auto.ori.med <- median(ao.df$L2R[!sex.idx], na.rm = TRUE)
            auto.rn.med <- median(fitted.l2r[!sex.idx], na.rm = TRUE)
            if (any(sex.idx)) {
                for (k in sex.chr) {
                    k.idx <- ao.df$chr == k
                    if (any(k.idx)) {
                        k.ori.diffmed <- median(ao.df$L2R.ori[k.idx], 
                                                na.rm = TRUE) - auto.ori.med
                        k.rn.diffmed <- median(fitted.l2r[k.idx], 
                                               na.rm = TRUE) - auto.rn.med
                        fitted.l2r[k.idx] <- fitted.l2r[k.idx] - 
                            k.rn.diffmed + k.ori.diffmed
                    }
                }
            }
            meta.b <- setmeta("gc.renorm", paste0(ren.res$renorm$pos, 
                                                  collapse = ","), meta.b)
        }
        rm(ren.res)
        ao.df[["L2R.GC"]] <- fitted.l2r - median(fitted.l2r, 
                                                 na.rm = TRUE)
        ao.df$L2R <- ao.df[["L2R.GC"]]
    }
    else {
        meta.b <- setmeta("gc.renorm", "FALSE", meta.b)
    }
    ao.df$L2R <- ao.df$L2R - median(ao.df$L2R, na.rm = TRUE)
    ao.df$AD[is.nan(ao.df$AD)] <- NA
    gaps <- which(diff(ao.df$pos) >= mingap)
    kends <- vapply(unique(ao.df$chrs), function(k) {
        max(which(ao.df$chrs == k))
    }, 1L)
    kbreaks <- sort(unique(c(gaps, kends)))
    ao.df$chrgap <- rep(seq_along(kbreaks), times = c(kbreaks[1], 
                                                      diff(kbreaks)))
    germ <- ao.df$CallF
    if (as.numeric(version$major) >= 4) {
        germ[germ %in% as.raw(c(8, 11))] <- as.raw(0)
        germ[germ != 0] <- as.raw(1)
    }
    else {
        germ[germ %in% c(8, 11)] <- 0
        germ[germ != 0] <- 1
    }
    tmsg("Building normalized object ...")
    my.ch <- sapply(unique(ao.df$chrs), function(x) {
        which(ao.df$chrs == x)
    })
    my.ascat.obj <- list(data = list(Tumor_LogR.ori = data.frame(sample = ao.df$L2R.ori, 
                                                                 row.names = ao.df$ProbeSetName), Tumor_LogR = data.frame(sample = ao.df$L2R, 
                                                                                                                          row.names = ao.df$ProbeSetName), Tumor_BAF = data.frame(sample = ao.df$BAF, 
                                                                                                                                                                                  row.names = ao.df$ProbeSetName), Tumor_AD = data.frame(sample = ao.df$AD, 
                                                                                                                                                                                                                                         row.names = ao.df$ProbeSetName), Tumor_LogR_segmented = NULL, 
                                     Tumor_BAF_segmented = NULL, Germline_LogR = NULL, Germline_BAF = NULL, 
                                     SNPpos = data.frame(chrs = ao.df$chr, pos = ao.df$pos, 
                                                         row.names = ao.df$ProbeSetName), ch = sapply(unique(ao.df$chr), 
                                                                                                      function(x) {
                                                                                                          which(ao.df$chr == x)
                                                                                                      }), chr = sapply(unique(ao.df$chrgap), function(x) {
                                                                                                          which(ao.df$chrgap == x)
                                                                                                      }), chrs = unique(ao.df$chr), samples = samplename, gender = as.vector(meta.b$predicted.gender), 
                                     sexchromosomes = sex.chr, failedarrays = NULL, additional = data.frame(RD.test = ao.df$ASignal, 
                                                                                                            RD.ref = ao.df$BSignal, LOR = ao.df$AD, LORvar = ao.df$LORvar, 
                                                                                                            stringsAsFactors = FALSE)), meta = list(basic = meta.b, 
                                                                                                                                                    affy = my.oschp$Meta), germline = list(germlinegenotypes = matrix(as.logical(germ), 
                                                                                                                                                                                                                      ncol = 1), failedarrays = NULL), CEL = list(ATChannelCel = affxparser::readCel(filename = ATChannelCel), 
                                                                                                                                                                                                                                                                  GCChannelCel = affxparser::readCel(filename = GCChannelCel)))
    colnames(my.ascat.obj$germline$germlinegenotypes) <- colnames(my.ascat.obj$data$Tumor_LogR) <- colnames(my.ascat.obj$data$Tumor_LogR.ori) <- colnames(my.ascat.obj$data$Tumor_BAF) <- samplename
    rownames(my.ascat.obj$germline$germlinegenotypes) <- ao.df$ProbeSetName
    my.ascat.obj$CEL$ATChannelCel$intensities <- as.integer(my.ascat.obj$CEL$ATChannelCel$intensities)
    my.ascat.obj$CEL$GCChannelCel$intensities <- as.integer(my.ascat.obj$CEL$GCChannelCel$intensities)
    gc()
    if (write.data) 
        print("attempting to save to RDS file...")
        # saveDir = paste0(out.dir, "/", samplename)
        RDSfilename = paste0(samplename, "_", arraytype, "_", genome, "_processed.RDS")
        
        # print(c("saving path: ", file.path(saveDir,RDSfilename)))
        saveRDS(my.ascat.obj, RDSfilename, compress = "bzip2")
    genopos <- ao.df$pos + cs$chromosomes$chr.length.toadd[ao.df$chrN]
    rm(ao.df)
    gc()
    if (plot) {
        tmsg("Plotting ...")
        kend <- genopos[vapply(my.ascat.obj$data$ch, max, 1L)]
        l2r.notna <- which(!is.na(my.ascat.obj$data$Tumor_LogR[, 
                                                               1]))
        l2r.rm <- runmed(my.ascat.obj$data$Tumor_LogR[, 1][l2r.notna], 
                         smo)
        l2r.ori.rm <- runmed(my.ascat.obj$data$Tumor_LogR.ori[, 
                                                              1][l2r.notna], smo)
        ### to check plot path
        folder_name = paste0(out.dir, "/", samplename)
        if(!dir.exists(folder_name)) dir.create(folder_name)
        print(paste0("plot path: ", out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"))
        png(paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_rawplot.png"), 1600, 1050)
        par(mfrow = c(3, 1))
        plot(genopos, my.ascat.obj$data$Tumor_LogR.ori[, 1], 
             pch = ".", cex = 3, col = "grey70", xaxs = "i", 
             yaxs = "i", ylim = c(-2, 2), main = paste0(samplename, 
                                                        " ", arraytype, " raw L2R profile (median-centered) / ", 
                                                        round(sum(abs(diff(l2r.ori.rm))), digits = 3)), 
             xlab = "Genomic position", ylab = "L2R")
        lines(genopos[l2r.notna], l2r.ori.rm, col = 1)
        abline(v = kend, col = 4, lty = 3, lwd = 2)
        abline(h = 0, col = 2, lty = 2, lwd = 2)
        plot(genopos, my.ascat.obj$data$Tumor_LogR[, 1], pch = ".", 
             cex = 3, col = "grey70", xaxs = "i", 
             yaxs = "i", ylim = c(-2, 2), main = paste0(samplename, 
                                                        " ", arraytype, " L2R profile (", 
                                                        l2r.level, ", median-centered)) / ", round(sum(abs(diff(l2r.rm))), digits = 3)), xlab = "Genomic position", ylab = "L2R")
        lines(genopos[l2r.notna], l2r.rm, col = 1)
        abline(v = kend, col = 4, lty = 3, lwd = 2)
        abline(h = 0, col = 2, lty = 2, lwd = 2)
        plot(genopos, my.ascat.obj$data$Tumor_BAF[, 1], pch = ".", 
             cex = 3, col = "grey75", xaxs = "i", 
             yaxs = "i", ylim = c(0, 1), main = paste0(samplename, 
                                                       " ", arraytype, " BAF profile"), 
             xlab = "Genomic position", ylab = "BAF")
        points(genopos[germ == 0], my.ascat.obj$data$Tumor_BAF[germ == 
                                                                   0, 1], pch = ".", cex = 3, col = 2)
        abline(v = kend, col = 4, lty = 3, lwd = 2)
        abline(h = 0.5, col = 2, lty = 2, lwd = 2)
        dev.off()
    }
    if (!oschp.keep) {
        tmsg("Removing temporary OSCHP file ...")
        file.remove(oscf)
    }
    message(tmsg("Done."))
    gc()
    if (return.data) 
        return(my.ascat.obj)
}









custom_Segment.ff = function (RDS.file = NULL, segmenter = "ASCAT", ASCATobj=NULL) 
{
    # if (is.null(RDS.file)) 
    #     stop(tmsg("A RDS file is needed !"), call. = FALSE)
    valid.segmenters <- c("ASCAT", "FACETS", "SEQUENZA")
    if (!(toupper(segmenter) %in% valid.segmenters)) 
        stop(tmsg(paste0("Segmenter should be one of : ", 
            paste0(valid.segmenters, collapse = ", "))), 
            call. = FALSE)
    # if (!file.exists(RDS.file)) 
    #     stop(tmsg(paste0("Could not find RDS file ", RDS.file, 
    #         " !")), call. = FALSE)
    # tmsg(paste0("Loading data from ", RDS.file, " ..."))
    # my.data <- readRDS(RDS.file)
    if (is.null(ASCATobj)) 
        stop(tmsg("An ASCAT object is needed !"), call. = FALSE)
    print("ASCAT object will be used, not RDS file.")
    my.data = ASCATobj
    do.call(paste0("Segment.", toupper(segmenter)), list(data = my.data, 
        out.dir = dirname(RDS.file)))
}









Segment.ASCAT = function (data = NULL, mingap = 5e+06, smooth.k = NULL, BAF.filter = 0.75, 
    homoCut = 0.05, penalty = 50, recenter = "l2r.centeredpeak", 
    calling.method = "mad", nrf = 0.5, SER.pen = 40, out.dir = getwd(), 
    return.data = FALSE, write.data = TRUE, plot = TRUE, force = FALSE) 
{
    `%do%` <- foreach::"%do%"
    calling.method <- tolower(calling.method)
    if (!is.list(data)) 
        stop(tmsg("data should be a list !"), call. = FALSE)
    if (!dir.exists(out.dir)) 
        stop(tmsg(paste0("Output directory [", out.dir, 
            "] does not exist !")), call. = FALSE)
    if (!(calling.method %in% c("mad", "density"))) 
        stop(tmsg("calling.method should be 'MAD' or 'density' !"), 
            call. = FALSE)
    if (calling.method == "mad" & is.null(nrf)) 
        stop(tmsg("If calling.method is set to 'MAD', nrf is required !"), 
            call. = FALSE)
    if (!is.null(SER.pen)) 
        if (!is.character(SER.pen)) 
            if (SER.pen <= 0) 
                stop(tmsg("SER.pen should be NULL, a character or an integer/float > 0 !"), 
                  call. = FALSE)
    samplename <- data$meta$basic$samplename
    tmsg(paste0("Sample : ", samplename))
    genome <- data$meta$basic$genome
    genome.pkg <- data$meta$basic$genome.pkg
    if (!genome.pkg %in% BSgenome::installed.genomes()) {
        if (genome.pkg %in% BSgenome::available.genomes()) {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), 
                call. = FALSE)
        }
        else {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), 
                call. = FALSE)
        }
    }
    tmsg(paste0("Loading ", genome.pkg, " ..."))
    suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
    BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
    cs <- chromobjector(BSg.obj)
    oridir <- getwd()
    if (write.data | plot) {
        odir <- paste0(out.dir, "/ASCAT/L2R")
        if (dir.exists(odir)) {
            if (force) {
                unlink(odir, recursive = TRUE, force = FALSE)
                dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
            }
            else stop(tmsg(paste0("A [", odir, "] dir already exists !")), 
                call. = FALSE)
        }
        else dir.create(path = odir, recursive = TRUE, showWarnings = FALSE)
    }
    else odir <- out.dir
    setwd(odir)
    data$meta$eacon <- c(data$meta$eacon, list(segmenter = "ASCAT", 
        mingap = mingap, BAF.filter = BAF.filter, BAF.segments.homo.limit = homoCut, 
        winsorize.k = if (is.null(smooth.k)) "NA" else smooth.k, 
        ASCAT.penalty = penalty, calling.method = calling.method, 
        calling.nrf = if (is.null(nrf)) "NA" else nrf, 
        small.events.rescue.PELT.penalty = if (is.null(SER.pen)) "NA" else SER.pen))
    if (!is.null(smooth.k)) {
        tmsg("Smoothing L2R outliers ...")
        cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), 
            vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), 
            MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
        l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
        cndf <- cndf[l2r.nona, ]
        cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", 
            method = "mad", k = smooth.k, tau = 1, verbose = FALSE)
        data$data$Tumor_LogR[l2r.nona, 1] <- cndf.wins[, 3, drop = FALSE]
        rm(list = c("cndf", "cndf.wins", "l2r.nona"))
    }
    tmsg("Filtering BAF...")
    if ("Tumor_BAF.unisomy" %in% names(data$data)) {
        called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & 
            !is.na(data$data$Tumor_BAF[, 1]) & !data$data$Tumor_BAF.unisomy[, 
            1])
    }
    else {
        called <- which(!data$germline$germlinegenotypes & !is.na(data$germline$germlinegenotypes) & 
            !is.na(data$data$Tumor_BAF[, 1]))
    }
    mBAF <- BAF2mBAF(data$data$Tumor_BAF[, 1])
    smoB <- round(length(called)/3300)
    if (smoB%%2 == 0) 
        smoB <- smoB + 1
    mBAF.rm <- runmed(mBAF[called], smoB)
    mBAF.diff <- abs(mBAF[called] - mBAF.rm)
    Bfiltered <- mBAF.diff <= quantile(mBAF.diff, BAF.filter)
    if (any(Bfiltered)) 
        data$germline$germlinegenotypes[called][!Bfiltered] <- TRUE
    rm(called, mBAF, smoB, mBAF.rm, mBAF.diff, Bfiltered)
    # if (plot) {
    #     png(paste0(samplename, ".Rorschach.png"), width = 980, 
    #         height = 980)
    #     EaCoN.Rorschard.plot(data = data)
    #     dev.off()
    # }
    if (!is.null(mingap)) {
        data$data$chr <- foreach::foreach(k = data$data$ch, .combine = "c") %do% 
            {
                gapz <- which(diff(data$data$SNPpos$pos[k]) >= 
                  mingap)
                return(unname(split(k, findInterval(k, k[gapz + 
                  1]))))
            }
    }
    tmsg("ASPCF segmentation ...")
    aspcf.res <- ASCAT::ascat.aspcf(ASCATobj = data$data, ascat.gg = data$germline, 
        penalty = penalty)
    aspcf.res$germline <- data$germline
    data$data <- aspcf.res
    rm(aspcf.res)
    pcftxt <- list.files(path = ".", pattern = "\\.PCFed\\.txt", 
        full.names = TRUE, recursive = FALSE, ignore.case = TRUE, 
        include.dirs = FALSE)
    if (length(pcftxt) > 0) 
        unlink(pcftxt)
    data$data$Tumor_BAF_segmented[[1]][(data$data$Tumor_BAF_segmented[[1]][, 
        1] < 0), 1] <- 0
    smo <- round(length(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[, 
        1]), 1])/550)
    if (smo%%2 == 0) 
        smo <- smo + 1
    if (recenter %in% c("l2r.mainpeak", "l2r.centeredpeak", 
        "l2r.median") | is.numeric(recenter)) {
        tmsg("Recentering ...")
        if (is.numeric(recenter)) {
            tmsg(paste0(" ... using direct value ", recenter, 
                "."))
            shifter <- recenter
        }
        else if (recenter == "l2r.median") {
            tmsg(" ... using median.")
            shifter <- stats::median(data$data$Tumor_LogR[, 1], 
                na.rm = TRUE)
        }
        else {
            lrrm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[, 
                1]), 1], smo)
            my.den <- stats::density(lrrm, bw = 0.015)
            if (recenter == "l2r.mainpeak") {
                tmsg(" ... using the main peak.")
                shifter <- my.den$x[which.max(my.den$y)]
            }
            else if (recenter == "l2r.centeredpeak") {
                tmsg(" ... using the most centered of populated peaks.")
                denf <- data.frame(x = my.den$x, y = my.den$y, 
                  sign = c(1, sign(diff(my.den$y))))
                denf$sign2 <- c(diff(denf$sign), 0)
                rrr <- rle(denf$sign2)
                repr <- data.frame(values = rrr$values, start = c(0, 
                  (cumsum(rrr$lengths[-c(length(rrr$lengths))]) + 
                    1)), end = cumsum(rrr$lengths), stringsAsFactors = F)
                npk <- which(repr$values == -2)
                if (length(npk) == 1) {
                  fx <- my.den$x[repr$start[npk[1]]]
                }
                else {
                  if (1 %in% npk) 
                    npk <- npk[-1]
                  if (nrow(repr) %in% npk) 
                    npk <- npk[nrow(repr)]
                  parea <- sapply(npk, function(r) {
                    sum(my.den$y[repr$start[r - 1]:repr$end[r + 
                      1]])/sum(my.den$y)
                  })
                  parea.rel <- parea/max(parea)
                  npk <- npk[parea > 0.1]
                  peaklist <- repr$start[npk]
                  shifter <- denf$x[peaklist[which.min(sapply(peaklist, 
                    function(p) {
                      abs(sum(denf$y[1:(p - 1)]) - sum(denf$y[(p + 
                        1):nrow(denf)]))/sum(denf$y)
                    }))]]
                }
            }
        }
        data$data$Tumor_LogR[, 1] <- data$data$Tumor_LogR[, 1] - 
            shifter
        data$data$Tumor_LogR_segmented <- data$data$Tumor_LogR_segmented - 
            shifter
        data$meta$eacon[["recenter-value"]] <- shifter
    }
    else if (is.null(recenter)) {
        tmsg("No recentering.")
    }
    else stop(tmsg("Invalid recentering method called !"), 
        call. = FALSE)
    tmsg("Smoothing L2R (for plots)...")
    cndf <- data.frame(Chr = rep(unlist(cs$chrom2chr[data$data$chrs]), 
        vapply(data$data$ch, length, 1L)), Position = unlist(data$data$ch), 
        MySample = data$data$Tumor_LogR[[1]], stringsAsFactors = FALSE)
    l2r.nona <- !is.na(data$data$Tumor_LogR[[1]])
    cndf <- cndf[l2r.nona, ]
    cndf.wins <- copynumber::winsorize(data = cndf, pos.unit = "bp", 
        method = "mad", k = 5, tau = 1, verbose = FALSE)
    data$data$Tumor_LogR_wins <- data$data$Tumor_LogR
    data$data$Tumor_LogR_wins[l2r.nona, ] <- cndf.wins[, 3, drop = FALSE]
    colnames(data$data$Tumor_LogR_wins) <- samplename
    rm(list = c("cndf", "cndf.wins", "l2r.nona"))
    if (!is.null(SER.pen)) {
        tmsg("Rescuing small events ...")
        seg.maxwidth <- 5e+06
        seg.maxn <- 500
        mydf <- data.frame(data$data$SNPpos, l2r = as.vector(data$data$Tumor_LogR[, 
            1]), idx.ori = 1:nrow(data$data$SNPpos))
        mydf$chrs <- as.character(mydf$chrs)
        mydf <- mydf[!is.na(mydf$l2r), ]
        chrends <- cumsum(rle(as.character(data$data$SNPpos$chrs[!is.na(data$data$Tumor_LogR[, 
            1])]))$lengths)
        if (is.character(SER.pen)) {
            seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, 
                penalty = SER.pen, method = "PELT", param.estimates = FALSE, 
                minseglen = 5)@cpts))
        }
        else if (is.numeric(SER.pen)) {
            if (SER.pen < 1) {
                seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, 
                  penalty = "Asymptotic", pen.value = SER.pen, 
                  method = "PELT", param.estimates = FALSE, 
                  minseglen = 5)@cpts))
            }
            else {
                SER.pen <- sum(abs(diff(runmed(mydf$l2r[!is.na(mydf$l2r)], 
                  smo))))/SER.pen
                seg.end <- try(suppressWarnings(changepoint::cpt.mean(data = mydf$l2r, 
                  penalty = "Manual", pen.value = SER.pen, 
                  method = "PELT", param.estimates = FALSE, 
                  minseglen = 5)@cpts))
            }
        }
        else stop(tmsg("SER.pen should be a character or a numeric !"), 
            call. = FALSE)
        if (is.character(seg.end)) {
            tmsg(" PELT segmentation failed with this combination of SER.pen and segmentLength options !")
            data$meta$eacon[["small.events.rescue.PELT.penalty"]] <- "ERROR"
        }
        else {
            seg.end <- sort(unique(c(seg.end, chrends)))
            seg.start <- c(1, seg.end[-length(seg.end)] + 1)
            seg.med <- vapply(1:length(seg.end), function(x) {
                return(median(mydf$l2r[seg.start[x]:seg.end[x]], 
                  na.rm = TRUE))
            }, 0.1)
            seg.width <- mydf$pos[seg.end] - mydf$pos[seg.start] + 
                1
            rescued <- which(seg.width < seg.maxwidth)
            tmsg(paste0(" Found ", length(rescued), "."))
            if (length(rescued) > seg.maxn) 
                tmsg("WARNING : Many small events found, profile may be noisy ! Consider using 'smooth.k', or for WES data, strengthen low depth filtering !")
            data$meta$eacon[["PELT-nseg"]] <- length(rescued)
            foreach::foreach(re = rescued, .combine = "c") %do% 
                {
                  interv <- mydf$idx.ori[seg.start[re]]:mydf$idx.ori[seg.end[re]]
                  data$data$Tumor_LogR_segmented[interv] <- median(data$data$Tumor_LogR[interv, 
                    1], na.rm = TRUE)
                  return()
                }
        }
    }
    bafpos <- data$data$SNPpos[rownames(data$data$Tumor_LogR_segmented) %in% 
        rownames(data$data$Tumor_BAF_segmented[[1]]), ]
    bafpos$ProbeSet <- rownames(data$data$Tumor_BAF_segmented[[1]])
    bafpos$BAF <- data$data$Tumor_BAF_segmented[[1]][, 1]
    bafpos$chrs <- as.character(bafpos$chrs)
    bafkend <- vapply(unique(bafpos$chrs), function(x) {
        max(which(bafpos$chrs == x))
    }, 1)
    bafbreaks <- cumsum(rle(bafpos$BAF)$lengths)
    baf.seg <- data.frame(end.idx = sort(unique(c(bafkend, bafbreaks))), 
        stringsAsFactors = FALSE)
    baf.seg$start.idx <- c(1, baf.seg$end.idx[-nrow(baf.seg)] + 
        1)
    baf.seg$length.idx <- baf.seg$end.idx - baf.seg$start.idx + 
        1
    baf.seg$start.probeset <- bafpos$ProbeSet[baf.seg$start.idx]
    baf.seg$end.probeset <- bafpos$ProbeSet[baf.seg$end.idx]
    baf.seg$chrA <- bafpos$chrs[baf.seg$start.idx]
    baf.seg$Chr <- unlist(cs$chrom2chr[baf.seg$chrA])
    baf.seg$Start <- bafpos$pos[baf.seg$start.idx]
    baf.seg$End <- bafpos$pos[baf.seg$end.idx]
    baf.seg$Width <- baf.seg$End - baf.seg$Start + 1
    baf.seg$Value <- bafpos$BAF[baf.seg$start.idx]
    baf.homocut <- homoCut
    baf.seg$Status <- "Unbalanced"
    baf.homo <- sort(which(baf.seg$Value <= baf.homocut))
    if (length(baf.homo) > 0) 
        baf.seg$Status[baf.homo] <- "Homo"
    baf.hetero <- sort(which(baf.seg$Value == 0.5))
    if (length(baf.hetero) > 0) 
        baf.seg$Status[baf.hetero] <- "Hetero"
    baf.seg.out <- baf.seg[, c(6, 7:12, 4:5, 1:3)]
    colnames(baf.seg.out) <- c("Chrom", "Chr", "Start", 
        "End", "Width", "BAF.Value", "BAF.Status", 
        "Start.FeatureName", "End.FeatureName", "Start.FeatureIdx", 
        "End.FeatureIdx", "Features")
    if (write.data) 
        write.table(baf.seg.out, paste0(samplename, ".SegmentedBAF.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
    tmsg("Calling L2R ...")
    l2r.segments <- foreach::foreach(k = 1:length(data$data$ch), 
        .combine = "rbind") %do% {
        segrle <- rle(data$data$Tumor_LogR_segmented[data$data$ch[[k]], 
            1])
        segdf <- data.frame(Probes = segrle$lengths, Value = segrle$values, 
            stringsAsFactors = FALSE)
        segdf$Chr <- unlist(cs$chrom2chr[data$data$chrs[k]])
        segdf$Chrom <- data$data$chrs[[k]]
        probecs <- cumsum(segdf$Probes)
        segdf$End.idx <- as.integer(cumsum(segdf$Probes) + data$data$ch[[k]][1] - 
            1)
        segdf$Start.idx <- as.integer(c(data$data$ch[[k]][1], 
            segdf$End.idx[-nrow(segdf)] + 1))
        segdf$Start <- as.integer(data$data$SNPpos$pos[segdf$Start.idx])
        segdf$End <- as.integer(data$data$SNPpos$pos[segdf$End.idx])
        segdf <- segdf[, c(3, 4, 7, 8, 2, 1, 6, 5)]
    }
    my.mad <- get.mad(data$data$Tumor_LogR[, 1])
    if (calling.method == "mad") {
        g.cut <- my.mad * nrf
        l.cut <- -g.cut
    }
    l2r.rm <- stats::runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[, 
        1]), 1], smo)
    if (calling.method == "density") {
        lr.den <- stats::density(l2r.rm, bw = 0.015)
        atf <- lr.den$y > max(lr.den$y)/20
        newx <- lr.den$x[atf]
        newy <- lr.den$y[atf]
        sigdiff <- sign(diff(newy))
        negz <- newx < 0
        l.idx <- length(which(negz)) - rle(rev(sign(diff(newy[negz]))))$lengths[1] + 
            1
        g.idx <- length(which(negz)) + rle(sign(diff(newy[!negz])))$lengths[1]
        g.cut <- newx[g.idx]
        l.cut <- newx[l.idx]
        if (abs(g.cut) > abs(l.cut)) 
            l.cut <- -g.cut
        else g.cut <- -l.cut
    }
    data$meta$eacon[["L2R-segments-gain-cutoff"]] <- g.cut
    data$meta$eacon[["L2R-segments-loss-cutoff"]] <- l.cut
    data$meta$eacon$MAD <- my.mad
    data$meta$eacon$SSAD <- sum(abs(diff(runmed(data$data$Tumor_LogR[!is.na(data$data$Tumor_LogR[, 
        1]), 1], smo))))
    gain.idx <- which(l2r.segments$Value > g.cut)
    loss.idx <- which(l2r.segments$Value < l.cut)
    normal.idx <- which(l2r.segments$Value >= l.cut & l2r.segments$Value <= 
        g.cut)
    tmsg("Writing CBS files ...")
    data$cbs$nocut <- data.frame(Samplename = samplename, l2r.segments[, 
        c(1, 3, 4, 6, 5)], stringsAsFactors = FALSE)
    colnames(data$cbs$nocut) <- c(samplename, "Chr", "Start", 
        "End", "Probes", "Log2Ratio")
    my.cbs.cut <- data$cbs$nocut
    my.cbs.cut$Log2Ratio[my.cbs.cut$Log2Ratio > l.cut & my.cbs.cut$Log2Ratio < 
        g.cut] <- 0
    data$cbs$cut <- my.cbs.cut
    rm(my.cbs.cut)
    if (write.data) {
        write.table(data$cbs$nocut, paste0(samplename, ".NoCut.cbs"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
        write.table(data$cbs$cut, paste0(samplename, ".Cut.cbs"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (plot) {
        tmsg("Plotting ...")
        l2r.seg.obj <- list(pos = l2r.segments, idx = list(gain = gain.idx, 
            loss = loss.idx, normal = normal.idx), cutval = c(l.cut, 
            g.cut))
        seg.col <- list(gain = "blue", outscale.gain = "midnightblue", 
            loss = "red", outscale.loss = "darkred", 
            normal = "black")
        l2r.chr <- unname(unlist(cs$chrom2chr[as.character(data$data$SNPpos$chrs)]))
        l2r.value <- data.frame(Chr = l2r.chr, Start = as.integer(data$data$SNPpos$pos), 
            End = as.integer(data$data$SNPpos$pos), Value = data$data$Tumor_LogR_wins[, 
                1], stringsAsFactors = FALSE)
        baf.value <- data.frame(Chr = l2r.chr, Start = as.integer(data$data$SNPpos$pos), 
            End = as.integer(data$data$SNPpos$pos), Value = data$data$Tumor_BAF[, 
                1], stringsAsFactors = FALSE)
        png(paste0(samplename, ".SEG.ASCAT.png"), width = 1850, 
            height = 980)
        par(mar = c(1, 6, 3, 1), mfrow = c(2, 1))
        EaCoN.l2rplot.geno(l2r = l2r.value, seg = l2r.seg.obj, 
            seg.col = seg.col, seg.type = "block", seg.normal = TRUE, 
            genome.pkg = genome.pkg, title = paste0(samplename, 
                " L2R"), ylim = c(-1.5, 1.5))
        EaCoN.bafplot.geno(baf = baf.value, seg = baf.seg, seg.type = "both", 
            genome.pkg = genome.pkg, title = paste0(samplename, 
                " BAF"))
        dev.off()
    }
    if (write.data) {
        tmsg("Saving data ...")
        saveRDS(data, paste0(samplename, ".SEG.ASCAT.RDS"), 
            compress = "bzip2")
    }
    setwd(oridir)
    tmsg("Done.")
    if (return.data) 
        return(data)
}












custom_ASCN.ff = function (ASCATobj=NULL) 
{
    # if (is.null(RDS.file)) 
    #     stop(tmsg("A RDS file is needed !"), call. = FALSE)
    # if (!file.exists(RDS.file)) 
    #     stop(tmsg(paste0("Could not find RDS file ", RDS.file, 
    #         " !")), call. = FALSE)
    # tmsg(paste0("Loading data from ", RDS.file, " ..."))
    # my.data <- readRDS(RDS.file)
    my.data <- ASCATobj
    segmenter <- toupper(my.data$meta$eacon$segmenter)
    tmsg(paste0("Provided data was segmented with ", 
        my.data$meta$eacon$segmenter, "."))
    out.dir <- sub(pattern = paste0(toupper(segmenter), "/L2R"), 
        replacement = "", x = dirname(RDS.file))
    if (out.dir == "") 
        out.dir <- "."
    do.call(paste0("ASCN.", toupper(segmenter)), list(data = my.data, 
        out.dir = out.dir, ...))
}









custom_ASCN.ASCAT =  function (data = NULL, gammaRange = c(0.35, 0.95), nsubthread = 1, 
    cluster.type = "PSOCK", out.dir = getwd(), force = FALSE, 
    ...) 
{
    if (!is.list(data)) 
        stop(tmsg("data should be a list !"), call. = FALSE)
    odir <- paste0(out.dir, "/ASCAT/ASCN")
    if (any(is.null(c(data$data$Tumor_LogR_segmented, data$data$Tumor_BAF_segmented)))) 
        stop(tmsg("No segmentation data found in the provided RDS file !"), 
            call. = FALSE)
    if (dir.exists(odir)) {
        if (force) {
            unlink(odir, recursive = TRUE, force = FALSE)
        }
        else stop(tmsg(paste0("Output directory [", out.dir, 
            "] already exists !")), call. = FALSE)
    }
    samplename <- data$meta$basic$samplename
    genome <- data$meta$basic$genome
    genome.pkg <- data$meta$basic$genome.pkg
    if (!genome.pkg %in% BSgenome::installed.genomes()) {
        if (genome.pkg %in% BSgenome::available.genomes()) {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " available but not installed. Please install it !")), 
                call. = FALSE)
        }
        else {
            stop(tmsg(paste0("BSgenome ", genome.pkg, " not available in valid BSgenomes and not installed ... Please check your genome name or install your custom BSgenome !")), 
                call. = FALSE)
        }
    }
    tmsg(paste0("Loading ", genome.pkg, " ..."))
    suppressPackageStartupMessages(require(genome.pkg, character.only = TRUE))
    BSg.obj <- getExportedValue(genome.pkg, genome.pkg)
    cs <- chromobjector(BSg.obj)
    tmsg("ASCN modeling (using ASCAT) ...")
    gammavec <- if (length(gammaRange) > 1) 
        seq(gammaRange[1], gammaRange[2], 0.05)
    else gammaRange
    oridirx <- getwd()
    cls <- parallel::makeCluster(spec = nsubthread, type = cluster.type, 
        outfile = "")
    doParallel::registerDoParallel(cls)
    gamma <- 0
    `%dopar%` <- foreach::"%dopar%"
    fit.val <- as.data.frame(foreach::foreach(gamma = gammavec, 
        .combine = "rbind", .inorder = TRUE) %dopar% {
        tmsg(paste0(" gamma = ", gamma))
        odirg <- paste0(odir, "/gamma", sprintf("%.2f", 
            gamma))
        dir.create(path = odirg, recursive = TRUE, showWarnings = FALSE)
        setwd(odirg)
        my.ascat.seg.ascn <- suppressWarnings(ASCAT::ascat.runAscat(ASCATobj = data$data, 
            gamma = gamma, ...))
        if (is.null(my.ascat.seg.ascn$nA)) {
            tmsg("  ASCAT could not find an optimal ploidy / cellularity from the data.")
            setwd(oridirx)
            unlink(odirg, recursive = TRUE)
            return(rep(NA, 8))
        }
        else {
            unlink(paste0(samplename, ".ASCATprofile.png"))
            unlink(paste0(samplename, ".sunrise.png"))
            tcn.tbl.ung <- dplyr::as.tbl(cbind(my.ascat.seg.ascn$segments, 
                nTotal = my.ascat.seg.ascn$segments$nMajor + 
                  my.ascat.seg.ascn$segments$nMinor, width = my.ascat.seg.ascn$segments$endpos - 
                  my.ascat.seg.ascn$segments$startpos + 1))
            tcn.tbl <- dplyr::group_by(tcn.tbl.ung, nTotal)
            tcn.tbl.prop <- dplyr::summarise(tcn.tbl, tot_width = sum(width))
            ascat.ploidy <- my.ascat.seg.ascn$ploidy
            median.ploidy <- limma::weighted.median(tcn.tbl.prop$nTotal, 
                tcn.tbl.prop$tot_width)
            baseline.ploidy <- tcn.tbl.prop$nTotal[which.max(tcn.tbl.prop$tot_width)]
            if (baseline.ploidy == 0) 
                baseline.ploidy <- tcn.tbl.prop$nTotal[order(tcn.tbl.prop$tot_width, 
                  decreasing = TRUE)[2]]
            weighted.ploidy <- sum(tcn.tbl.prop$nTotal * (tcn.tbl.prop$tot_width/sum(tcn.tbl.prop$tot_width)))
            my.ascat.seg.ascn$ploidy <- list(ascat = as.numeric(my.ascat.seg.ascn$ploidy), 
                median = median.ploidy, most.width = baseline.ploidy, 
                width.weighted = weighted.ploidy)
            # saveRDS(my.ascat.seg.ascn, paste0(samplename, ".ASCN.ASCAT.RDS"), 
            #     compress = "bzip2")

            # outfile <- paste0(samplename, ".gamma", gamma, 
            #     ".cn")
            # outdf <- my.ascat.seg.ascn$segments
            # outdf$chrom <- outdf$chr
            # outdf$chr <- unlist(cs$chrom2chr[outdf$chrom])
            # outdf$Width <- outdf$endpos - outdf$startpos + 1
            # outdf$TCN <- outdf$nMajor + outdf$nMinor
            # outdf <- outdf[, c(1, 2, 7, 3, 4, 8, 9, 5, 6)]
            # colnames(outdf)[1:5] <- c(samplename, "Chr", 
            #     "Chrom", "Start", "End")
            # write.table.fast(x = outdf, file = outfile)

            # outfile <- paste0(samplename, ".gamma", gamma, 
            #     "_model.txt")
            # modeldf <- data.frame(key = c("Sample", "Gamma", 
            #     "Goodness.of.Fit", "Psi", "Ploidy.ASCAT", 
            #     "Ploidy.Median", "Ploidy.Most.Width", 
            #     "Ploidy.Width.weighted", "Cellularity"), 
            #     value = c(samplename, gamma, unname(my.ascat.seg.ascn$goodnessOfFit), 
            #       unname(my.ascat.seg.ascn$psi), my.ascat.seg.ascn$ploidy$ascat, 
            #       my.ascat.seg.ascn$ploidy$median, my.ascat.seg.ascn$ploidy$most.width, 
            #       my.ascat.seg.ascn$ploidy$width.weighted, unname(my.ascat.seg.ascn$aberrantcellfraction)), 
            #     stringsAsFactors = FALSE)
            # write.table.fast(x = modeldf, file = outfile, header = FALSE)


            ylim <- 6
            segments.genostart <- cs$chromosomes$chr.length.toadd[outdf$Chr] + 
                my.ascat.seg.ascn$segments$startpos
            segments.genoend <- cs$chromosomes$chr.length.toadd[outdf$Chr] + 
                my.ascat.seg.ascn$segments$endpos
            ink <- cs$chromosomes$chrN %in% outdf$Chr
            # png(paste0(samplename, ".ASCN.ASCAT.png"), 
            #     width = 1850, height = 980)
            # par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
            # plot(0, 0, type = "n", xlim = c(0, cs$genome.length), 
            #     xaxs = "i", ylim = c(0, (ylim + 0.1)), 
            #     main = paste0(samplename, " Allele-Specific Copy Number (Gamma = ", 
            #       gamma, ")\n Cellularity = ", my.ascat.seg.ascn$aberrantcellfraction, 
            #       " // Ploidy = (A:", round(my.ascat.seg.ascn$ploidy$ascat, 
            #         digits = 2), ", M:", round(my.ascat.seg.ascn$ploidy$median, 
            #         digits = 2), ", MW:", round(my.ascat.seg.ascn$ploidy$most.width, 
            #         digits = 2), ", WW:", round(my.ascat.seg.ascn$ploidy$width.weighted, 
            #         digits = 2), ") // Goodness of fit = ", 
            #       round(my.ascat.seg.ascn$goodnessOfFit, digits = 2), 
            #       "% // Psi = ", my.ascat.seg.ascn$psi, 
            #       " // nSeg = ", nrow(my.ascat.seg.ascn$segments), 
            #       " // Predicted gender = ", data$data$gender), 
            #     xlab = "Genomic position", ylab = "ASCN", 
            #     xaxt = "n", cex.main = 2)
            # abline(v = c(segments.genostart, segments.genoend), 
            #     col = "grey90", lwd = 1)
            # segments(segments.genostart, my.ascat.seg.ascn$segments$nMajor + 
            #     0.05, segments.genoend, my.ascat.seg.ascn$segments$nMajor + 
            #     0.05, pch = ".", col = 2, lwd = 6)
            # segments(segments.genostart, my.ascat.seg.ascn$segments$nMinor - 
            #     0.05, segments.genoend, my.ascat.seg.ascn$segments$nMinor - 
            #     0.05, pch = ".", col = 3, lwd = 6)
            # up.nMajor <- which(my.ascat.seg.ascn$segments$nMajor > 
            #     ylim)
            # up.nMinor <- which(my.ascat.seg.ascn$segments$nMinor > 
            #     ylim)
            # dn.nMajor <- which(my.ascat.seg.ascn$segments$nMajor == 
            #     0)
            # dn.nMinor <- which(my.ascat.seg.ascn$segments$nMinor == 
            #     0)
            # if (length(up.nMajor) > 0) 
            #     segments(segments.genostart[up.nMajor], (ylim + 
            #       0.2) + 0.05, segments.genoend[up.nMajor], (ylim + 
            #       0.2) + 0.05, pch = ".", col = "orange", 
            #       lwd = 6)
            # if (length(up.nMinor) > 0) 
            #     segments(segments.genostart[up.nMinor], (ylim + 
            #       0.2) - 0.05, segments.genoend[up.nMinor], (ylim + 
            #       0.2) - 0.05, pch = ".", col = "lightgreen", 
            #       lwd = 6)
            # if (length(dn.nMajor) > 0) 
            #     segments(segments.genostart[dn.nMajor], 0.05, 
            #       segments.genoend[dn.nMajor], 0.05, pch = ".", 
            #       col = "darkred", lwd = 6)
            # if (length(dn.nMinor) > 0) 
            #     segments(segments.genostart[dn.nMinor], -0.05, 
            #       segments.genoend[dn.nMinor], -0.05, pch = ".", 
            #       col = "darkgreen", lwd = 6)
            # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, 
            #     lty = 3, lwd = 2)
            # abline(h = 0:ylim, col = "grey75", lty = 3)
            # try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, 
            #     labels = cs$chromosomes$chrom[ink], pos = 1, 
            #     cex = 1))
            # graphics::plot(0, 0, type = "n", xlim = c(0, 
            #     cs$genome.length), xaxs = "i", ylim = c(0, 
            #     (ylim + 0.1)), main = paste0(samplename, " Total Copy Number"), 
            #     xlab = "Genomic position", ylab = "TCN", 
            #     xaxt = "n", cex.main = 2)
            # abline(v = c(segments.genostart, segments.genoend), 
            #     col = "grey90", lwd = 1)
            # abline(h = my.ascat.seg.ascn$ploidy$median, col = "red", 
            #     lty = 2)
            # abline(h = my.ascat.seg.ascn$ploidy$most.width, col = "cyan", 
            #     lty = 2)
            # nTotal <- my.ascat.seg.ascn$segments$nMajor + my.ascat.seg.ascn$segments$nMinor
            # up.nTotal <- which(nTotal > ylim)
            # dn.nTotal <- which(nTotal == 0)
            # segments(segments.genostart, nTotal, segments.genoend, 
            #     nTotal, pch = ".", col = 4, lwd = 6)
            # if (length(up.nTotal) > 0) 
            #     segments(segments.genostart[up.nTotal], (ylim + 
            #       0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 
            #       0.2) + 0.05, pch = ".", col = "cyan", 
            #       lwd = 6)
            # if (length(dn.nTotal) > 0) 
            #     segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 
            #       0, pch = ".", col = "midnightblue", 
            #       lwd = 6)
            # abline(v = cs$chromosomes$chr.length.sum[ink], col = 1, 
            #     lty = 3, lwd = 2)
            # abline(h = 0:ylim, col = "grey75", lty = 3)
            # try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, 
            #     labels = cs$chromosomes$chrom[ink], pos = 1, 
            #     cex = 1))
            # dev.off()

            
            segments.posN <- unlist(cs$chrom2chr[my.ascat.seg.ascn$segments_raw$chr])
            segments.genostart <- cs$chromosomes$chr.length.toadd[segments.posN] + 
                my.ascat.seg.ascn$segments_raw$startpos
            segments.genoend <- cs$chromosomes$chr.length.toadd[segments.posN] + 
                my.ascat.seg.ascn$segments_raw$endpos
            # png(paste0(samplename, ".rawprofile.png"), 
            #     width = 1850, height = 980)
            # par(mar = c(1, 4, 5, 1), mfrow = c(2, 1))
            # plot(0, 0, type = "n", xlim = c(0, cs$genome.length), 
            #     xaxs = "i", ylim = c(0, (ylim + 0.1)), 
            #     main = paste0(samplename, " RAW Allele-Specific Copy Number (Gamma = ", 
            #       gamma, ")"), xlab = "Genomic position", 
            #     ylab = "ASCN", xaxt = "n", cex.main = 2)
            # abline(v = c(segments.genostart, segments.genoend), 
            #     col = "grey90", lwd = 1)
            # segments(segments.genostart, my.ascat.seg.ascn$segments_raw$nAraw + 
            #     0.05, segments.genoend, my.ascat.seg.ascn$segments_raw$nAraw + 
            #     0.05, pch = ".", col = 2, lwd = 6)
            # segments(segments.genostart, my.ascat.seg.ascn$segments_raw$nBraw - 
            #     0.05, segments.genoend, my.ascat.seg.ascn$segments_raw$nBraw - 
            #     0.05, pch = ".", col = 3, lwd = 6)
            # up.nMajor <- which(my.ascat.seg.ascn$segments_raw$nAraw > 
            #     ylim)
            # up.nMinor <- which(my.ascat.seg.ascn$segments_raw$nBraw > 
            #     ylim)
            # dn.nMajor <- which(my.ascat.seg.ascn$segments_raw$nAraw == 
            #     0)
            # dn.nMinor <- which(my.ascat.seg.ascn$segments_raw$nBraw == 
            #     0)
            # if (length(up.nMajor) > 0) 
            #     segments(segments.genostart[up.nMajor], (ylim + 
            #       0.2) + 0.05, segments.genoend[up.nMajor], (ylim + 
            #       0.2) + 0.05, pch = ".", col = "orange", 
            #       lwd = 6)
            # if (length(up.nMinor) > 0) 
            #     segments(segments.genostart[up.nMinor], (ylim + 
            #       0.2) - 0.05, segments.genoend[up.nMinor], (ylim + 
            #       0.2) - 0.05, pch = ".", col = "lightgreen", 
            #       lwd = 6)
            # if (length(dn.nMajor) > 0) 
            #     segments(segments.genostart[dn.nMajor], 0.05, 
            #       segments.genoend[dn.nMajor], 0.05, pch = ".", 
            #       col = "darkred", lwd = 6)
            # if (length(dn.nMinor) > 0) 
            #     segments(segments.genostart[dn.nMinor], -0.05, 
            #       segments.genoend[dn.nMinor], -0.05, pch = ".", 
            #       col = "darkgreen", lwd = 6)
            # abline(v = cs$chromosomes$chr.length.sum, col = 1, 
            #     lty = 3, lwd = 2)
            # abline(h = 0:ylim, col = "grey75", lty = 3)
            # try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, 
            #     labels = cs$chromosomes$chrom[ink], pos = 1, 
            #     cex = 1))
            # graphics::plot(0, 0, type = "n", xlim = c(0, 
            #     cs$genome.length), xaxs = "i", ylim = c(0, 
            #     (ylim + 0.1)), main = paste0(samplename, " RAW Total Copy Number"), 
            #     xlab = "Genomic position", ylab = "TCN", 
            #     xaxt = "n", cex.main = 2)
            # abline(v = c(segments.genostart, segments.genoend), 
            #     col = "grey90", lwd = 1)
            # abline(h = my.ascat.seg.ascn$ploidy$median, col = 2, 
            #     lty = 2)
            # abline(h = my.ascat.seg.ascn$ploidy$most.width, col = "purple", 
            #     lty = 2)
            # nTotal <- my.ascat.seg.ascn$segments_raw$nAraw + 
            #     my.ascat.seg.ascn$segments_raw$nBraw
            # up.nTotal <- which(nTotal > ylim)
            # dn.nTotal <- which(nTotal == 0)
            # segments(segments.genostart, nTotal, segments.genoend, 
            #     nTotal, pch = ".", col = 4, lwd = 6)
            # if (length(up.nTotal) > 0) 
            #     segments(segments.genostart[up.nTotal], (ylim + 
            #       0.2) + 0.05, segments.genoend[up.nTotal], (ylim + 
            #       0.2) + 0.05, pch = ".", col = "cyan", 
            #       lwd = 6)
            # if (length(dn.nTotal) > 0) 
            #     segments(segments.genostart[dn.nTotal], 0, segments.genoend[dn.nTotal], 
            #       0, pch = ".", col = "midnightblue", 
            #       lwd = 6)
            # abline(v = cs$chromosomes$chr.length.sum, col = 1, 
            #     lty = 3, lwd = 2)
            # abline(h = 0:ylim, col = "grey75", lty = 3)
            # try(text(x = cs$chromosomes$mid.chr.geno[ink], y = ylim, 
            #     labels = cs$chromosomes$chrom[ink], pos = 1, 
            #     cex = 1))
            # dev.off()


            cnpTotal <- my.ascat.seg.ascn$nA + my.ascat.seg.ascn$nB
            my.xlim <- c(min(cnpTotal, na.rm = TRUE) - 0.5, max(cnpTotal, 
                na.rm = TRUE) + 0.5)
            # png(paste0(samplename, ".TCNvsL2R.png"), width = 980, 
            #     height = 980)
            # par(mar = c(4, 4, 5, 1))
            # graphics::plot(cnpTotal, data$data$Tumor_LogR_segmented, 
            #     main = paste0(samplename, "\nCoherence of estimated total copy number (TCN) versus post-ASPCF segmented log2(ratio) (L2R)"), 
            #     xlab = "TCN", ylab = "L2R", xaxs = "i", 
            #     xlim = my.xlim, ylim = c(-1.5, 1.5), type = "n")
            # for (k in sort(unique(cnpTotal))) {
            #     my.yval <- data$data$Tumor_LogR_segmented[cnpTotal == 
            #       k]
            #     my.min <- min(my.yval, na.rm = TRUE)
            #     my.max <- max(my.yval, na.rm = TRUE)
            #     rect(xleft = my.xlim[1], ybottom = my.min, xright = my.xlim[2], 
            #       ytop = my.max, border = adjustcolor(k + 1, 
            #         alpha.f = 0.25), col = adjustcolor(k + 1, 
            #         alpha.f = 0.25))
            #     segments(k, my.min, k, my.max, col = k + 1)
            # }
            # points(cnpTotal, data$data$Tumor_LogR_segmented, 
            #     pch = ".", cex = 5, col = cnpTotal + 1)
            # abline(h = 0, lty = 2)
            # dev.off()


            # png(paste0(samplename, ".Rorschach.clown.png"), 
            #     width = 980, height = 980)
            # EaCoN.Rorschard.plot(data = data, cnpTotal = cnpTotal)
            # dev.off()


            tmsg(paste0("    ", round(my.ascat.seg.ascn$goodnessOfFit, 
                digits = 3), " / ", my.ascat.seg.ascn$psi))
            setwd(oridirx)
            return(unname(c(gamma, unlist(my.ascat.seg.ascn$ploidy, 
                use.names = FALSE), my.ascat.seg.ascn$aberrantcellfraction, 
                my.ascat.seg.ascn$goodnessOfFit, my.ascat.seg.ascn$psi)))
        }
    }, stringsAsFactors = FALSE)
    parallel::stopCluster(cls)
    rownames(fit.val) <- gammavec
    colnames(fit.val) <- c("gamma", "ploidy.ASCAT", 
        "ploidy.Median", "ploidy.Most.width", "ploidy.Width.weighted", 
        "aberrant.cell.fraction", "GoF", "psi")
    if (any(!is.na(fit.val$gamma))) {
        fit.val[, 1] <- gammavec
        gammaOpt.idx <- which.max(fit.val$GoF)
        gammaOpt <- fit.val$gamma[gammaOpt.idx]
        write.table(fit.val, file = paste0(odir, "/", samplename, 
            ".gammaEval.txt"), sep = "\t", quote = FALSE, 
            row.names = FALSE)
        # png(paste0(odir, "/", samplename, ".gammaEval.png"), 
        #     width = 1850, height = 980)
        # par(mfrow = c(3, 1), mar = c(2, 4, 3, 1), cex = 1)
        # graphics::plot(fit.val$gamma, fit.val$GoF, xlab = "Gamma", 
        #     ylab = "Goodness of fit", main = "Goodness of fit curve", 
        #     type = "b", pch = 20)
        # points(fit.val$gamma[gammaOpt.idx], fit.val$GoF[gammaOpt.idx], 
        #     pch = 20, col = 2)
        # abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
        # abline(h = fit.val$GoF[gammaOpt.idx], lty = 2, col = 2)
        # graphics::plot(fit.val$gamma, fit.val$psi, xlab = "Gamma", 
        #     ylab = "Psi", main = "Psi curve", type = "b", 
        #     pch = 20)
        # points(fit.val$gamma[gammaOpt.idx], fit.val$psi[gammaOpt.idx], 
        #     pch = 20, col = 2)
        # abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
        # ploidy.mat <- as.matrix(fit.val[, 2:5])
        # plo.ymax <- max(ploidy.mat, na.rm = TRUE)
        # graphics::plot(fit.val$gamma, fit.val$ploidy.ASCAT, xlab = "Gamma", 
        #     ylab = "Ploidy", main = "Ploidy : ASCAT (A=black), Median (M=red), Most width (MW=cyan), Width-weighted (WW=yellow)", 
        #     type = "b", pch = 20, col = "black", 
        #     ylim = c(0, plo.ymax), lwd = 2)
        # abline(h = seq.int(2, plo.ymax, 2), lty = 2, col = "grey50")
        # lines(fit.val$gamma, fit.val$ploidy.Median, type = "b", 
        #     pch = 20, col = "red", lwd = 2)
        # lines(fit.val$gamma, fit.val$ploidy.Most.width, type = "b", 
        #     pch = 20, col = "cyan", lwd = 2)
        # lines(fit.val$gamma, fit.val$ploidy.Width.weighted, type = "b", 
        #     pch = 20, col = "yellow", lwd = 2)
        # abline(v = fit.val$gamma[gammaOpt.idx], lty = 2, col = 2)
        # dev.off()
    }
    else {
        tmsg("WARNING : ASCN failed for all evaluated gamma values !")
    }
    return(my.ascat.seg.ascn)
}