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
    ## process CEL files on linux since apt.oncoscan.process() is 16-bit on windows
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
    pbn_217251 = as.vector(my.oschp$Genotyping$Calls$ProbeSetName)
    CNdf = my.oschp$ProbeSets$CopyNumber
    pbn_217611 = as.vector(CNdf$ProbeSetName)
    missingProbes = setdiff(pbn_217611, pbn_217251)
    
    # to avoid the dataframes from having different numbers of rows. some rows (grossly 300 of them) are removed because they have NULL values of log ratio.
    my.oschp$AlgorithmData$MarkerABSignal = filter(my.oschp$AlgorithmData$MarkerABSignal, !(ProbeSetName %in% missingProbes))
    my.oschp$ProbeSets$AllelicData = filter(my.oschp$ProbeSets$AllelicData, !(ProbeSetName %in% missingProbes))
    my.oschp$ProbeSets$CopyNumber = filter(my.oschp$ProbeSets$CopyNumber, !(ProbeSetName %in% missingProbes))
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
        # print(" destination of .RDS: ")
        # print(paste0(samplename, "_", arraytype, "_", genome, "_processed.RDS"))
        # saveRDS(my.ascat.obj, paste0(out.dir, "/", samplename, "/", samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "bzip2")
        saveRDS(my.ascat.obj, paste0(samplename, "_", arraytype, "_", genome, "_processed.RDS"), compress = "bzip2")
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
                                                        l2r.level, ", median-centered)) / ", round(sum(abs(diff(l2r.rm))), 
                                                                                                   digits = 3)), xlab = "Genomic position", 
             ylab = "L2R")
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



















