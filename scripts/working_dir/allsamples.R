## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
resDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
res_GI_dir = file.path(resDir, "GI_all_methods")
setwd(working_dir)
## open working directory in Files tab
options("max.print"=100)
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(dplyr)
library("GGally")



to2colDf = function(GI_table, agilentClass=NULL) {
        ## provide only columns with GI values + possibly column agilentClass, which must be the last.
        vecGI=NULL
        vecGrps=NULL
        if(is.null(agilentClass)) {
            nbCols = dim(GI_table)[2]
        } else {
            nbCols = dim(GI_table)[2] - 1
        }
        for (i in 1:nbCols) {
            if(is.null(vecGI)) {
                #### first col of future df: result of first package
                vecGI = c(GI_table[,i])
                #### 2nd col: group
                vecGrps = rep(colnames(GI_table[i]), dim(GI_table)[1])
                vecGrpsNum = rep(i, dim(GI_table)[1])
                #### add col: color
                if(!is.null(agilentClass)) {
                    vecCol = as.vector(GI_table$AgilentClass)
                } else {
                    vecCol = rep(NULL, dim(GI_table)[1])
                }
            } else {
                vecGI = c(vecGI, GI_table[,i])
                vecGrps = c(vecGrps, rep(colnames(GI_table[i]), dim(GI_table)[1]))
                vecGrpsNum = c(vecGrpsNum, rep(i, dim(GI_table)[1]))
                if(!is.null(agilentClass)) {
                    vecCol = c(vecCol, as.vector(GI_table$AgilentClass))
                } else {
                    vecCol = c(vecCol, rep(NULL, dim(GI_table)[1]))
                }
            }
        }
        if(all(is.null(vecCol))) {
            df2Cols = data.frame(GI=vecGI, pkg=vecGrps, grpNum=vecGrpsNum)
        } else {
            df2Cols = data.frame(GI=vecGI, pkg=vecGrps, grpNum=vecGrpsNum, color = vecCol)
        }
        return(df2Cols)  
    }

addAgilentClass = function(GI_table) {
        GI_table = GI_table %>%
            mutate(AgilentClass = case_when(GI_Agilent<10 ~ "faible",
                # (GI_Agilent>=9.99) & (GI_Agilent<=10.01) ~ "intermediate",
                GI_Agilent>10 ~ "haut"))
        return(GI_table)
    }

plotOnePkg = function(GI_table, pkg) {
    GI_col = paste0("GI_", pkg)
    dfOnePkg = dplyr::select(GI_table, c(sample, GI_col))
    colnames(dfOnePkg) = c("sample", "GI")
    dfOnePkg = dfOnePkg[order(dfOnePkg$GI), ]
    plot(dfOnePkg$GI, xaxt="n", main = paste0(pkg, " GI"), xlab="", ylab="GI")
    axis(1,at=1:length(dfOnePkg$GI),labels=dfOnePkg$sample, las = 2)
    
}

############################### main
if (sys.nframe() == 0){
    GI_table_path = file.path(res_GI_dir, "GI_all_methods.txt")
    GI_table = read.table(GI_table_path, h=T)
    GI_table = dplyr::filter(GI_table, sample!="17-VV")


    ## correlation plots 
    GI_cols = dplyr::select(GI_table, c("GI_Agilent", "GI_oncoscanR", "GI_rCGH", "GI_CGHcall", "GI_ASCAT"))
    grpNames = colnames(GI_cols)
    grpNames = substr(grpNames, 4, nchar(grpNames))
    ggpairs(GI_cols, lower=list(continuous=wrap("smooth", colour="blue")),
            diag=list(continuous=wrap("barDiag", fill="dark blue")),
            upper=list(corSize=6), axisLabels='show', cardinality_threshold=20, columnLabels=grpNames)

    
    ## remove aberrant values
    GI_table_edited = GI_table
    GI_table_edited = dplyr::mutate(GI_table_edited, GI_ASCAT = replace(GI_ASCAT, GI_ASCAT>200, NA))
    # GI_table_edited = dplyr::mutate(GI_table_edited, GI_rCGH = replace(GI_rCGH, GI_rCGH>80, NA))
    GI_table_edited = addAgilentClass(GI_table_edited)
    GI_cols_edited = dplyr::select(GI_table_edited, c("GI_Agilent", "GI_oncoscanR", "GI_rCGH", "GI_CGHcall", "GI_ASCAT"))
    
    ## correlation plots but aberrant values removed
    # ggpairs(GI_cols_edited, lower=list(continuous=wrap("smooth", colour="blue")),
    #         diag=list(continuous=wrap("barDiag", fill="dark blue")),
    #         upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)

    
    ## plot distribution  of values for each group
    GI2col = to2colDf(dplyr::select(GI_table_edited,  c("GI_Agilent", "GI_oncoscanR", "GI_rCGH", "GI_CGHcall", "GI_ASCAT", "AgilentClass")), "AgilentClass")
    ## add grp for lines
    unit = c(GI_table_edited$sample)
    GI2col$lineGrp = rep(unit, 5)
    f <- ggplot(GI2col, aes(pkg, GI, fill = factor(color)))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=2, stroke=NA) + 
        # geom_line() +
        xlab("") + ylab("Genomic Index") +  labs(fill="groupe selon Agilent")
    ## the same, with lines joining same points
    library(ggrepel)
    grpNames = unique(GI2col$pkg)
    grpNames = substr(grpNames, 4, nchar(grpNames))
    gg = ggplot(GI2col)
    ## uncomment next line to draw lines between points of same sample
    # gg = gg + geom_line(aes(x=grpNum, y=GI, group=lineGrp, color=color))
    gg = gg + geom_point(data=GI2col, aes(y=GI,x=grpNum,color=color), size=4, alpha=0.5, position = position_jitter(0.1))
    gg = gg + scale_x_continuous(labels=grpNames)
    gg = gg + theme_bw() + xlab(NULL) + guides(color = guide_legend(title = "Groupe selon Agilent"))
    gg = gg + geom_abline(intercept = 10, slope=0, alpha=0.5)
    gg = gg + geom_abline(intercept = 21, slope=0, alpha=0.5)
    colors <- c("faible" = "#00BFC4", "haut" = "#F8766D")
    gg = gg + scale_color_manual(values = colors)
    # gg = gg + geom_label_repel(aes(y=GI,x=grpNum,label=lineGrp), box.padding = 0.5, point.padding = 0.1, segment.color = 'grey50')
    gg
    
    
    ### geom_label_repel
    if (FALSE) {
        nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv", sep = ",")
        nbaplot <- ggplot(nba, aes(x= MIN, y = PTS)) + geom_point(color = "blue", size = 3)
        nbaplot + geom_label_repel(aes(label = Name), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +theme_bw()
    }
    
    #### test for adding lines
    if (FALSE) {
        val=c(5,6,9,8,8,7,7,3,2,6,5,9)
        x = as.data.frame(val)
        x$ref = c(1,1,1,1,2,2,2,2,3,3,3,3)
        x$indivs = c("a","b","c","d","a","b","c","d","a","b","c","d")
        gg = ggplot(x)
        gg = gg + geom_point(data=x, aes(y=val,x=ref), size=5, alpha=0.3)
        gg = gg + geom_line(aes(x=ref, y=val, color=indivs))
        gg
    }
    
    
    ## plot pkgs individually
    row.names(GI_table) = GI_table$sample
    o = GI_table$GI_oncoscanR
    agilent = GI_table$GI_Agilent
    ggplot(data=GI_table, aes(x=GI_Agilent, y=GI_oncoscanR)) + theme_bw()
    plot(GI_table$GI_Agilent, GI_table$GI_oncoscanR, ylab="GI oncoscanR", xlab="GI Agilent")
    
    plotOnePkg(GI_table, "oncoscanR")
    plotOnePkg(GI_table, "CGHcall")
    plotOnePkg(GI_table, "rCGH")
    plotOnePkg(dplyr::filter(GI_table_edited, GI_rCGH <80), "rCGH")
    plotOnePkg(GI_table_edited, "ASCAT")
    plotOnePkg(GI_table, "Agilent")


    ## compare ASCAT-estimated purity vs HES-estimated purity
    purityTablePath = file.path(res_GI_dir, "cellularity.txt")
    purityTable = read.table(purityTablePath, h=T)
    purityTable = purityTable[order(purityTable$purity_HES), ]
    purityTable$estimPurity_ASCAT = as.numeric(purityTable$estimPurity_ASCAT)*100
    plot(y=purityTable$estimPurity_ASCAT, x=purityTable$purity_HES, xlab = "HES-estimated purity", ylab = "ASCAT-estimated purity", main = "comparison of purity estimates", xlim=c(0,100), ylim=c(0,100))

    
    
    ## plot runTime
    # GI_table = addAgilentClass(GI_table)
    runTime2col = to2colDf(dplyr::select(GI_table, contains("runTime")))
    colnames(runTime2col) = c("runTime", "pkg")
    grpNames = unique(runTime2col$pkg)
    grpNames = substr(grpNames, 9, nchar(grpNames))
    f <- ggplot(runTime2col, aes(pkg, runTime))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=1) +
        xlab("") + ylab("Temps de calcul") +
        theme_bw()
    

    ### plot ROC curves
    library(pROC)
    subset = dplyr::filter(GI2col, pkg=="GI_oncoscanR")
    resRoc = roc(subset$color, subset$GI)
    gg = ggroc(resRoc, colour = 'steelblue', size=2, legacy.axes=TRUE)
    gg + theme_bw() + ggtitle("oncoscanR")
    print(resRoc)
    
    subset = dplyr::filter(GI2col, pkg=="GI_rCGH")
    resRoc = roc(subset$color, subset$GI)
    gg = ggroc(resRoc, colour = 'steelblue', size=2, legacy.axes=TRUE)
    gg + theme_bw() + ggtitle("rCGH")
    print(resRoc)
    
    
    subset = dplyr::filter(GI2col, pkg=="GI_CGHcall")
    resRoc = roc(subset$color, subset$GI)
    gg = ggroc(resRoc, colour = 'steelblue', size=2, legacy.axes=TRUE)
    gg + theme_bw() + ggtitle("CGHcall")
    print(resRoc)
    
    subset = dplyr::filter(GI2col, pkg=="GI_ASCAT")
    resRoc = roc(subset$color, subset$GI)
    gg = ggroc(resRoc, colour = 'steelblue', size=2, legacy.axes=TRUE)
    gg + theme_bw() + ggtitle("ASCAT")
    print(resRoc)
    
    ## ROC curves
    data(aSAH)
    
    ## Basic
    roc(aSAH$outcome, aSAH$s100b)
    roc(outcome ~ s100b, aSAH)
    
    ## smoothing
    x = roc(outcome ~ s100b, aSAH, smooth=TRUE) 
    plot(x, legacy.axes=TRUE)
    
    
    ### compute Spearman cor coeff
    
    pkg = "oncoscanR"
    cor.test(x=GI_table[paste0("GI_", pkg)][[1]], y=GI_table$GI_Agilent, method="pearson")
    
    pkg = "rCGH"
    cor.test(x=GI_table[paste0("GI_", pkg)][[1]], y=GI_table$GI_Agilent, method="pearson")
    
    pkg = "CGHcall"
    cor.test(x=GI_table[paste0("GI_", pkg)][[1]], y=GI_table$GI_Agilent, method="pearson")
    
    pkg = "ASCAT"
    cor.test(x=GI_table[paste0("GI_", pkg)][[1]], y=GI_table$GI_Agilent, method="pearson")
    
    
    
    
    ###### plot WGV
    ## imports
    source(file.path(working_dir, "CGHcall_functions.R")) # for plotSegTable()
    source(file.path(working_dir ,"rCGH_functions.R")) # To use removePointsForQuickPlotting()
    
    ## colors
    customColors = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
    
    alreadyGoodPos=FALSE
    sample = "6-VJ"
    pathProbeData = "C:/Users/e.bordron/Desktop/CGH-scoring/data/working_data/from_laetitia/premiers_E_recus/all_probeset"
    ## 1: load probe data for current sample
    loadSegTable = function(sample, pkg, resDir) {
        ## load file
        pkgDir = file.path(resDir, pkg, "segTables")
        segtabPath = paste0(pkgDir, "/", sample, ".tsv")
        segTable = read.table(segtabPath, sep='\t', h=TRUE)
        ## adapt columns according to pkg
        if(pkg=="ASCAT") {
            # segTable$logR = segTable$nAraw + segTable$nBraw 
            # segTable = dplyr::select(segTable, -c("nMajor", "nMinor", "nAraw", "nBraw"))
            # colnames(segTable) = c("sample" ,"chrom", "loc.start", "loc.end", "CN")
            colnames(segTable) = c("chrom", "loc.start", "loc.end", "CN", "nbProbes")
            chrom_as_str = unlist(stringr::str_split(segTable$chrom, "chr"))
            chrom_as_str = chrom_as_str[nzchar(chrom_as_str)]
            segTable$chrom = chrom_as_str
            
        } else if(pkg=="rCGH") {
            colnames(segTable) = c("chrom", "loc.start", "loc.end", "CN", "nbProbes") 
            
        } else if(pkg=="CGHcall") {
            colnames(segTable) = c("chrom", "loc.start", "loc.end", "CN") 
            
        } else if(pkg=="OncoscanR") {
            colnames(segTable) = c("sample", "chrom", "loc.start", "loc.end", "CN")
        }
        return(segTable)
    }
    resDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
    pkgs = c("OncoscanR", "rCGH", "CGHcall", "ASCAT")
    
    ## 2: load raw probes data
    probeset_txt_file = paste0(pathProbeData, "/", sample, ".probeset.txt")
    rawPrbData = read.table(probeset_txt_file, h=TRUE, sep='\t')
    colnames(rawPrbData) = c("ProbeName","ChrNum","ChrStart","Log2Ratio","WeightedLog2Ratio","Allele.Difference","NormalDiploid","BAF")
    rawPrbData = dplyr::filter(rawPrbData, ChrNum<23)
    rawPrbData = removePointsForQuickPlotting(rawPrbData)
    rawPrbData = getNewPos(rawPrbData)
    rawPrbData = rawPrbData[c(1:4, length(rawPrbData))]
    ## 3: plot both

    for (pkg in pkgs) {
        print(c("===pkg===: ", pkg))
        segTable = loadSegTable(sample, pkg, resDir)
        if(pkg=="OncoscanR") {alreadyGoodPos=TRUE}
        plot(y=rawPrbData$Log2Ratio, x=rawPrbData$absPos, pch = 20, cex=0.01, col="dark grey", xlab = "", ylab = "log Ratio", main = paste0(pkg, "  ", sample), ylim = c(-2,2))
        print(c("segTable: ", segTable))
        plotSegTableForWGV(segTable, sample, savePlot=FALSE, genGrid=FALSE, segColor="dark red", alreadyGoodPos=alreadyGoodPos) # segtables must have these columns: chrom, loc.start, loc.end, CN
        alreadyGoodPos = FALSE
    }
    
}

layout_matrix <- matrix(c(1:5), ncol = 1)
layout(layout_matrix)
layout(1)
plot(1)