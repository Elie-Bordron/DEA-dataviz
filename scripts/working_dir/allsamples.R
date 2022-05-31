## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
resDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
res_GI_dir = file.path(resDir, "GI_all_methods")
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(dplyr)
library("GGally")


# uniquePath = "C:/Users/warew/Desktop/CGH_scoring"


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
                #### create df with one col: result of first package
                vecGI = c(GI_table[,i])
                #### add col: group
                vecGrps = rep(colnames(GI_table[i]), dim(GI_table)[1])
                #### add col: color
                if(!is.null(agilentClass)) {
                    vecCol = as.vector(GI_table$AgilentClass)
                } else {
                    vecCol = rep(NULL, dim(GI_table)[1])
                }
            } else {
                vecGI = c(vecGI, GI_table[,i])
                vecGrps = c(vecGrps, rep(colnames(GI_table[i]), dim(GI_table)[1]))
                if(!is.null(agilentClass)) {
                    vecCol = c(vecCol, as.vector(GI_table$AgilentClass))
                } else {
                    vecCol = c(vecCol, rep(NULL, dim(GI_table)[1]))
                }
            }
        }
        if(all(is.null(vecCol))) {
            df2Cols = data.frame = data.frame(GI=vecGI, pkg=vecGrps)
        } else {
            df2Cols = data.frame = data.frame(GI=vecGI, pkg=vecGrps, color = vecCol)
        }
        return(df2Cols)  
    }

addAgilentClass = function(GI_table) {
        GI_table = GI_table %>%
            mutate(AgilentClass = case_when(GI_Agilent<10 ~ "low",
                # (GI_Agilent>=9.99) & (GI_Agilent<=10.01) ~ "intermediate",
                GI_Agilent>10 ~ "high"))
        return(GI_table)
    }

plotOnePkg = function(GI_table, pkg) {
    GI_col = paste0("GI_", pkg)
    dfOnePkg = dplyr::select(GI_table, c(sample, GI_col))
    colnames(dfOnePkg) = c("sample", "GI")
    dfOnePkg = dfOnePkg[order(dfOnePkg$GI), ]
    plot(dfOnePkg$GI, xaxt="n", main = paste0(pkg, " GI"))
    axis(1,at=1:length(dfOnePkg$GI),labels=dfOnePkg$sample, las = 2)
}

############################### main
if (interactive()) { 
    GI_table_path = file.path(res_GI_dir, "GI_all_methods.txt")
    GI_table = read.table(GI_table_path, h=T)
    ## R base correlation plots for all comparisons
    pairs(~ GI_oncoscanR + GI_CGHcall + GI_rCGH + GI_ASCAT + GI_Agilent, data = GI_table)

    ## correlation plots 
    GI_cols = dplyr::select(GI_table, c("GI_oncoscanR", "GI_CGHcall", "GI_rCGH", "GI_ASCAT", "GI_rCGH", "GI_Agilent"))
    ggpairs(GI_cols, lower=list(continuous=wrap("smooth", colour="blue")),
            diag=list(continuous=wrap("barDiag", fill="dark blue")),
            upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)

    ## remove aberrant values
    GI_table_edited = GI_table
    GI_table_edited = dplyr::mutate(GI_table_edited, GI_ASCAT = replace(GI_ASCAT, GI_ASCAT>200, NA))
    # GI_table_edited = dplyr::mutate(GI_table_edited, GI_rCGH = replace(GI_rCGH, GI_rCGH>80, NA))
    GI_table_edited = addAgilentClass(GI_table_edited)
    GI_cols_edited = dplyr::select(GI_table_edited, c("GI_oncoscanR", "GI_CGHcall", "GI_rCGH", "GI_ASCAT", "GI_rCGH", "GI_Agilent"))
    
    ## correlation plots but aberrant values removed
    ggpairs(GI_cols_edited, lower=list(continuous=wrap("smooth", colour="blue")),
            diag=list(continuous=wrap("barDiag", fill="dark blue")),
            upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)

    ## plot distribution  of values for each group
    GI2col = to2colDf(dplyr::select(GI_table_edited,  c("GI_oncoscanR", "GI_CGHcall", "GI_rCGH", "GI_ASCAT", "GI_rCGH", "GI_Agilent", "AgilentClass")), "AgilentClass")
    f <- ggplot(GI2col, aes(pkg, GI, fill = factor(color)))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=2, stroke=NA) +
        xlab("") + ylab("Genomic Index") +  labs(fill="")
    ## the same, zoomed in 
    f <- ggplot(GI2col, aes(pkg, GI, fill = factor(color)))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=0.5, stroke=NA) +
        xlab("") + ylab("Genomic Index") + labs(fill="") + ylim(0,30)
    
    
    
    
    ## plot pkgs individually
    row.names(GI_table) = GI_table$sample
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
    plot(y=purityTable$estimPurity_ASCAT, x=purityTable$purity_HES, xlab = "HES-estimated purity", ylab = "ASCAT-estimated purity", main = "comparison of purity estimates")

    
    
    ## plot runTime
    # GI_table = addAgilentClass(GI_table)
    runTime2col = to2colDf(dplyr::select(GI_table, contains("runTime")))
    colnames(runTime2col) = c("runTime", "pkg")
    f <- ggplot(runTime2col, aes(pkg, runTime))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=1) +
        xlab("") + ylab("Processing time (s)")
}



