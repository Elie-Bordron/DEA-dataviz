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
if (sys.nframe() == 0){
    GI_table_path = file.path(res_GI_dir, "GI_all_methods.txt")
    GI_table = read.table(GI_table_path, h=T)
    GI_table = dplyr::filter(GI_table, sample!="17-VV")
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
    # ggpairs(GI_cols_edited, lower=list(continuous=wrap("smooth", colour="blue")),
    #         diag=list(continuous=wrap("barDiag", fill="dark blue")),
    #         upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)

    ## plot distribution  of values for each group
    GI2col = to2colDf(dplyr::select(GI_table_edited,  c("GI_oncoscanR", "GI_CGHcall", "GI_rCGH", "GI_ASCAT", "GI_rCGH", "GI_Agilent", "AgilentClass")), "AgilentClass")
    ## add grp for lines
    unit  = c(GI_table_edited$sample)
    GI2col$lineGrp = rep(unit, 5)
    f <- ggplot(GI2col, aes(pkg, GI, fill = factor(color)))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=2, stroke=NA) + 
        # geom_line() +
        xlab("") + ylab("Genomic Index") +  labs(fill="")
    ## the same, with lines joining same points
    library(ggrepel)
    grpNames = unique(GI2col$pkg)
    
    gg = ggplot(GI2col)
    ## uncomment next line to draw lines between from same sample
    # gg = gg + geom_line(aes(x=grpNum, y=GI, group=lineGrp, color=color))
    gg = gg + geom_point(data=GI2col, aes(y=GI,x=grpNum,color=color), size=5, alpha=0.5, position = position_jitter(0.1))
    gg = gg + scale_x_continuous(labels=grpNames)
    gg = gg + theme_bw() + xlab(NULL)
    gg = gg + geom_label_repel(aes(y=GI,x=grpNum,label=lineGrp), box.padding = 0.2, point.padding = 0.1, segment.color = 'grey50')
    gg
    
    library(ggplot2)
    
    nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv", sep = ",")
    
    nbaplot <- ggplot(nba, aes(x= MIN, y = PTS)) + 
        geom_point(color = "blue", size = 3)
    
    ### geom_label_repel
    nbaplot + 
        geom_label_repel(aes(label = Name),
                         box.padding   = 0.35, 
                         point.padding = 0.5,
                         segment.color = 'grey50') +
        theme_classic()
    
    #### test 
    val=c(5,6,9,8,8,7,7,3,2,6,5,9)
    x = as.data.frame(val)
    x$ref = c(1,1,1,1,2,2,2,2,3,3,3,3)
    x$indivs = c("a","b","c","d","a","b","c","d","a","b","c","d")
    gg = ggplot(x)
    gg = gg + geom_point(data=x, aes(y=val,x=ref), size=5, alpha=0.3)
    gg = gg + geom_line(aes(x=ref, y=val, color=indivs))
    gg
    
    
    # Change Color of a R ggplot Jitter
    
    # Importing the ggplot2 library
    library(ggplot2)
    
    # Creating basic Jitter
    ggplot(ChickWeight, aes(x = Diet, y = weight)) + 
        geom_jitter(aes(colour = Diet))

        
    plot(ChickWeight$Diet, ChickWeight$weight)    
    
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
    purityTable$estimPurity_ASCAT = as.numeric(purityTable$estimPurity_ASCAT)*100
    plot(y=purityTable$estimPurity_ASCAT, x=purityTable$purity_HES, xlab = "HES-estimated purity", ylab = "ASCAT-estimated purity", main = "comparison of purity estimates", xlim=c(0,100), ylim=c(0,100))

    
    
    ## plot runTime
    # GI_table = addAgilentClass(GI_table)
    runTime2col = to2colDf(dplyr::select(GI_table, contains("runTime")))
    colnames(runTime2col) = c("runTime", "pkg")
    f <- ggplot(runTime2col, aes(pkg, runTime))
    f + geom_dotplot(binaxis = "y", stackdir = "centerwhole", binwidth=1) +
        xlab("") + ylab("Processing time (s)")
}



