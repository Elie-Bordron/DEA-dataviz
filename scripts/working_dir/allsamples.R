## set working directory
working_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/scripts/working_dir"
resDir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results"
setwd(working_dir)
## open working directory in Files tab
rstudioapi::filesPaneNavigate(working_dir)
#loading libraries
library(dplyr)
library("GGally")





if (!interactive()) { ## main
    GI_table_path = file.path(resDir, "GI_all_methods/GI_all_methods.txt")
    GI_table = read.table(GI_table_path, h=T)
    pairs(~ GI_oncoscanR + GI_CGHcall + GI_rCGH + GI_ASCAT + GI_Agilent, data = GI_table)
    GI_sorted = sort(GI_table$GI_CGHcall)
    plot(GI_sorted)
    
    ## correlation plots 
    ggpairs(GI_table[, 2:6], lower=list(continuous=wrap("smooth", colour="blue")),
        diag=list(continuous=wrap("barDiag", fill="dark blue")),
        upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)
    
    ## correlation plots but aberrant values removed
    GI_table_edited = GI_table
    GI_table_edited = dplyr::mutate(GI_table_edited, GI_ASCAT = replace(GI_ASCAT, GI_ASCAT>400, NA))
    ggpairs(GI_table_edited[, 2:6], lower=list(continuous=wrap("smooth", colour="blue")),
        diag=list(continuous=wrap("barDiag", fill="dark blue")),
        upper=list(corSize=6), axisLabels='show', cardinality_threshold=20)
    
    plot(GI_table_edited$GI_oncoscanR, GI_table_edited$GI_CGHcall)
    # 
    # 
    # ggpairs(iris[, 1:4], lower=list(continuous="smooth", params=c(colour="blue")),
    #         diag=list(continuous="bar", params=c(colour="blue")), 
    #         upper=list(params=list(corSize=6)), axisLabels='show')
    # 
}













