saveGI_ResToFile = function(currPkgDfGI,pkgName,addColsToSave=NULL,GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods") {
    ###### input dataframe currPkgDfGI's rows must exactly match allGiTable$sample column.
    ## load all GI results file
    GI_filePath =  file.path(GI_dir, "GI_all_methods.txt")
    # print(c("GI_filePath: ", GI_filePath))
    allGiTable = read.table(GI_filePath, h=T)
    ## write results of current package in it
    if(!is.null(addColsToSave)) {
        colsToSave=c("GI", "nbAlter", "nbChr", "runTime", addColsToSave)
    } else {
        colsToSave=c("GI", "nbAlter", "nbChr", "runTime")
    }
    # print(paste0("saving columns ", colsToSave))
    if(all(row.names(currPkgDfGI)==allGiTable$sample)) {
        # allGiTable_copy = allGiTable
        resList = list()
        for(col in colsToSave) {
            currColName = paste0(col,"_",pkgName)
            print(currColName)
            allGiTable[currColName] <- currPkgDfGI[col]
        }
    } else {
        stop("row names don't match")
    }
    ## save the modified table in a file
    write.table(allGiTable, GI_filePath, sep="\t", row.names=FALSE, quote=F)
}

