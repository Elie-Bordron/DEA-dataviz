saveGI_ResToFile = function(currPkgDfGI,pkgName,addColsToSave=NULL,GI_dir = "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results/GI_all_methods") {
    ###### input dataframe currPkgDfGI's rows must exactly match allGiTable$sample column.
    ## load all GI results file
    GI_filePath =  file.path(GI_dir, "GI_all_methods.txt")
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
            # resList[currColName]<- currPkgDfGI[col]
            # print(c("resList: ", resList))
            # print(allGiTable["qkdkqdjh"])
            # allGiTable$paste0("nbAlter_",pkgName) = currPkgDfGI$nbAlter
            # allGiTable$paste0("nbChr_",pkgName) = currPkgDfGI$nbChr
            # allGiTable$paste0("runTime_",pkgName) = currPkgDfGI$runTime
        }
    } else {
        stop("row names don't match")
    }
    ## save the modified table in a file
    print(c("allGiTable: ", allGiTable))
    # print(c("allGiTable_copy: ", allGiTable_copy))
    write.table(allGiTable, GI_filePath, sep="\t", row.names=FALSE, quote=F)
}

