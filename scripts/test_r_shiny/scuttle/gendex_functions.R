buildFileName <- function(res_dir="", prefix="", extension=".csv") {
    if(prefix!="")prefix = paste0(prefix, "_")
    paste0(res_dir, "/", prefix, "out_", Sys.Date(), extension)
}


addSuffix = function(elmt, suffix="default_suffix") {
    return(paste0(elmt, suffix))
}

addSuffixToStrVec = function(strVec, suffix) {
    for (i in 1:length(strVec)){
        strVec[i] = addSuffix(strVec[i], suffix)
    }
    return(strVec)
}

