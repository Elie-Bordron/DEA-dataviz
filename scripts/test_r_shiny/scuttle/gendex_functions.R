buildFileName <- function(res_dir="", prefix="", extension=".csv") {
    if(prefix!="")prefix = paste0(prefix, "_")
    paste0(res_dir, "/", prefix, "out_", Sys.Date(), extension)
}