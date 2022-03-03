library(oncoscanR)
setwd( "C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/results")

res = oncoscanR::workflow_oncoscan.run("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/2-AD/2-ADREC.RC.OSCHP.segments_FULL_LOCATION.txt", "F")
# res = oncoscanR::workflow_oncoscan.run("C:/Users/e.bordron/Desktop/CGH-scoring/M2_internship_Bergonie/data/working_data/chas_example.txt", "M")

armlevelCN = res$armlevel
losses = armlevelCN$LOSS
gains = armlevelCN$GAIN

nb_alter = length(losses) + length(gains)
nb_chr = 7 #calculated visually 

calcGI = function(nb_alter, nb_chr) {
  GI=(nb_alter^2)/nb_chr
  print(c("GI: ", GI))
  return(GI)
}

calcGI(nb_alter, nb_chr)

