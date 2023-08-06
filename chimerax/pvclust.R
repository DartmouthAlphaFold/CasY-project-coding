#!/usr/bin/Rscript

# pvclust - UPGMA/WPGMA/Ward

cat("\014")
rm(list = ls())
try(dev.off(dev.list()["RStudioGD"]), silent=TRUE)

library(pvclust)
library(openxlsx)

input_path <- "C:/XResearch/Matchmaking/selected/alignment.xlsx"
sheets <- openxlsx::getSheetNames(input_path)
data <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=input_path, colNames=TRUE, rowNames=TRUE)
names(data) <- sheets

output_path <- "C:/XResearch/Matchmaking/selected/"
setwd(output_path)

set.seed(0)
pdf("pvclust.pdf")

startTime <- Sys.time()
i = 1
for (df in data){
  df_name <- names(data)[i]
  print(df_name)
  
  # convert to uncondensed matrix
  df[upper.tri(df)] <- t(df)[upper.tri(df)]
  
  # create the dendrogram
  pv_result <- pvclust(df, nboot = 10000,  parallel=TRUE, 
                       method.hclust = "ward.D2", method.dist = "correlation")
  
  # plot
  plot(pv_result, cex = 0.8, main=paste("Cluster Dendrogram\n", df_name))
  pvrect(pv_result, alpha = 0.95)
  
  # check statistics by node
  x <- seplot(pv_result, identify=FALSE, main = paste("p-value vs. SE plot\n", df_name))
  print(pv_result, which = x, digits=5)
  
  # pick significant clusters
  clusters <- pvpick(pv_result, alpha = 0.95)
  print(clusters)
  
  i <- i+1
}  

dev.off()
endTime <- Sys.time()
print(endTime - startTime)


# test graph
rmsd <- data$`RMSD (All pairs)`
rmsd[upper.tri(rmsd)] <- t(rmsd)[upper.tri(rmsd)]
res <- pvclust(rmsd, nboot = 10000,  parallel=TRUE, 
                     method.hclust = "ward.D2", method.dist = "correlation")
pdf("rmsd.pdf", width=13, height=10)
plot(res, cex=0.8, cex.pv=0.6, 
     print.num=FALSE, print.pv=TRUE,
     main="Cluster Dendrogram of RMSD")
pvrect(res, alpha = 0.95)
dev.off()
