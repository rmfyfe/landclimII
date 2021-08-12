#######################################################################################
## LOOPING REVEALS                                                                   ##
## REVEALS-in-R written by Petr Kunes available at https://github.com/petrkunes/LRA  ##
## "wrapper" script: Ralph Fyfe                                                      ##    ##
#######################################################################################

#set working directory 
#wd must include r script and all data and list files
setwd("C:/temp/github_training/landclimII/run.REVEALS.files")

#load R function
source("LRA-REVEALS-release.R")

### SET INPUT VARIABLES
#PPE files for Europe
avg <- read.table("av-RPP.all.landclimII.csv", row.names = 1, header = T, sep = ",")
alvc <- read.table("avar-RPP.all.landclimII.csv", sep = ",", row.names = 1, header = T,check.names=F)
#list of grid cell files
gc.file.list <- read.csv("GC.list.csv", header=F)
#taxa names for writing out results
headings <- read.csv("taxa.csv", header=F)

### ### ### ### ### ### ###
### DO NOT EDIT BELOW   ###
### ### ### ### ### ### ###


#set up output files (for 25 time windows used in LANDCLIMII)

for (j in 1:25){
  outfile1 <- paste0("TW.", j, ".mean.csv")
  write.table(headings, file=outfile1, append=T, quote=F, row.names=F, col.names=F, sep=",")
  outfile2 <- paste0("TW.", j, ".SE.csv")
  write.table(headings, file=outfile2, append=T, quote=F, row.names=F, col.names=F, sep=",")
}

#run loop for performing REVEALS
count <- 1
for (i in 1:nrow(gc.file.list)){
  grid <- as.character(gc.file.list[i,])
  print(paste("grid cell", count, "of", nrow(gc.file.list), ":", grid))
  sitelist <- read.csv(grid, header=F)
  
  #generate list of sites and data
  largesites <- list()
  for(i in 1:nrow(sitelist)){
    name <- sitelist[i,]
    sitedata <- read.table(as.character(sitelist[i,]), sep=",", row.names=1, check.names=F, header=T)
    largesites[[name]] <- sitedata
  }
  #run RV
  estimates <- REVEALS(largesites, avg, alvc, 3, 100000, dwm = "gpm neutral")
  
  #write data for grid cell to file
  outfileMean <- paste0("output.mean.", grid)
  outfileSE <- paste0("output.SE.", grid)
  write.csv(100*estimates$Mean_V, outfileMean)
  write.csv(100*estimates$Mean_SE, outfileSE)

  #write data for time windows to file
  for (j in 1:25){
    outfile1 <- paste0("TW.", j, ".mean.csv")
    write.table(t(c(grid, 100*estimates$Mean_V[,j])), file=outfile1, append=T, quote=F, row.names=F, col.names=F, sep=",")
    outfile2 <- paste0("TW.", j, ".SE.csv")
    write.table(t(c(grid, 100*estimates$Mean_SE[,j])), file=outfile2, append=T, quote=F, row.names=F, col.names=F, sep=",")
  }
  count <- count + 1
}

