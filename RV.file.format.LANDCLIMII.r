###################################################################
## aggregate data and prepare files for LRA (REVEALS) LANDCLIMII ##
###################################################################

#requires:
#
# (1) correct directory structure
#     dir                  [root for all lookup files]
#     dir\\input           [all unformatted count data]
# (2) lookup files (in dir)
#     LANDCLIM.ALLEUROPE.taxa.LUT.080219.csv [2 cols with old and new var names]
#     run.all.taxa.100219.csv                [single column of all PPE taxa to use, INCLUDES age.cal.BP]
#     bins.csv                               [single column of break points to define time windows]
#     TWlookup.csv                           [2 cols: col1 numbers 1:n; col2 mid-point of TW]
#     site.metadata.csv                      [minimum 3 cols: col1="filenamecsv" col2="radius" col3="model"]
# (3) count data (in dir\\input)
#     rows as samples, cols as taxa
#     formatted with col1 age.cal.BP col2:n un-harmonised taxa

##############################################################
## step 1: set working directories and load essential files ##
##############################################################

#define the working directory here
wd <- "C:\\temp\\github_training\\landclimII"

#check/build directory structure
input1 <- paste0(wd, "\\input\\")

outfolder1 <- paste0(wd, "\\temp.output01\\")
if(dir.exists(outfolder1) == TRUE){
  print("outfolder1 directory exists")
} else {
  dir.create(outfolder1)
  print("outfolder1 directory created")
}

outfolder2 <- paste0(wd, "\\temp.output02\\")
if(dir.exists(outfolder2) == TRUE){
  print("outfolder2 directory exists")
} else {
  dir.create(outfolder2)
  print("outfolder2 directory created")
}

outfolder3 <- paste0(wd, "\\LRA.ready.files\\")
if(dir.exists(outfolder3) == TRUE){
  print("LRA.ready.files directory exists")
} else {
  dir.create(outfolder3)
  print("LRA.ready.files directory created")
}

##read in necessary formatting files  
setwd(wd)

#read taxon harmonisation look up table
lookup <- read.csv("LANDCLIM.ALLEUROPE.taxa.LUT.080219.csv", header=T, check.names=F)
#list of RPP taxa
taxon <- read.csv("run.all.taxa.100219.csv", header=F)
#break points between time windows
bins <- as.numeric(unlist(read.csv("bins.csv", header=F)))
#lookup table for time windows
TWlookup <- read.csv("TWlookup.csv", header=T)
#site metadata including site size (radius) and site type (1=bog 2=lake)
site.info <- read.csv("site.metadata.csv", header=T)

###################################################################
## step 2: harmonise & aggregate counts using taxon lookup table ##
###################################################################

setwd(input1)
file.list <- list.files()
count = 1

for (i in file.list){
  print(paste("doing", i, count, "of", length(file.list)))
  output <- paste0(outfolder1, "rpp.excl.", i)
  
  #read data
  poll.data<-read.csv(i, header=T, check.names=F)
  
  #harmonise taxa
  poll.harm <- setNames(poll.data, lookup$RPP.LANDCLIMII[match(names(poll.data), lookup$oldvar)])
  
  #write, then re-read, then delete a temporary file
  write.csv(t(poll.harm), "poll.harm.csv")
  poll.harm <- read.csv("poll.harm.csv", header=T)
  file.remove("poll.harm.csv")
  
  #aggregate counts by taxa
  #remove NA columns
  poll.harm <- poll.harm[, colSums(is.na(poll.harm)) != nrow(poll.harm)]
  #aggregate
  poll.agg <- aggregate(. ~ X, poll.harm, sum)
  
  #write out table
  write.table(poll.agg, sep=",", output, row.names=F, col.names=F)
  count = count + 1
} 

###################################################
## step 3: add columns for PPE taxa not present  ##
###################################################

setwd(outfolder1)
file.list <- list.files()
count = 1

#number of taxa are taken from taxon table and used in write.table line

for (i in file.list){
  print(paste("doing", i, count, "of", length(file.list)))
  output.name <- paste0(outfolder2, "ready_", i)
  data <- read.csv(i, header=F)
  poll.joined<-merge(x=taxon, y=data, by = "V1", all=TRUE)
  poll.joined[is.na(poll.joined)] <- 0
  
  write.table(t(poll.joined[1:nrow(taxon),]), output.name, sep=",", row.names = F, col.names = F)
  count = count + 1
}

##############################################
## step 4: aggregate counts by time windows ##
##############################################

setwd(outfolder2)

site.list <- list.files()
count <- 1

for (i in site.list){
  print(paste("doing", i, count, "of", length(file.list)))
  #set up output and input data for site
  site.name <- gsub(".csv", "", i)
  site.name <- gsub("ready_rpp.excl.", "", site.name)
  radius <- site.info$radius[site.info$filename == site.name]
  model <- site.info$model[site.info$filename == site.name]
  output.name <- paste0(outfolder3, "LRAready_", site.name, ".csv")
  
  #read in pollen data
  pollen.uf <- read.csv(i, header=T)
  
  #extract sample ages
  age <- pollen.uf[,1]
  
  #match sample ages to time bins
  ages.binned <- cbind(age, findInterval(age, bins))
  
  #add time bins to original pollen data
  pollen.bound <- as.data.frame(cbind(ages.binned[,2], pollen.uf))
  names(pollen.bound)[1] <- "bins"
  
  #aggregate taxa counts by time windows
  poll.ready <- aggregate(. ~ bins, pollen.bound, sum)
  
  #rename TWs to midpoint values
  TW.code <- as.data.frame(poll.ready[,1])
  TW.code[] <- TWlookup$class[match(poll.ready[,1], TWlookup$TW)]
  
  #bind all together
  poll.ready3 <- cbind(TW.code, radius=radius, model=model, poll.ready[,3:ncol(poll.ready)])
  names(poll.ready3)[1] <- "TW"
  
  #add cols for TW not present in original data
  TWcols <- as.data.frame(TWlookup$class)
  names(TWcols)[1] <- "TW"
  poll.final <- merge(x=TWcols, y=poll.ready3, by = "TW", all=TRUE)
  poll.final[is.na(poll.final)] <- 0
  
  #write results output file
  write.table(t(poll.final[1:25,]), output.name, sep=",", row.names = T, col.names = F)
  count <- count + 1
}

#################################################################
## step 5: remove temporary files (optional!!)                 ##  
## useful to keep until in case checks are needed!!            ##
#################################################################

setwd(wd)
unlink(substr(outfolder1, 1, nchar(outfolder1)-1), recursive = TRUE)
unlink(substr(outfolder2, 1, nchar(outfolder2)-1), recursive = TRUE)

####################################################
## step 6: generate file.lists by grid cell       ##
## only necessary if using REVEALS-in-R in a loop ##
####################################################

#requires site metadata csv to have col LCGRID_ID

setwd(wd)

#check / make folder for storing csv lists of files by gridcell
outfolder <- paste0(wd, "\\file.lists\\")

if(dir.exists(outfolder) == TRUE){
  print("file list directory exists")
} else {
  dir.create(outfolder)
  print("file list directory created")
}

#uses object site.info read in step 1

for(i in unique(site.info$LCGRID_ID)) {
  outfile <- paste0(outfolder, "site.list.", i, ".csv")
  for(row in 1:nrow(site.info)){
    if(i == site.info[row, "LCGRID_ID"]){
      write.table(paste0("LRAready_", site.info[row,2], ".csv"),
                  file=outfile, append=T, quote=F, row.names=F, col.names=F, sep=",")
    } 
  }
}

#generate list of csv grid file lists
write.table(list.files(outfolder), "GC.list.csv", row.names=F, col.names=F, sep=",")

##################
##      END     ##
##################
