LANDCLIMII project r-scripts  
code used for preparation and analysis in LANDCLIMII project  
  
Code is prepared in blocks  
  
*RV.file.format.LANDCLIMII.r* prepares unformatted pollen count data in csv files to an LRAready csv format  
requires:  
-LANDCLIM.ALLEUROPE.taxa.LUT.080219.csv (synonym table for harmonising taxonomy)  
-bins.csv (used to define start and end points of time windows for aggregating samples)  
-run.all.taxa.100219.csv (used to add rows for RPP taxa where not present)  
-TWlookup.csv (used to add columns for timewindows where not present)  
-site.metadata.csv (used to add site size and type to count data)  
-pollen count data (in csv file format, example datasets provided in folder "input")  

*Kunes_REVEALS_iterate_across_sitelist.r* controls the application of the REVEALS in R function (Abraham et al 2014: REF1) to LRA ready files, by iterating across groups of sites (grouped into 1 degree grid cells for LANDCLIMII project). Selection of grid cells controlled by a list of all grid cells.  
requires (all in the working directory):  
-*LRA-REVEALS-release.R* (Kunes' REVEALS-in-R script, available at REF2)  
-LRAready_sitename.csv files (product of *RV.file.format.LANDCLIMII.r* script)  
-site.list.xxx.csv files (one per grid cell, xxx as grid cell name/code)  
-GC.list.csv (list of all site.list.xxx.csv files)  
-taxa.csv (list of RPP taxa for formatting output files)  
-av-RPP.all.landclimII.csv (RPP and fall speed of pollen for taxa)  
-avar-RPP.all.landclimII.csv (variance-covariance matrix of RPP values)  
  
outputs csv files for  
(a) grid cells where rows = taxa and cols = time windows (mean and standard errors separately)  
(b) time windows where rows = grid cells and cols = taxa REVEALS estimates (mean and standard error separately)  

REFERENCES  
REF1: Abraham, V., Oušková, V., & Kuneš, P. 2014. Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. PLoS ONE 9: e100117.  
REF2: https://github.com/petrkunes/LRA  

