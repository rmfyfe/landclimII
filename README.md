# landclimII
code used for preparation and analysis in LANDCLIMII project  
Code is prepared in blocks  
  
RV.file.format.LANDCLIMII.r prepares unformatted pollen count data in csv files to an LRAready csv format  
uses:  
-LANDCLIM.ALLEUROPE.taxa.LUT.080219.csv (synonym table for harmonising taxonomy)  
-bins.csv (used to define start and end points of time windows for aggregating samples)
-run.all.taxa.100219.csv (used to add rows for RPP taxa where not present)
-TWlookup.csv (used to add columns for timewindows where not present)
-site.metadata.csv (used to add site size and type to count data)  
-pollen count data (in csv file format, example datasets provided in folder "input")

