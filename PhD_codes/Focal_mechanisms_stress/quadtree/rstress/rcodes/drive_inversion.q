#-----------------------------------------------------------------------------
# processing a set of .csv files with strike/dip/rake/error
# build and load routines
source("first.q")
reload()
#-----------------------------------------------------------------------------
# EDIT CODE IN THIS SECTION --------------------------------------------------
#-----------------------------------------------------------------------------

# define file locations
# data files
#inputdir <- "indata/temporal_1/"
#inputdir <- "indata/temp_overlap/"
inputdir <- "indata/"

# output files
#outputdir <- "outdata_temporal_1"
outputdir <- "outdata/"
# lookup tables used by the stress inversion routines
tabdir <- "tab/"

# list of earthquake files to process
# each input file has a header row, and then strike/dip/rake/error
# one per row, all values in degrees
# error=error of strike/dip/rake in degrees (e.g. 10)
#filenames <- c("crete.csv","dixie.csv")



# clusters
filenames <- c("cluster_3.csv",
"cluster_4.csv",
"cluster_5.csv",
"cluster_6.csv",
"cluster_8.csv",
"cluster_9.csv",
"cluster_10.csv",
"cluster_11.csv",
"cluster_12.csv",
"cluster_13.csv",
"cluster_14.csv",
"cluster_15.csv",
"cluster_16.csv")

#filenames <- c("cRed1Green2.csv","cYellow3Blue5.csv","cOrange4.csv","cPurple6.csv","cCyan7.csv","cGrey8.csv","cPink9.csv")

#filenames <- c("Cluster1red.csv","Cluster2pink.csv","Cluster3yellow.csv","Cluster4green.csv","Cluster5cyan.csv","Cluster6blue.csv","Cluster7purple.csv","Cluster8brown.csv","Cluster9seagreen.csv","Cluster10grey.csv","Cluster11lightblue.csv","Cluster12orange.csv","Cluster13black.csv","Cluster14tan.csv","Cluster15coral.csv","Cluster16rosybrown.csv","Cluster17plum.csv","Cluster18darkorange.csv","Cluster19magenta.csv","Cluster20chocolate.csv")
#filenames <-
#c("red1.csv","orange2.csv","yellow3.csv","green4.csv","cyan5.csv","blue6.csv","purple7.csv","brown8.csv","grey9.csv")

#filenames <-
#c("01red.csv","02green.csv","03yellow.csv","04orange.csv","05blue.csv","06purple.csv","07cyan.csv","08grey.csv","09brown.csv","10pink.csv")

#filenames <-
#c("west_revision.csv","basin_revision.csv","east_revision.csv")

#filenames <-
#c("Lauramay.csv")
#filenames <-
#c("Western.csv","CentralBasin.csv","Eastern.csv","Basin_Eastern.csv","All_Regions.csv")

#filenames <-
#c("above_plate_cluster.csv","upperhalf_plate_cluster.csv","lowerhalf_plate_cluster.csv","below_plate_cluster.csv")

#filenames <-
#c("170less_cluster.csv","170more_cluster.csv")

#filenames <-
#c("170less_lowerhalf_cluster.csv", "170less_upperhalf_cluster.csv")

#filenames <-
#c("200less_lowerhalf_cluster.csv")

#filenames <-
#c("upperhalf_more10picks.csv", "lowerhalf_more10picks_less170.csv")

#filenames <-
#c("upperhalf_5pick33err.csv", "lowerhalf_5pick33err_170less.csv", "lowerhalf_5pick33err.csv")


# keep this value F (FALSE) unless a lot of post processing information
# is needed (the full pdf) -- it takes a lot of space.
save.all <- T
#-----------------------------------------------------------------------------
# END OF CODE TO EDIT --------------------------------------------------------
#-----------------------------------------------------------------------------

# number of files to process
nfiles <- length(filenames)
for(i in 1:nfiles) {

  # Read the file in standard format with strike/dip/rake/error in degrees
  # Get file name
  fin <- filenames[i]
  cat(paste("Reading file",fin,"at",date(),"..."))
  # Full path to input file
  ff <- paste(inputdir,fin,sep="")
  # Extract filename stem - to build other filenames later
  fname <- unlist(strsplit(fin,split=".",fixed=T))[1]
  # Read the data - skip=1 means skip one row before reading the data
  df <- read.table(ff,skip=1,header=F,sep=",",)
  names(df) <- c("strike","dip","rake","err")
  cat(paste(nrow(df),"earthquakes read.\n"))
  
  # Do the inversion - this creates a number of output files in outputdir
  # NB - set save=F to stop R dumping all of the output objects
  #      (the output text files are probably enough for most purposes)
  stress.out <- stress.inversion(fin,
                                 inputdir=inputdir,outputdir=outputdir,
                                 tabdir=tabdir,
                                 verbose=T,save=save.all,opar=opar)
}

#-----------------------------------------------------------------------------
