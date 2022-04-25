#-----------------------------------------------------------------------------
# processing a set of .csv files with strike/dip/rake/error
# build and load routines
source("first.q")
reload()

##for(cc in c("k1","k2","k3","k4","k1.3d","k2.3d","k3.3d","k4.3d")) {
##for(cc in c("k1","k2","k3","k4")) {###,"k1.3d","k2.3d","k3.3d","k4.3d")) {
##for(cc in c("k1.3d","k2.3d","k3.3d","k4.3d")) {###,"k1.3d","k2.3d","k3.3d","k4.3d")) {
for(cc in c("c1","c2","c3","c4","c5","c6")) {
  dirname <- paste("../stress.monitoring/data/",cc,"/",sep="")
  ##filenames <- system(paste("ls ",dirname,"sdr_cluster_",cc,"_????",sep=""),intern=T)
  filenames <- system(paste("ls ",dirname,"sdr_stress_cluster_",cc,"_????",sep=""),intern=T)
  filenames <- sapply(strsplit(filenames,dirname),function(x) x[2])
  ncluster <- length(filenames)
  # define file locations
  # data files
  inputdir <- dirname
  # output files
  outputdir <- dirname
  # lookup tables used by the stress inversion routines
  tabdir <- "tab/"

  nfiles <- length(filenames)
  for(i in 1:nfiles) {
  
    # Read the file in standard format with strike/dip/rake/error in degrees
    # Get file name
    fin <- filenames[i]
    cat(paste("Reading file",fin,"at",date(),"..."))
    # Full path to input file
    ff <- paste(inputdir,fin,sep="")
    # Extract filename stem - to build other filenames later
    ####fname <- unlist(strsplit(fin,split=".",fixed=T))[1]
    fname <- filenames[i]
    # Read the data - skip=1 means skip one row before reading the data
    df <- read.table(ff,skip=1,header=F,sep=",",)
    names(df) <- c("strike","dip","rake","err")
    cat(paste(nrow(df),"earthquakes read.\n"))
    
    # Do the inversion - this creates a number of output files in outputdir
    # NB - set save=F to stop R dumping all of the output objects
    #      (the output text files are probably enough for most purposes)
    stress.out <- stress.inversion(fin,fname=fname,
                                   inputdir=inputdir,outputdir=outputdir,
                                   tabdir=tabdir,
                                   verbose=T,save=F,opar=opar)
  }

}

for(cc in (c("c1","c2"))) {
  for(pp in c("e","l")) {
    
    dirname <- paste("../stress.monitoring/data/canty_",cc,"_",pp,"/",sep="")
    filenames <- system(paste("ls ",dirname,"sdr_stress_cluster_canty_",cc,"_",pp,"_????",sep="")
                        ,intern=T)
    filenames <- sapply(strsplit(filenames,dirname),function(x) x[2])
    ncluster <- length(filenames)
    # define file locations
    # data files
    inputdir <- dirname
    # output files
    outputdir <- dirname
    # lookup tables used by the stress inversion routines
    tabdir <- "tab/"

  nfiles <- length(filenames)
  for(i in 1:nfiles) {
  
    # Read the file in standard format with strike/dip/rake/error in degrees
    # Get file name
    fin <- filenames[i]
    cat(paste("Reading file",fin,"at",date(),"..."))
    # Full path to input file
    ff <- paste(inputdir,fin,sep="")
    # Extract filename stem - to build other filenames later
    ####fname <- unlist(strsplit(fin,split=".",fixed=T))[1]
    fname <- filenames[i]
    # Read the data - skip=1 means skip one row before reading the data
    df <- read.table(ff,skip=1,header=F,sep=",",)
    names(df) <- c("strike","dip","rake","err")
    cat(paste(nrow(df),"earthquakes read.\n"))
    
    # Do the inversion - this creates a number of output files in outputdir
    # NB - set save=F to stop R dumping all of the output objects
    #      (the output text files are probably enough for most purposes)
    stress.out <- stress.inversion(fin,fname=fname,
                                   inputdir=inputdir,outputdir=outputdir,
                                   tabdir=tabdir,
                                   verbose=T,save=F,opar=opar)
  }

}}

#-----------------------------------------------------------------------------
