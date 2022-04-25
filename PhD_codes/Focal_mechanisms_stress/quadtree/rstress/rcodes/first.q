#############################################################################
# Run this file first - everything else follows
# source("first.q"); reload()
#############################################################################
reload <- function(object.list=F,
                   source.list=c("first.q","opar.q",
                                 "stress_util.q","stress_plot.q",
                                 "stress_func.q","stress_data.q",
                                 "stress_jtfn.q","stress_devl.q"),
                   package.list=c("boot"),  ###c("adapt","boot")
                   fortran.lib.dir="./",
                   fortran.libname=c("stress_f77","libadapt"),
                   make.flags="FFLAGS=-O3\\ -Wuninitialized\\ -g",
                   recompile=T) {
   windows <- .Platform$OS.type=="windows"
   llextension <- .Platform$dynlib.ext
   if(windows && !(fortran.lib.dir=="./" || fortran.lib.dir=="")) {
     stop("Under windows all .f files must be in the current directory")
   }
   if(windows) fortran.lib.dir <- ""
   cat("Required packages...\n")
   for(pname in package.list) {
     cat(paste("Package...",pname,"..."))
     eval(parse(text=paste("require(",pname,")")))
     cat("...OK\n")
   }
   cat("Reloading files...\n")
   for(fname in source.list) {
      cat(paste("R source code...",fname,"..."))
      eval(parse(text=paste("source(\"",fname,"\")",sep="")))
      cat("...OK\n")
   }
   if(recompile) {
      cat("Compiling Fortran...\n")
      if(!is.null(make.flags)) {
        Sys.setenv(MAKEFLAGS=make.flags)
        cat(paste("MAKEFLAGS =",Sys.getenv("MAKEFLAGS"),"\n"))
      }
      cmd <- paste("R CMD SHLIB ",
                   paste(fortran.lib.dir,fortran.libname,".f",
                         sep="",collapse=" "),
                   sep="")
      print(cmd)
      system(cmd)
   }
   if(object.list) {
     cat("List of exported objects:\n")
     if(windows) {
       system(paste("nm -g ",fortran.lib.dir,fortran.libname,".o",sep=""))
     } else {
       system(paste("nm -g ",fortran.lib.dir,fortran.libname,llextension,sep=""))
     }
   }
   cat("Loading Fortran...\n")
   #for(object in fortran.additional.libs) dyn.load(object)
   if(windows) {
      # a single .dll is created under windows
      libn <- paste(fortran.lib.dir,fortran.libname[1],llextension,sep="")
      try(dyn.unload(libn),silent=T)
      dyn.load(libn)
    } else {
      # under unix there are several - each to be loaded separately
      # apparently this is no longer a difference between windows
      # and unix - builds a single .so to load
      libs <- paste(fortran.lib.dir,fortran.libname,llextension,sep="")
      for(libn in libs[1]) {
        try(dyn.unload(libn),silent=T)
        print(libn)
        dyn.load(libn)
      }
   }
   #cat("Making Michael inversion programs...\n")
   #make.michael.cc(srcdir="michael-algorithm/src/",bindir="bin/")
   invisible()
}

#############################################################################
