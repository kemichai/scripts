#############################################################################
# Functions written by JT
#############################################################################
single.page.summary <- function(invresult,psummary,
                                store.stereonet.file=NULL,opar=NULL) {
  # Make a single-page summary of the posterior distribution

  # First extract coordinate vectors
  cat("Extracting coordinate vectors\n")
  nuvec <- invresult$epar.vecs$nu
  phivec <- invresult$epar.vecs$phi
  thetavec <- invresult$epar.vecs$theta
  psivec <- invresult$epar.vecs$psi
  mean.rmat <- as.earth(psummary$mean.epar)
  mean.rmat <- mean.rmat$rmat.gs
  map.rmat  <- as.earth(psummary$map.epar)
  map.rmat <- map.rmat$rmat.gs
  
  # Extract axis positions
  cat("Extracting axis positions\n")
  dd <- list(NULL,NULL,NULL)
  dd[[1]] <- invresult
  dd[[2]] <- retabulate.qflist(c(0,0,0),dd[[1]],exchange=1)	# Retab for S2
  dd[[3]] <- retabulate.qflist(c(0,0,0),dd[[1]],exchange=2)	# Retab for S3

  # Compute z-axis density
  cat("Computing z-axis density\n")
  ddz <- list(NULL,NULL,NULL)
  ddz[[1]] <- apply(dd[[1]]$qvals0,c(2,3),sum)	# Marginalise over nu, psi (S1)
  ddz[[2]] <- apply(dd[[2]]$qvals0,c(2,3),sum)	# Ditto (S2)
  ddz[[3]] <- apply(dd[[3]]$qvals0,c(2,3),sum)	# Ditto (S3)
  for(i in (1:3)) {
     # smooth each function in phi
     # grid size in phi and theta
     nphi <- dd[[i]]$nepar[2]
     ntheta <- dd[[i]]$nepar[3]
     # grid spacing in phi
     dphi <- dd[[i]]$epar.vecs$phi[2]-dd[[i]]$epar.vecs$phi[1]
     # for each row scale according to the scale phi.smooth.scale
     phi.smooth.scale <- dphi/2 #= half of grid spacing at the equator
     s <- trunc( phi.smooth.scale/dphi/sin(dd[[i]]$epar.vecs$theta) )
     s <- pmin(s,trunc(nphi/2)) # s=half-width of triangular boxcar filter
     # for the first and last columns take the mean to enforce continuity
     ddz[[i]][,1] <- mean(ddz[[i]][,1])
     ddz[[i]][,ntheta] <- mean(ddz[[i]][,ntheta])
     # now smooth each column
     ddz[[i]] <- array.column.boxcar(ddz[[i]],s,wrap=T)
  }

  # Produce stereogram
  cat("Producing stereogram\n")
  par(opar)
  par(mfrow=c(2,2))
  par(mar=c(0.5,0.5,0.5,0.5))
  cols <- c("red","green","blue")               # Define plotting colours
                                                # for max,mid,min
  stress.stereonet(invresult,psummary,
                   draw.grid=T,draw.map.rmat=T,
                   cols=cols,
                   store.stereonet.file=store.stereonet.file)
  
##   ststore <- list(NA,NA,NA)
##   ####postscript(file="output.ps",horizontal=F,onefile=T,width=8,height=8)
##   start.stereogram(plot.centre=T)               # Construct basemap
##   stereogram.grid()                             # Add grid
##   ststore[[1]] <- contour.stereogram(phivec,thetavec,ddz[[1]],add=T,col=cols[1]) # S1
##   ststore[[2]] <- contour.stereogram(phivec,thetavec,ddz[[2]],add=T,col=cols[2]) # S2
##   ststore[[3]] <- contour.stereogram(phivec,thetavec,ddz[[3]],add=T,col=cols[3]) # S3
##   # Note that in the following we must interpret the first (third) column
##   # of the rotation matrix as S3 (S1) and colour it accordingly
##   points.rmat.stereogram(mean.rmat[,3],pch=19,col=cols[1])
##   points.rmat.stereogram(mean.rmat[,2],pch=19,col=cols[2])
##   points.rmat.stereogram(mean.rmat[,1],pch=19,col=cols[3])
##   points.rmat.stereogram(map.rmat[,3],pch=21,col=cols[1],bg="white")
##   points.rmat.stereogram(map.rmat[,2],pch=21,col=cols[2],bg="white")
##   points.rmat.stereogram(map.rmat[,1],pch=21,col=cols[3],bg="white")
##   radius.stereogram(psummary$mean.phi.shmax,col="black",lwd=1,lty="dashed")
##   radius.stereogram(psummary$mean.phi.shmax+pi,col="black",lwd=1,lty="dashed")

##   # Write out the stereogram contour line segments
##   if(!is.null(store.stereonet.file)) {
##     # store stereogram file
##     for(i in 1:3) {
##       # break up any line segments with NA values before writing
##       ofile <- paste(store.stereonet.file,"_S",i,".dat",sep="")
##       cat(paste("Writing",ofile,"\n"))
##       sink(ofile)
##       ststore[[i]] <- delist1(lapply(ststore[[i]],split.na))
##       lapply(ststore[[i]], function(x) {
##          cat("####\n")
##          x <- matrix(x,ncol=2)
##          write.table(x,append=T,sep=",",row.names=F,col.names=F)
##          invisible()
##       })
##       sink()
##     }
##   }
  
  # Illustrate azimuths
  par(mar=c(5,4,1.5,1.5))
  azphi.density <- psummary$azphi.density
  shmax.density <- psummary$phi.shmax
  x <- shmax.density$x*180/pi
  y <- shmax.density$y
  y <- y/mean(y)
  ymax <- max(y)
  ymax <- max(unlist(sapply(azphi.density,function(aa) max(aa$y/mean(aa$y)))))
  plot(x,y,col="black",type="l",yaxs="i",
  	xlim=c(0,180),ylim=c(0,1.1*ymax),
  	ylab="Posterior density",xlab="Azimuth")
  for (i in 1:3) {
    x <- azphi.density[[i]]$x*180/pi
    y <- azphi.density[[i]]$y
    y <- y/mean(y)
    lines(x,y,col=cols[4-i])
  }

  # Write on key
  plot(NA,NA,xlab="",ylab="",axes=F,xlim=c(0,1),ylim=c(0,1))
  legend(0,1.0,xjust=0,yjust=1.0,
         col=c(cols,"black"),lty=c(1,1,1,2),bty="n",
         legend=c(expression(S["max"]),expression(S["mid"]),
                  expression(S["min"]),expression(S["Hmax"])))
  legend(0,0,xjust=0,yjust=0,
         pch=c(19,21),pt.bg="white",col="black",
         legend=c("Posterior mean","Posterior mode"),bty="n")

  # Illustrate stress ratio
  plot(nuvec,psummary$nu.density,type="l",yaxs="i",
       ylim=1.1*c(0,max(psummary$nu.density)),
       xlab=expression(nu),
       ylab="Posterior density")
  
  ###dev.off()

  par(opar)
  invisible()
}


bivariate.plots <- function(psummary) {
  # Now do bivariate plots
  cat("Making bivariate plots\n")
  par(mfrow=c(6,6),mar=c(4,4,0.5,0.5))
  phivec <- phivec*180/pi
  thetavec <- thetavec*180/pi
  psivec <- psivec*180/pi
  
  # X-Y contour plot of phi vs. theta
  cat("    phi vs. theta\n")
  zz <- apply(dd[[1]]$qvals0,c(2,3),sum)
  par(mfg=c(1,4))
  contour(phivec,thetavec,zz,
  	xlab=expression(phi),ylab=expression(theta),
  	cex.lab=1.5)
  # X-Y contour plot of psi vs. theta
  cat("    psi vs. theta\n")
  zz <- apply(dd[[2]]$qvals0,c(4,3),sum)
  par(mfg=c(1,5))
  contour(psivec,thetavec,zz,
  	xlab=expression(psi),ylab=expression(theta),
          cex.lab=1.5)
  # X-Y contour plot of nu vs. theta
  cat("    nu vs. theta\n")
  zz <- apply(dd[[3]]$qvals0,c(1,3),sum)
  par(mfg=c(1,6))
  contour(nuvec,thetavec,zz,
          xlab=expression(nu),ylab=expression(theta),
          cex.lab=1.5)
  # X-Y contour plot of psi vs. phi
  cat("    psi vs. phi\n")
  zz <- apply(dd[[1]]$qvals0,c(4,2),sum)
  par(mfg=c(2,5))
  contour(psivec,phivec,zz,
          xlab=expression(psi),ylab=expression(phi),
          cex.lab=1.5)
  # X-Y contour plot of nu vs. phi
  cat("    nu vs. phi\n")
  zz <- apply(dd[[2]]$qvals0,c(1,2),sum)
  par(mfg=c(2,6))
  contour(nuvec,phivec,zz,
          xlab=expression(nu),ylab=expression(phi),
          cex.lab=1.5)
  # X-Y contour plot of nu vs. psi
  cat("    nu vs. psi\n")
  zz <- apply(dd[[3]]$qvals0,c(1,4),sum)
  par(mfg=c(3,6))
  contour(nuvec,psivec,zz,
          xlab=expression(nu),ylab=expression(psi),
          cex.lab=1.5)

  invisible()
}
#############################################################################
