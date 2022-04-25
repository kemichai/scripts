#############################################################################
# Plotting routines
#############################################################################

# plot device stuff
set.square <- function() {
  # set the plot region to be square
  #par(pty="s")
  par(asp=1)
  #fin.min <- min(par()$fin)
  #par(fin=rep(fin.min,2))
  #par(mar=c(0,0,1,0))
  par(mar=c(0,0,0,0))
  invisible()
}
end.square <- function() {
  par(pty="m")
  par(mar=c(5,4,4,2))
  invisible()
}

# low level functions
draw.circle <- function(r=1, centre=c(0,0), npts=101, ...) {
   # draw a circle of radius r centred at c(0,0)
   draw.arc(r=r, centre=centre, npts=npts, theta1=0, theta2=2*pi, ...)
   invisible()
}
draw.arc <- function(r=1, centre=c(0,0), theta1=0, theta2=2*pi, npts=101, ...) {
   # draw an arc of radius r centred at c(0,0) from theta1 to theta2
   theta <- seq(from=theta1, to=theta2, length=npts)
   lines( centre[1]+r*cos(theta), centre[2]+r*sin(theta), ...)
   invisible()
}
draw.ellipse <- function(r=1, a=1, b=1, 
                             thetaa=0,
                             rel.thetaa=T, # theta1 and theta2 measured from thetaa?
                             centre=c(0,0), npts=101, ...) {
  # draw an ellipse with principal axes (a,b) centered
  # at centre with the principal axis at angle thetaa to the horizontal.
  # r is an overall scale
  draw.ellipse.arc(r=r, a=a, b=b, thetaa=thetaa, rel.thetaa=rel.thetaa,
                   centre=centre, npts=npts, ...)
  invisible()
}
draw.ellipse.arc <- function(r=1, a=1, b=1, theta1=0, theta2=2*pi,
                             thetaa=0,
                             rel.thetaa=T, # theta1 and theta2 measured from thetaa?
                             centre=c(0,0), npts=101, ...) {
   theta <- seq(from=theta1, to=theta2, length=npts)
   if(rel.thetaa) theta <- theta + thetaa
   xp <- r*a*cos(theta)
   yp <- r*b*sin(theta)
   lines( centre[1] + cos(thetaa)*xp - sin(thetaa)*yp ,
          centre[2] + sin(thetaa)*xp + cos(thetaa)*yp, ...)
   invisible()
}

draw.box <- function(xcen,ycen,r1=1,r2=1,angle=0,border=NA,col="grey",
                     label=NULL, cex=1.0, ...) {
  # draw a box centred at (xcen,ycen) with sides (r1,r2)
  # rotated clockwise by angle and filled with colour col
  # put text label in the centre - also rotated
  x1 <- c(0,1,1,0,0)-0.5
  y1 <- c(0,0,1,1,0)-0.5
  x2 <- r1*x1
  y2 <- r2*y1
  theta <- pi/180*angle
  amat <- array(c(cos(theta),sin(theta)*c(1,-1),cos(theta)), dim=c(2,2))
  xyp <- amat %*% rbind(x2,y2)
  polygon(xcen+xyp[1,],ycen+xyp[2,],border=border,col=col,...)
  if(!is.null(label)) {
    srt.store <- par()$srt
    par(srt=angle)
    text(xcen,ycen,label=label,cex=cex)
    par(srt=srt.store)
  }
  invisible()
}
#############################################################################
# plotting directions - phi is a bearing (i.e. angle clockwise from North)

plot.phi <- function(x,y,phi,scale=100,add=F,cex=0.4,pch=16,...) {
  # plot points and azimuths
  if(!add) {
    plot(range(x)+scale*c(-1,1), range(y)+scale*c(-1,1), type="n", ...)
  }
  points(x,y, cex=cex,pch=pch)
  vv <- scale*cbind(sin(phi),cos(phi))
  arrows(x-vv[,1],y-vv[,2], x+vv[,1],y+vv[,2], length=0)
  invisible()
}

plot.phi.interval <- function(x,y,phi,phi.interval,
                           scale=100,add=F,cex=0.4,pch=16,
                           fill.col="grey",...) {
  # plot points and azimuths with errors
  dim(phi.interval) <- c(length(phi.interval)/2,2)
  if(!add) {
    plot(range(x)+scale*c(-1,1), range(y)+scale*c(-1,1), type="n", ...)
  }
  # draw interval estimate
  vv1 <- scale*cbind(sin(phi.interval[,1]),cos(phi.interval[,1]))
  vv2 <- scale*cbind(sin(phi.interval[,2]),cos(phi.interval[,2]))
  ppx <- rbind(x-vv1[,1], x+vv1[,1], x+vv2[,1], x-vv2[,1], x-vv1[,1], NA)
  ppy <- rbind(y-vv1[,2], y+vv1[,2], y+vv2[,2], y-vv2[,2], y-vv1[,2], NA)
  polygon(as.vector(ppx), as.vector(ppy), col=fill.col)
  # draw point estimate
  vv <- scale*cbind(sin(phi),cos(phi))
  arrows(x-vv[,1],y-vv[,2], x+vv[,1],y+vv[,2], length=0)
  # draw points
  points(x,y, cex=cex,pch=pch)
  invisible()
}
#############################################################################
# Drawing on circles
start.circle <- function(r=1,rmax=1.4) {
  par(mar=c(1,1,1,1))
  plot(NA,NA,xlim=c(-rmax,rmax),ylim=c(-rmax,rmax),axes=F,
       xlab="",ylab="",asp=1)
  points(0,0,pch="+",cex=1.5)
  draw.arc(r=1,centre=c(0,0),theta1=0,theta2=2*pi)
  invisible()
}
points.circle <- function(theta,r=1,centre=c(0,0), pch=16, ...) {
  points( centre[1]+r*cos(theta), centre[2]+r*sin(theta), pch=pch, ...)
  invisible()
}
radius.circle <- function(theta,r=1,centre=c(0,0), ...) {
  sapply(theta,
         function(theta,r,...) {
           lines(centre[1]+c(0,r*cos(theta)), centre[2]+c(0,r*sin(theta)), ...)
         }, r=r, ...)
  invisible()
}
lines.f.circle <- function(theta,fval,scale=1,r=1,centre=c(0,0), ...) {
  lines( centre[1]+(r+scale*fval)*cos(theta),
         centre[2]+(r+scale*fval)*sin(theta), ...)
  invisible()
}
points.f.circle <- function(theta,fval,scale=1,r=1,centre=c(0,0), pch=16, ...) {
  points( centre[1]+(r+scale*fval)*cos(theta),
          centre[2]+(r+scale*fval)*sin(theta), pch=pch, ...)
  invisible()
}


#############################################################################
# Stereogram
start.stereogram <- function(scale=1.1,convention="geo",plot.centre=T,lwd=2,...) {
   # stereogram: std: x is LEFT, y is UP, z is OUT OF PAPER
   #             geo: x (North) is UP,  y (East) is LEFT, z (Down) is INTO PAPER
   par(mar=c(0,0,0,0))
   blank.plot(scale=scale)
   draw.circle(lwd=lwd,...)
   if(plot.centre) {
      # draw a cross at the centre
      rp <- 0.04
      lines( c(0,0), rp*c(-1,1), lwd=max(1,lwd-1),... )
      lines( rp*c(-1,1), c(0,0), lwd=max(1,lwd-1),... )
   }
   invisible()
}
plot.centre.stereogram <- function(lwd=2,...) {
   # draw a cross at the centre of a stereogram
   rp <- 0.04
   lines( c(0,0), rp*c(-1,1), lwd=max(1,lwd-1),... )
   lines( rp*c(-1,1), c(0,0), lwd=max(1,lwd-1),... )
   invisible()
}
points.xy.stereogram <- function(x, y, convention="geo", pch=16, ...) {
  # draw (x,y) points on a stereogram
  if(convention=="std") {
    points(x,y, pch=pch, ...)
  } else {
    points(y,x, pch=pch, ...) 
  }
  invisible()
}
lines.xy.stereogram <- function(x, y, convention="geo", ...) {
  # draw (x,y) lines on a stereogram
  if(convention=="std") {
    lines(x,y, ...)
  } else {
    lines(y,x, ...) 
  }
  invisible()
}
polygon.xy.stereogram <- function(x, y, convention="geo", ...) {
  # draw a polygon on a stereogram
  if(convention=="std") {
    polygon(x,y, ...)
  } else {
    polygon(y,x, ...) 
  }
  invisible()
}
convert.phiv2.xy.stereogram <- function(phiv2,
                                        method="lambert",
                                        reverse=T) {
  # convert (phi,theta) values to (x,y) coordinates
  dim(phiv2) <- c(length(phiv2)/2,2)
  # check that (0<theta<pi/2) before starting
  idx <- !is.na(phiv2[,2]) & phiv2[,2]>pi/2  # theta>pi/2?
  if(reverse) {
     # reverse direction of any vectors where theta>pi/2
     phiv2[idx,1] <- phiv2[idx,1] + pi
     phiv2[idx,2] <- pi - phiv2[idx,2]
  }
  # now compute locations on the diagram
  if(method=="lambert") { ## method=="lambert"
     r <- sin( phiv2[,2]/2 )*sqrt(2.0) # sin(theta/2)*sqrt(2)
  } else { ## method=="wulff"
     r <- tan( phiv2[,2]/2 ) # tan(theta/2)
  }
  x <- r*cos(phiv2[,1])   # r*cos(phi)
  y <- r*sin(phiv2[,1])   # r*sin(phi)
  xy <- cbind(x,y)
  if(length(xy)==2) dim(xy) <- NULL
  return(xy)
}
convert.xy.phiv2.stereogram <- function(xy) {
  # convert (x,y) coordinates values to phiv2=(phi,theta)
  dim(xy) <- c(length(xy)/2,2)
  phi <- atan2( xy[,2], xy[,1] )
  r <- sqrt(xy[,1]^2 + xy[,2]^2)
  theta <- 2*atan(r)
  phiv2 <- cbind(phi,theta)
  if(length(phiv2)==2) dim(phiv2) <- NULL
  return(phiv2)
}
points.phiv2.stereogram <- function(phiv2, convention="geo", pch=16,
                                    reverse=F, clip=T, ...) {
  # draw(phiv2)=(phi,theta) on a stereogram
  dim(phiv2) <- c(length(phiv2)/2,2)
  ##if(clip) phiv2[phiv2[,2]<=pi/2,] <- NA
  xy <- convert.phiv2.xy.stereogram(phiv2, reverse=reverse)
  dim(xy) <- c(length(xy)/2,2)
  points.xy.stereogram(xy[,1], xy[,2], convention=convention, pch=pch, ...)
  invisible()
}
text.phiv2.stereogram <- function(phiv2, convention="geo", 
                                  reverse=F, clip=T, ...) {
  # draw text at (phiv2)=(phi,theta) on a stereogram
  dim(phiv2) <- c(length(phiv2)/2,2)
  ##if(clip) phiv2[phiv2[,2]<=pi/2,] <- NA
  xy <- convert.phiv2.xy.stereogram(phiv2, reverse=reverse)
  dim(xy) <- c(length(xy)/2,2)
  if(convention=="std") {
     text(x=xy[,1], y=xy[,2], ...)
  } else {
     text(x=xy[,2], y=xy[,1], ...)
  }
  invisible()
}
lines.phiv2.stereogram <- function(phiv2, convention="geo",
                                   no.cross=F, ...) {
  # draw a line specified by (phiv2)=(phi,theta) on a stereogram
  xy <- convert.phiv2.xy.stereogram(phiv2)
  dim(xy) <- c(length(xy)/2,2)
  if(no.cross) {
    # prevent line segments crossing the stereogram from side to
    # side due to wrapping
    # set to NA any (x,y) pairs immediately after or before such cases
    n <- nrow(xy)
    diffxy <- abs(xy[-n,]-xy[-1,])
    dim(diffxy) <- c(n-1,2)
    dxdy <- apply(diffxy, 1, function(x) sqrt(sum(x^2)))
    idx <- c(dxdy>1.5,FALSE)
    nidx <- length(idx[idx])
    if(nidx>0) {
       rr <- ifelse(idx,2,1)
       xy1 <- cbind( rep(xy[,1],rr), rep(xy[,2],rr) ) 
       nxy <- nrow(xy)
       idx1 <- (1:nxy)[idx] + 1:nidx
       xy1[idx1,] <- c(NA,NA)
       xy <- xy1
    }
  }
  lines.xy.stereogram(xy[,1], xy[,2], convention=convention, ...)
  if(convention=="std") {
     return(invisible(xy))  # changed from invisible() ##!!== RA 18-Jun-09
  } else {
     return(invisible(cbind(xy[,2],xy[,1])))
  }
}
radius.stereogram <- function(phi, ...) {
  # draw a radius on a stereogram
  sapply(phi, function(phi,...) {
     lines.phiv2.stereogram(rbind(c(0,0),c(phi,pi/2)), ...)
   }, ...)
  invisible()  
}
polygon.phiv2.stereogram <- function(phiv2, convention="geo",
                                     xcen=0, ycen=0,
                                     gscale=1, xscale=1, yscale=1, 
                                     ...) {
  # draw a polygon specified by (phiv2)=(phi,theta) on a stereogram
  xy <- convert.phiv2.xy.stereogram(phiv2)
  dim(xy) <- c(length(xy)/2,2)
  # NB - note that (x,y) are reversed on a stereogram
  xy[,1] <- ycen + yscale*gscale*xy[,1]
  xy[,2] <- xcen + xscale*gscale*xy[,2]
  polygon.xy.stereogram(xy[,1], xy[,2], convention=convention, ...)
  invisible()
}
points.vecs.stereogram <- function(vecs, convention="geo",
                                   pch=16, clip=T, ...) {
  # draw unit vectors as points on a stereogram
  dim(vecs) <- c(length(vecs)/3,3)
  phiv2 <- extract.phiv2.nvec(vecs)
  dim(phiv2) <- c(length(phiv2)/2,2)
  if(clip) phiv2[vecs[,3]<0,] <- NA
  points.phiv2.stereogram(phiv2, convention=convention, pch=pch, ...)
  invisible()
}
points.rmat.stereogram <- function(rmat, convention="geo",
                                   pch=paste(1:3), cex=1.5, clip=T, revvecs=T,
                                   col=c("red","black"), ...) {
  # draw unit vectors of a rotation matrix as points on a stereogram
  if(revvecs) {
    vecs <- t(cbind(rmat,-rmat))
  } else {
    vecs <- t(rmat)
  }
  if(length(col)==2) {
     col <- ifelse(vecs[,3]>0, col[1], col[2])
  } else if(length(col)==6) {
     col <- col
  } else if(length(col)==3) {
     col <- rep(col,2)
  } else if(length(col)==1) {
     col <- col
  }
  points.vecs.stereogram(vecs, convention=convention,
                         pch=pch, cex=cex, clip=clip, col=col, ...)
  invisible()
}
points.phiv.stereogram <- function(phiv, convention="geo", pch=paste(1:3),
                                   cex=1.5, clip=T, revvecs=F,
                                   col=c("red","black"), ...) {
  # draw unit vectors of rotation matrices specified by phiv
  # as points on a stereogram
  dim(phiv) <- c(length(phiv)/3,3)
  tmp <- apply(phiv, 1,
         function(phiv, convention, pch, cex, clip, revvecs, col, ...) {
           rmat <- rotmat.phiv(phiv)
           points.rmat.stereogram(rmat, convention=convention,
                                  pch=pch, cex=cex, clip=clip, revvecs=revvecs,
                                  col=col, ...)
           return(0)
         }, convention=convention, pch=pch, cex=cex, clip=clip,
            revvecs=revvecs, col=col, ...)
  invisible()
}
points.epar.stereogram <- function(epar, convention="geo",
                                   pch=c("m","2","M"),
                                   cex=1.5, clip=T,
                                   col=c("red","black"), ...) {
  # draw unit vectors of a rotation matrices specified by epar
  # as points on a stereogram
  dim(epar) <- c(length(epar)/4,4)
  phiv <- epar[,-1]
  points.phiv.stereogram(phiv, convention=convention, pch=pch,
                         cex=cex, clip=clip, col=col, ...)
  invisible()
}
points.sdr.stereogram <- function(sdr, convention="geo", pch=paste(1:3),
                                  cex=1.5, clip=T,
                                  col=c("red","black"), ...) {
  # draw unit vectors of a rotation matrices specified by sdr
  # as points on a stereogram
  phiv <- convert.sdr.phiv(sdr)
  points.phiv.stereogram(phiv, convention=convention,
                         pch=pch, cex=cex, clip, col=col, revvecs=T, 
                         ...)
  invisible()
}
lines.vecs.stereogram <- function(vecs, convention="geo", clip=T, ...) {
  # draw lines connecting a set of unit vectors on a stereogram
  dim(vecs) <- c(length(vecs)/3,3)
  phiv2 <- extract.phiv2.nvec(vecs)
  dim(phiv2) <- c(length(phiv2)/2,2)
  if(clip) phiv2[vecs[,3]<0,] <- NA
  lines.phiv2.stereogram(phiv2, convention=convention, ...)
  invisible()
}
polygon.vecs.stereogram <- function(vecs, convention="geo", clip=F,
                                    xcen=0, ycen=0,
                                    gscale=1, xscale=1, yscale=1, 
                                    ...) {
  # draw polygon using a set of unit vectors on a stereogram
  dim(vecs) <- c(length(vecs)/3,3)
  phiv2 <- extract.phiv2.nvec(vecs)
  dim(phiv2) <- c(length(phiv2)/2,2)
  if(clip) phiv2[vecs[,3]<0,] <- NA
  polygon.phiv2.stereogram(phiv2, convention=convention,
                           xcen=xcen, ycen=ycen,
                           gscale=gscale, xscale=xscale, yscale=yscale, 
                           ...)
  invisible()
}
contour.stereogram <- function(phivec, thetavec, z, convention="geo",
                               no.cross=T, 
                               add=F, nlevels=10, levels=NULL,
                               zrange=NULL, plot.centre=F,
                               ...) {
  # Draw a contour plot for a stereogram
  if(!add) start.stereogram(convention=convention, plot.centre=plot.centre)
  if(!is.null(zrange)) {
    z[z<zrange[1]] <- zrange[1]
    z[z>zrange[2]] <- zrange[2]
  }
  if(is.null(levels)) levels <- pretty(range(z,na.rm=T))
  ccl <- contourLines(phivec, thetavec, z, nlevels=nlevels, levels=levels)
  cclxy <- lapply(ccl, function(c1, convention, no.cross, ...) {
      lines.phiv2.stereogram(cbind(c1$x,c1$y), convention=convention,
                             no.cross=no.cross, ...)
    }, convention=convention, no.cross=no.cross, ...)
  invisible(cclxy) # changed from invisible() ##!!== RA 18-Jun-09
}
draw.great.circle.arc.stereogram <- function(vec1,vec2, bvec1,bvec2, do.plot=T,
                                             npts=101, convention="geo",
                                             zero.tol=1.e-15, ...) {
  # draw a great circle arc from vec1 to vec2 on a stereogram
  # where bvec1, bvec2 form an orthogonal basis for the great circle
  lambda1 <- atan2( dot(vec1,bvec2), dot(vec1,bvec1) )
  lambda2 <- atan2( dot(vec2,bvec2), dot(vec2,bvec1) )
  if(abs(lambda2-lambda1)>pi) {
    if(lambda1<lambda2) {
      lambda1 <- lambda1+2*pi
    } else {
      lambda2 <- lambda2+2*pi
    }
  }
  #lambda.test <- 0.5*(lambda1+lambda2)
  #vec.test <- cos(lambda.test)*bvec1 + sin(lambda.test)*bvec2
  #if(vec.test[3]<0) {
  #  if(lambda1<lambda2) {
  #    lambda1 <- lambda1+2*pi
  #  } else {
  #    lambda2 <- lambda2+2*pi
  #  }
  #}
  lambda <- seq(from=lambda1, to=lambda2, length=npts)
  vecs <- outer(cos(lambda),bvec1) + outer(sin(lambda),bvec2)
  # catch the case where vec1=-vec2 and vecs has gone negative
  if(max(abs(vec1+vec2))<zero.tol && mean(vecs[,3])<0) vecs <- -vecs[npts:1,]
  # eliminate any roundoff errors:
  vecs[1,] <- vec1
  vecs[npts,] <- vec2
  # draw the arc
  if(do.plot) lines.vecs.stereogram(vecs, convention=convention, ...)
  invisible(vecs)
}
draw.great.circle.stereogram <- function(phiv2, convention="geo",
                                         npts=101, clip=T, ...) {
  # Draw a great circle with a pole at (phi,theta)=phiv2
  # unit vector towards the pole

  # Make sure the pole has theta<pi/2, otherwise fix it:
  if(phiv2[2]>pi/2) {
    phiv2[1] <- phiv2[1]+pi
    phiv2[2] <- pi-phiv2[2]
  }
  # Plot the pole
  #points.phiv2.stereogram(phiv2,convention=convention,col="red")

  nvec <- nvec.phiv2(phiv2)
  # scale this outwards to the edge: svec lies in the xy plane
  svec <- c(nvec[1:2],0)
  nn.svec <- sqrt(sum(svec^2))
  if(nn.svec==0) {
    svec <- c(1,0,0)
    nn.svec <- 1
  }
  svec <- svec/nn.svec
  
  # form two unit vectors perpendicular to nvec
  uvec <- cross.product(nvec, svec)
  nn.uvec <- sqrt(sum(uvec^2))
  if(nn.uvec==0) {
     uvec <- cross.product(nvec, c(0,1,0))
     nn.uvec <- sqrt(sum(uvec^2))
  }
  uvec <- uvec/nn.uvec
  vvec <- cross.product(nvec,uvec)
  # set of points on a great circle
  ptheta1 <- ifelse(clip,    0,    0)
  ptheta2 <- ifelse(clip,   pi, 2*pi)
  ptheta <- seq(from=ptheta1, to=ptheta2, length=npts)
  vecs <- outer(cos(ptheta),uvec) + outer(sin(ptheta),vvec)
  lines.vecs.stereogram(vecs, convention=convention, ...)
  invisible(vecs)
}
stereogram.grid <- function(type=2,ng=9,nf=3,col="grey",add=T, ...) {
  # add a grid to a stereogram
  # ignores add=T
  if(type==1) {
     idx <- c(1,2,3)
     tmin <- 1
     phig <- seq(from=0, to=2*pi, length=4*ng+1)
     phip <- seq(from=0, to=2*pi, length=4*nf*ng+1)
     thetag <- seq(from=0, to=pi/2, length=ng+1)
     thetap <- seq(from=0, to=pi/2, length=nf*ng+1)
  } else if(type==2) {
     idx <- c(3,1,2)
     tmin <- 2
     phig <- seq(from=0, to=pi, length=2*ng+1)
     phip <- seq(from=0, to=pi, length=2*nf*ng+1)
     thetag <- seq(from=0, to=pi, length=2*ng+1)
     thetap <- seq(from=0, to=pi, length=2*nf*ng+1)
  } else {
     idx <- c(3,1,2)
     tmin <- 2
     phig <- seq(from=0, to=pi, length=2*ng+1)
     phip <- seq(from=0, to=pi, length=2*nf*ng+1)
     thetag <- seq(from=0, to=pi, length=2*ng+1)
     thetap <- seq(from=0, to=pi, length=2*nf*ng+1)
  }
  phi <- phip
  for(i in 2:(length(thetag)-1)) {
     theta <- thetag[i]
     v1 <- cbind( sin(theta)*cbind( cos(phi),sin(phi) ), cos(theta) )[,idx]
     lines.vecs.stereogram(v1, col=col, ...)
  }
  theta <- thetap
  for(i in tmin:(length(phig)-1)) {
     phi <- phig[i]
     v1 <- cbind(sin(theta)*cbind( rep(cos(phi),length(thetap)),
                                   rep(sin(phi),length(thetap))),
                 cos(theta) )[,idx]
     lines.vecs.stereogram(v1, col=col, ...)
  }
  invisible()   
}
#############################################################################
# Beach Ball Plot
beach.ball <- function(sdr,
                       xcen=0, ycen=0, gscale=1, xscale=1, yscale=1, 
                       convention="geo", add=F,
                       col="black", zero.tol=1.e-15,
                       plot.centre=F, npts=101,
                       ...) {
  if(!add) {
    start.stereogram(convention=convention,plot.centre=plot.centre)
  } else {
    draw.circle(r=gscale, centre=c(xcen,ycen))
  }
  rmat1 <- rotmat.sdr(sdr) # = [u,a,n]
  # make sure that avec and nvec (vectors 2 and 3) in this matrix
  # point upwards: i.e. no negative signs in the third row in cols 2 and 3
  # first make sure nvec=rmat1[,3] points upwards
  #if(rmat1[3,3]<0) print("flip 1,3")
  if(rmat1[3,3]<0) rmat1 <- rmat1%*%invaxis.2()
  # now check that avec=rmat1[,2] points upwards
  # if not then exchange axes 1 and 3
  #if(rmat1[3,2]<0) print("exchange 1,3")
  if(rmat1[3,2]<0) rmat1 <- rmat1%*%exchaxis.2()
  
  # note sign of uvec
  uvec.up <- rmat1[3,1]>=0
  # note whether avec is horizontal
  avec.horizontal <- is.zero(rmat1[3,2],zero.tol)
  # horizontal unit vector to boundary in direction of nvec
  b.nvec <- c(rmat1[1:2,3],0)
  ctmp <- sqrt(sum(b.nvec^2))
  if(ctmp==0) {
    b.nvec <- rmat[,3]
  } else {
    b.nvec <- b.nvec/ctmp
  }
  # horizontal unit vector to boundary in direction of uvec
  b.uvec <- c(rmat1[1:2,1],0)
  ctmp <- sqrt(sum(b.uvec^2))
  if(ctmp==0) {
    b.uvec <- rmat[,1]
  } else {
    b.uvec <- b.uvec/ctmp
  }
    
  # calculate boundary points
  bmat1 <- apply(rmat1, 2,
              function(rvec) {
                if(all(rvec==c(0,0,1))) {
                  uvec <- c(1,0,0)
                } else {
                  uvec <- c(rvec[2],-rvec[1],0)/sqrt(sum(rvec[1:2]^2))
                }
                return(uvec)
              })
  
  # calculate the third perpendicular vector
  vmat1 <- apply(rbind(rmat1,bmat1), 2,
                 function(rb) {
                   vvec <- cross.product(rb[1:3], rb[4:6])
                 })

  # Identify points for line drawing
  rvec <- rmat1[,2] # crossing point avec (2)
  b1 <- +bmat1[,1]  # boundary point on circle perp to uvec (1)
  b2 <- -bmat1[,3]  # boundary point on circle perp to nvec (3)
  # arcs between points:
  # crossing point (avec) to boundary on circle perp to uvec (1)
  v1 <- draw.great.circle.arc.stereogram(rvec,  b1,  bmat1[,1],vmat1[,1], 
                                         do.plot=F, npts=npts)
  if(avec.horizontal && !uvec.up) {
     # along circle edge from circle perp to uvec to circle perp to nvec
     # need to make sure we go via the correct points
     if(all(is.zero(rvec-b1,zero.tol))) {
        v2.1 <- draw.great.circle.arc.stereogram(
                                        b1,  -b.nvec,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2.2 <- draw.great.circle.arc.stereogram(
                                       -b.nvec,  b2,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2 <- rbind(v2.1,v2.2)
     } else if(all(is.zero(rvec-b2,zero.tol))) {
        v2.1 <- draw.great.circle.arc.stereogram(
                                        b1,  b.uvec,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2.2 <- draw.great.circle.arc.stereogram(
                                        b.uvec,  b2,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2 <- rbind(v2.1,v2.2)
     } else {
        v2 <- draw.great.circle.arc.stereogram(
                                    b1,  b2,  c(1,0,0), c(0,1,0),
                                    do.plot=F, npts=npts)
     }
  } else {
     # along circle edge from circle perp to uvec to circle perp to nvec
     v2 <- draw.great.circle.arc.stereogram(
                                 b1,  b2,  c(1,0,0), c(0,1,0),
                                 do.plot=F, npts=npts)
  }
  # back from boundary to avec along circle perp to nvec (3)
  v3 <- draw.great.circle.arc.stereogram(
                              b2,rvec,  bmat1[,3],vmat1[,3],
                              do.plot=F, npts=npts)
  # shade region
  poly.pts <- rbind(v1,v2,v3)
  polygon.vecs.stereogram(poly.pts, convention=convention,
                          xcen=xcen, ycen=ycen,
                          gscale=gscale, xscale=xscale, yscale=yscale, 
                          col=col, border=T, lwd=2)

  # arcs between points:
  # crossing point (avec) to boundary on circle perp to uvec (1)
  v1 <- draw.great.circle.arc.stereogram(
                              rvec,  -b1,  bmat1[,1],vmat1[,1],
                              do.plot=F, npts=npts)
  if(avec.horizontal && !uvec.up) {
     # along circle edge from circle perp to uvec to circle perp to nvec
     # need to make sure we go via the correct points
     if(all(is.zero(rvec+b1,zero.tol))) {
        v2.1 <- draw.great.circle.arc.stereogram(
                                       -b1,  -b.nvec,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2.2 <- draw.great.circle.arc.stereogram(
                                       -b.nvec, -b2,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2 <- rbind(v2.1,v2.2)
     } else if(all(is.zero(rvec+b2,zero.tol))) {
        v2.1 <- draw.great.circle.arc.stereogram(
                                       -b1,  b.uvec,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2.2 <- draw.great.circle.arc.stereogram(
                                        b.uvec, -b2,  c(1,0,0), c(0,1,0),
                                        do.plot=F, npts=npts)
        v2 <- rbind(v2.1,v2.2)
     } else {
        v2 <- draw.great.circle.arc.stereogram(
                                -b1,  -b2,  c(1,0,0),c(0,1,0),
                                 do.plot=F, npts=npts)
     }
  } else {
     v2 <- draw.great.circle.arc.stereogram(
                                  -b1,  -b2,  c(1,0,0),c(0,1,0),
                                   do.plot=F, npts=npts)
  }
  # back from boundary to avec along circle perp to nvec (3)
  v3 <- draw.great.circle.arc.stereogram(
                             -b2, rvec,  bmat1[,3],vmat1[,3],
                              do.plot=F, npts=npts)
  # shade region
  poly.pts <- rbind(v1,v2,v3)
  polygon.vecs.stereogram(poly.pts, convention=convention,
                          xcen=xcen, ycen=ycen,
                          gscale=gscale, xscale=xscale, yscale=yscale, 
                          col=col, border=T, lwd=2)

  invisible()
}

#############################################################################
# plots on the unit sphere

start.unitsphere <- function(alpha.view=55*pi/180,beta.view=45*pi/180) {
   # start plot
   psize <- 1.2
   tsize <- 1.5
   asize <- 1.4
   par(pty="s",xpd=NA)
   plot(psize*c(-1,1),psize*c(-1,1), typ="n", axes=F, xlab="", ylab="",asp=1)
   # draw axes
   xyz <- asize*diag(3)
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   apply(pxpypz[,1:2], 1, function(pxpy) lines(c(0,pxpy[1]),c(0,pxpy[2]),
                                               lwd=2))
   xyz <- tsize*diag(3)
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   # label axes
   text(pxpypz[1,1], pxpypz[1,2], labels="X", cex=1.5)
   text(pxpypz[2,1], pxpypz[2,2], labels="Y", cex=1.5)
   text(pxpypz[3,1], pxpypz[3,2], labels="Z", cex=1.5)
   # draw crosses on axes
   cxl <- 0.08
   a1 <- c(1,0,0,1,0,0)
   a2 <- cxl*c(0,-1,0,0,1,0)
   a3 <- cxl*c(0,0,-1,0,0,1)
   idx <- c(c(3,1,2),3+c(3,1,2))
   for(i in 1:3) {
     a1 <- a1[idx]; a2 <- a2[idx]; a3 <- a3[idx]
     lines.unitsphere(t(array(a1+a2, dim=c(3,2))),
                      clip=F,alpha.view=alpha.view,beta.view=beta.view,lwd=2)
     lines.unitsphere(t(array(a1+a3, dim=c(3,2))),
                      clip=F,alpha.view=alpha.view,beta.view=beta.view,lwd=2)
   }
   
   draw.circle(lwd=2)
   par(pty="m")
   invisible()
}
end.unitsphere <- function(alpha.view=55*pi/180,beta.view=45*pi/180) {
  par(pty="m")
  invisible()
}
xyz.to.pxpypz <- function(xyz,alpha.view=55*pi/180,beta.view=45*pi/180) {
  dim(xyz) <- c(length(xyz)/3,3)
  sin.alpha.view <- sin(alpha.view)
  cos.alpha.view <- cos(alpha.view)
  sin.beta.view <- sin(beta.view)
  cos.beta.view <- cos(beta.view)
  rmat <- array(c( -sin.beta.view,cos.beta.view,0,
   -cos.alpha.view*cos.beta.view,-cos.alpha.view*sin.beta.view,sin.alpha.view,
    sin.alpha.view*cos.beta.view, sin.alpha.view*sin.beta.view,cos.alpha.view),
                  dim=c(3,3))
  pxpypz <- t( t(rmat)%*%t(xyz) )
  if(length(pxpypz)==3) dim(pxpypz) <- NULL
  return(pxpypz)
}
contour.unitsphere <- function(phivec, thetavec, ftab, 
                               add=F, nlevels=10,
                               clip=T,
                               alpha.view=55*pi/180,beta.view=45*pi/180,
                               ...) {
  # Draw a contour plot for a unitsphere
  if(!add) start.unitsphere(alpha.view,beta.view)
  ccl <- contourLines(phivec,thetavec,ftab,nlevels=nlevels)
  lapply(ccl, function(c1,clip,alpha.view,beta.view,...) {
      lines.phi.theta.unitsphere(c1$x,c1$y, clip=clip,
                             alpha.view=alpha.view, beta.view=beta.view, ...)
    }, clip=clip,
       alpha.view=alpha.view, beta.view=beta.view, ...)
  invisible()
}

lines.unitsphere <- function(xyz,clip=T,
                             alpha.view=55*pi/180,beta.view=45*pi/180,...) {
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   if(clip) {
     idx <- (pxpypz[,3]<0)
     pxpypz[idx,1] <- NA
     pxpypz[idx,2] <- NA
   }
   lines(pxpypz[,1],pxpypz[,2],...)
   invisible()
}
lines.phi.theta.unitsphere <- function(phi,theta,clip=T,
                                   alpha.view=55*pi/180,beta.view=45*pi/180,
                                   ...) {
   xyz <- cbind(sin(theta)*cbind(cos(phi),sin(phi)),cos(theta))
   lines.unitsphere(xyz,clip=clip,alpha.view=alpha.view,beta.view=beta.view,
                    ...)
   invisible()
}
points.unitsphere <- function(xyz,clip=T,
                              alpha.view=55*pi/180,beta.view=45*pi/180,...) {
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   dim(pxpypz) <- c(length(pxpypz)/3,3)
   if(clip) {
     idx <- (pxpypz[,3]<0)
     pxpypz[idx,1] <- NA
     pxpypz[idx,2] <- NA
   }
   points(pxpypz[,1],pxpypz[,2],...)
   invisible()
}
points.phi.theta.unitsphere <- function(phi,theta,clip=T,
                                    alpha.view=55*pi/180,beta.view=45*pi/180,
                                    ...) {
   xyz <- cbind(sin(theta)*cbind(cos(phi),sin(phi)),cos(theta))
   points.unitsphere(xyz,clip=clip,alpha.view=alpha.view,beta.view=beta.view,
                     ...)
   invisible()
}
points.phiv2.unitsphere <- function(phiv2,clip=T,
                                    alpha.view=55*pi/180,beta.view=45*pi/180,
                                    ...) {
   points.phi.theta.unitsphere(phiv2[,1],phiv2[,2],clip=T,
                                    alpha.view=55*pi/180,beta.view=45*pi/180,
                                    ...) 
   invisible()
}
points.phiv.unitsphere <- function(phiv,
                                   pch=rep(16,3),
                                   col=c("red","blue","green"),
                                   lwd=c(1,1,1),
                                   add=T, ...) {

  # plot the directions of a set of right-hand-sets
  dim(phiv) <- c(length(phiv)/3,3)
  if(!add) start.unitsphere()
  if(length(pch)==1) pch <- rep(pch,3)
  if(length(col)==1) col <- rep(col,3)
  if(length(lwd)==1) lwd <- rep(lwd,3)
  re <- rotmat.phiv.asrow(phiv)
  pv1 <- extract.phi.theta.unitvec(re[,  1:3])
  pv2 <- extract.phi.theta.unitvec(re[,3+1:3])
  pv3 <- extract.phi.theta.unitvec(re[,6+1:3])

  points.phi.theta.unitsphere(pv1[,1],pv1[,2],
                          col=col[1], pch=pch[1], lwd=lwd[1], ...)
  points.phi.theta.unitsphere(pv2[,1],pv2[,2],
                          col=col[2], pch=pch[2], lwd=lwd[2], ...)
  points.phi.theta.unitsphere(pv3[,1],pv3[,2],
                          col=col[3], pch=pch[3], lwd=lwd[3], ...)
  invisible()
}
vector.unitsphere <- function(xyz,clip=T,
                              alpha.view=55*pi/180,beta.view=45*pi/180,
                              length=0.12, lwd=2, angle=25, ...) {
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   dim(pxpypz) <- c(length(pxpypz)/3,3)
   if(clip) {
     idx <- (pxpypz[,3]<0)
     pxpypz[idx,1] <- NA
     pxpypz[idx,2] <- NA
   }
   zero <- rep(0,dim(pxpypz)[1])
   arrows(zero, zero, pxpypz[,1], pxpypz[,2],
          length=length, lwd=lwd, angle=angle, ...)
   invisible()
}
polygon.unitsphere <- function(xyz,clip=T,
                               alpha.view=55*pi/180,beta.view=45*pi/180,...) {
   pxpypz <- xyz.to.pxpypz(xyz,alpha.view=alpha.view,beta.view=beta.view)
   if(clip) {
     idx <- (pxpypz[,3]<0)
     pxpypz[idx,1] <- NA
     pxpypz[idx,2] <- NA
   }
   polygon(pxpypz[,1],pxpypz[,2],...)
   invisible()
}
spin.contour.unitsphere <- function(phivec, thetavec, ftab,
                        nstep=30, ncycle=2, speed=30,
                        alpha.view=55*pi/180, beta.view=45*pi/180,
                        under=F, slight=F, ...) {
   # spin a unitsphere plot
   par(mfrow=c(1,1))
   if(slight) alpha.view <- 80*pi/180
   if(under) alpha.view <- pi-alpha.view
   for(i in 0:(ncycle*nstep)) {
      contour.unitsphere(phivec, thetavec, ftab, 
                         alpha.view=alpha.view,
                         beta.view=beta.view+2*pi*i/nstep, ...)
      system(paste("sleep",1/speed))
   }
   invisible()
}
spin.unitsphere <- function(func, nstep=30, ncycle=2, speed=30,
                            alpha.view=55*pi/180, beta.view=45*pi/180,
                            under=F, slight=F, ask=F, ...) {
   # spin a unitsphere plot
   save.ask <- par()$ask
   par(ask=ask)
   par(mfrow=c(1,1))
   if(slight) alpha.view <- 80*pi/180
   if(under) alpha.view <- pi-alpha.view
   for(i in 0:(ncycle*nstep)) {
      func(alpha.view=alpha.view, beta.view=beta.view+2*pi*i/nstep, ...)
      system(paste("sleep",1/speed))
   }
   par(ask=save.ask)
   invisible()
}
#############################################################################
angelier.plot <- function(sdr,add=F,type=1,
                          legend=T,cex.legend=0.7,col="black") {
  # make plots like Figure 3a and 3b of Angelier 1982
  if(!add) start.stereogram()
  dim(sdr) <- c(length(sdr)/3,3)
  if(is.character(type) || type==2) {
    # plot principal axes
    if(is.character(type)) {
      type <- tolower(type)
    } else {
      type <- c("u","n","b","p","t")
    }
    type[type=="s"] <- "u"
    type[type=="a"] <- "b"
    apply(sdr,1,function(sdr1) {
         rmat1 <- rotmat.sdr(sdr1)
         if(any(type=="s")) points.vecs.stereogram(rbind(rmat1[,1],-rmat1[,1]),pch=3,col=col)  # slip
         if(any(type=="n")) points.vecs.stereogram(rbind(rmat1[,3],-rmat1[,3]),pch=10,col=col) # fault normal
         if(any(type=="b")) points.vecs.stereogram(rbind(rmat1[,2],-rmat1[,2]),pch=8,col=col)  # B axis
         paxis <- normalise(rmat1[,1]-rmat1[,3])
         if(any(type=="p")) points.vecs.stereogram(rbind(paxis,-paxis),pch=15,col=col) # P axis
         taxis <- normalise(rmat1[,1]+rmat1[,3])
         if(any(type=="t")) points.vecs.stereogram(rbind(taxis,-taxis),pch=22,col=col) # T axis
       })
    if(legend) {
      ii <- match(type, c("u","n","b","p","t"))
      pch1 <- c(3,10,8,15,22)[ii]
      leglab <- c("slip","normal","B","P","T")[ii]
      usr <- par()$usr
      par(usr=c(0,1,0,1))
      legend(1,0, xjust=1,yjust=0,
             pch=pch1, 
             legend=leglab, 
             pt.cex=cex.legend,cex=cex.legend,bty="n")
      par(usr=usr)
    }
  } else if(type==1) {
    # plot fault planes as arcs, and show the slip vectors
    apply(sdr,1,function(sdr1) {
         rmat1 <- rotmat.sdr(sdr1)
         if(rmat1[3,1]<0) rmat1 <- cbind(-rmat1[,1],rmat1[,2],-rmat1[,3])
         uvec <- rmat1[,1]
         bvec <- cross.product(rmat1[,3],c(0,0,1))
         bvec <- bvec/sqrt(sum(bvec^2))
         points.vecs.stereogram(uvec,col=col)
         draw.great.circle.arc.stereogram(bvec,-bvec, rmat1[,1],rmat1[,2],
                                          col=col)
       })
  } else {
    stop("type not recognised")
  }
  invisible()
}
#############################################################################
