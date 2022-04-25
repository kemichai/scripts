#############################################################################
# Functions for earthquake analysis
#############################################################################
# util06.q
#############################################################################
# required libraries
#require()

#############################################################################
# Functions for creating and viewing eps files

gsview <- function(filename) {
   gvcmd <- ifelse(.Platform$OS.type=="windows","gsview32","gv")
   system(paste(gvcmd,filename))
}
ask.to.continue <- function() {
  cat("Hit return to continue: ")
  invisible(scan(quiet=T))
}

open.eps <- function(epsfile, width=7.0, aspect=1, height=NULL) {
   # start an .eps file
   if(is.null(height)) height <- width*aspect
   postscript(file=epsfile, onefile=F, paper="special",
              horizontal=F, width=width, height=height)
   invisible(epsfile)
}
close.eps <- function(epsfile=NULL,view=F) {
   gvcmd <- ifelse(.Platform$OS.type=="windows","gsview32","gv")
   # close an .eps file
   if(names(dev.cur())=="postscript") {
      dev.off()
      if(view) system(paste(gvcmd,epsfile))
   } else {
      stop("Current device is not postscript")
   }
   invisible()
}
view.eps <- function(epsfile) {
   gvcmd <- ifelse(.Platform$OS.type=="windows","gsview32","gv")
   system(paste(gvcmd,epsfile))
}
blank.plot <- function(scale=1, asp=1,
                       x=scale*c(-1,1), y=scale*c(-1,1), ...) {
   plot(x=x,y=y,type="n",
        axes=F,xlab="",ylab="",asp=asp,...)
   invisible()
}
   

optimiser.message <- function(object, method="nlm", look=F) {
  if(look) {
    if(!is.null(object$optim)) {
      object <- object$optim
      method <- "optim"
    } else if(!is.null(object$nlm)) {
      object <- object$nlm
      method <- "nlm"
    } else {
      object <- NULL
      method <- ""
    }
  }
  if(method=="optim") {
    val <- object$convergence
    retval <- switch(paste(val),
          "0"="indicates successful convergence.",
          "1"="indicates that the iteration limit 'maxit' had been reached.",
         "10"="indicates degeneracy of the Nelder-Mead simplex.",
         "51"="indicates a warning from the \"L-BFGS-B\" method;\n see component \'message\' for further details.",
         "52"="indicates an error from the \"L-BFGS-B\" method;\n see component \'message\' for further details."
                     )
  } else if(method=="nlm") {
    val <- object$code
    retval <- switch(paste(val),
          "1"="relative gradient is close to zero, current iterate is\nprobably solution.",
          "2"="successive iterates within tolerance, current iterate is\n probably solution.",
          "3"="last global step failed to locate a point lower than\n\'estimate\'.  Either \'estimate\' is an approximate local\nminimum of the function or \'steptol\' is too small.",
          "4"="iteration limit exceeded.",
          "5"="maximum step size 'stepmax' exceeded five consecutive\ntimes.  Either the function is unbounded below, becomes\nasymptotic to a finite value from above in some\ndirection or 'stepmax' is too small."
             )
  } else {
    val <- ""
    retval <- paste("Method",method,"not known.")
  }
  retval <- paste("Convergence of",method,":",val,"\n",retval)
  return(retval)
}
#############################################################################
# Messing with lists and other stuff
split.na <- function(xy) {
  # split an [n x 2] array at rows where there are NA values
  # outputting a list
  nxy <- nrow(xy)
  idx <- (1:nxy)[apply(xy,1,function(x) any(is.na(x)))]
  if(length(idx)>0) {
     ii <- cbind(c(1,idx+1), c(idx-1,nxy))
     xy <- as.list(by(ii,1:nrow(ii),function(ii) xy[ii[1,1]:ii[1,2],]))
   } else {
     xy <- list(xy)
   }
   return(xy)
}
delist1 <- function(xx) {
  # take a 2 level list, and make it a 1 level list
  rx <- list()
  i <- 0
  for(j in 1:length(xx)) {
    if(is.list(xx[[j]])) {
       for(k in 1:length(xx[[j]])) {
         i <- i+1
         rx[[i]] <- xx[[j]][[k]]
       }
    } else {
      i <- i+1
      rx[[i]] <- xx[[j]]
    }
  }
  return(rx)
}
#############################################################################
# using Fortran functions

# utilities
f77.ftest <- function(x) {
   n <- 2
   nfun <- 1
   retval <- rep(0,nfun)
   flist <- .Fortran("ftest", as.integer(n),
                              as.double(x),
                              as.integer(nfun),
                              as.double(retval))
   return(flist[[4]])
}
f77.sfloor <- function(x) {
   retval <- 0
   flist <- .Fortran("sfloor", as.double(x), as.integer(retval))
   return(flist[[2]])
}

#############################################################################
# vector and matrix operations

# vector products
dot <- function(a,b) sum(a*b)
cross.product <- function(a,b) c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])

# normalise a unit vector
normalise <- function(x) x/sqrt(sum(x^2))

# tests on matrices
is.zero <- function(x,zero.tol=0) {
  return(all(abs(x)<=zero.tol))
}
is.orthogonal <- function(rmat, zero.tol=1.0e-15) {
  # returns true if the matrix rmat is orthogonal
  all( abs(rmat%*%t(rmat)-diag(c(1,1,1))) <= zero.tol )
}
is.symmetric <- function(rmat, zero.tol=1.0e-15) {
  # returns true if the matrix rmat is symmetric
  all( abs(rmat-t(rmat)) <= zero.tol )
}
sum.mat.trap.single <- function(a,i) {
  # sum matrix a over (single) index i using the trapezium rule
  d <- dim(a)
  dk <- (1:length(d))[-i]
  d0 <- d[i]  # size of dimension i
  d1 <- d[-i] # size of other dimensions
  ix <- rep("",length(d))
  ix[i] <- 1
  ix1 <- paste(ix,collapse=",")
  ix[i] <- d0
  ix2 <- paste(ix,collapse=",")
  eval(parse(text=paste("a0 <- 0.5*(a[",ix1,"]+a[",ix2,"])",sep="")))
  a1 <- apply(a,dk,sum)
  return(a1-a0)
}
sum.mat.trap <- function(a,margin=NULL) {
  # sum matrix a over all indices except margin using the trapezium rule
  d <- dim(a)
  if(is.null(d)) d <- length(a)
  if(is.null(margin)) {
    idx <- 1:length(d)
    if(length(idx)==0) {
      return(a)
    } else {
      return(sum(scale.mat.edges(a,idx)))
    }
  } else {
    idx <- (1:length(d))[-margin]
    if(length(idx)==0) {
      return(a)
    } else {
      return(apply(scale.mat.edges(a,idx),margin,sum))
    }
  }
}
scale.mat.edges <- function(a,idx=NULL,scale=0.5) {
  # multiply the first and last entries of indices idx by 0.5
  d <- dim(a)
  if(is.null(d)) d <- length(a)
  if(is.null(idx)) idx <- 1:length(d)
  for(i in idx) {
     d0 <- d[i]
     ix <- rep("",length(d))
     ix[i] <- 1
     ix1 <- paste(ix,collapse=",")
     ix[i] <- d0
     ix2 <- paste(ix,collapse=",")
     eval(parse(text=paste("a[",ix1,"] <- scale*a[",ix1,"]",sep="")))
     eval(parse(text=paste("a[",ix2,"] <- scale*a[",ix2,"]",sep="")))
  }
  return(a)
}

# orthonormalisation
gram.schmidt.3d <- function(rmat,idx=1:3) {
  # orthonormalise all vectors (columns of rmat)
  # assumes 3x3 matrix
  normalise <- function(x) x/sqrt(sum(x^2))
  
  rmat[,idx[1]] <- normalise(rmat[,idx[1]])
  rmat[,idx[2]] <-(rmat[,idx[2]]
                   - dot(rmat[,idx[1]],rmat[,idx[2]])*rmat[,idx[1]])
  rmat[,idx[2]] <- normalise(rmat[,idx[2]])
  rmat[,idx[3]] <-(rmat[,idx[3]]
                   - (dot(rmat[,idx[1]],rmat[,idx[3]])*rmat[,idx[1]]
                     +dot(rmat[,idx[2]],rmat[,idx[3]])*rmat[,idx[2]]) )
  rmat[,idx[3]] <- normalise(rmat[,idx[3]])
  return(rmat)
}

# Fortran versions
# Scalar product
f77.scalar.product <- function(a,b) {
   n <- length(a)
   if(length(b)!=n) stop("Vectors must be the same length")
   adotb <- 0
   flist <- .Fortran("scalar_product_",
                     as.double(a), as.double(b),
                     as.integer(n),
                     as.double(adotb))
   return(flist[[4]])
}
f77.mat.transpose <- function(x) {
   n <- nrow(x)
   m <- ncol(x)
   tmat <- array(0,dim=c(m,n))
   flist <- .Fortran("mat_transpose_", as.double(x),
                                       as.integer(n),
                                       as.integer(m),
                                       as.double(tmat))
   tmat <- array(flist[[4]], dim=c(m,n))
   return(tmat)
}   
f77.sqmat.transpose.inplace <- function(x) {
   n <- sqrt(length(x))
   flist <- .Fortran("sqmat_transpose_inplace_", as.double(x), as.integer(n))
   x <- array(flist[[1]], dim=c(n,n))
   return(x)
}   
f77.matrix.multiply <- function(a,b) {
  n1 <- dim(a)[1]
  n12 <- dim(a)[2]
  if(dim(b)[1]!=n12) stop("Arrays not conformable")
  n2 <- dim(b)[2]
  outmat <- array(0,dim=c(n1,n2))
  flist <- .Fortran("matrix_multiply_",
                    as.double(a), as.double(b),
                    as.integer(n1), as.integer(n12), as.integer(n2),
                    as.double(outmat))
  outmat <- array(flist[[6]],dim=c(n1,n2))
  return(outmat)
}

sqrt.matrix <- function(amat) {
  # square root of a (symmetric, positive definite) matrix
  ee <- eigen(amat)
  idx <- rev(order(ee$values))
  rmat <- ee$vectors[,idx]
  dmat <- diag(sqrt(ee$values[idx]))
  smat <- rmat %*% dmat %*% t(rmat)
  return(smat)
}
polar.decomposition <- function(amat) {
  # polar decomposition of a matrix
  bmat <- t(amat) %*% amat
  kmat <- sqrt.matrix(bmat)
  mmat <- amat %*% solve(kmat)
  return(list(polar=mmat,elliptical=kmat))
}

#############################################################################
# Averaging directional data

mean.azimuth <- function(phi,w=NULL) {
  # mean of azimuthal data
  # (optional weights w)
  n <- length(phi)
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  atan2( sum(w*sin(phi)), sum(w*cos(phi)) )
}

mean.azimuth.axis <- function(phi,w=NULL) {
  # mean of axial data
  # (optional weights w)
  n <- length(phi)
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  0.5*atan2( sum(w*sin(2*phi)), sum(w*cos(2*phi)) )
}
mean.phiv2 <- function(phiv2,w=NULL) {
  # mean of spherical polar coordinates (phi,theta)
  # (optional weights w)
  n <- length(phiv2)/2
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  mvec <- apply(nvec.phiv2(phiv2),2,
                function(x,w) sum(w*x), w=w)
  extract.phiv2.nvec(normalise(mvec))
}
mean.phiv2.axis <- function(phiv2,w=NULL) {
  # mean of axial spherical polar coordinates (phi,theta)
  # (optional weights w)
  n <- length(phiv2)/2
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  nvecs <- nvec.phiv2(phiv2)
  amat <- t(nvecs) %*% diag(w) %*% nvecs
  ee <- eigen(amat)
  i <- match(max(ee$values),ee$values)
  extract.phiv2.nvec(ee$vectors[,i])
}
mean.phiv <- function(phiv,w=NULL) {
  # mean of Euler angle coordinates (phi,theta,psi)
  # (optional weights w)
  n <- length(phiv)/3
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  smat <- rotmat.phiv.asrow(phiv)
  amat <- array(apply(smat,2,function(x,w) sum(w*x),w=w),dim=c(3,3))
  rmat <- polar.decomposition(amat)$polar
  if(det(rmat)<0) rmat[,3] <- -rmat[,3]
  extract.phiv.rotmat(rmat)
}
mean.phiv.axis <- function(phiv,w=NULL) {
  # mean of axial Euler angle coordinates (phi,theta,psi)
  # (optional weights w)
  n <- length(phiv)/3
  dim(phiv) <- c(n,3)
  if(is.null(w)) {
    w <- 1/n*rep(1,n)
  } else {
    w <- w/sum(w)
  }
  smat <- rotmat.phiv.asrow(phiv)
  xmat <- t(apply(smat,1,function(x) {
                          dim(x) <- c(3,3)
                          retval <- c( x[,1]%*%t(x[,1]), 
                                       x[,2]%*%t(x[,2]),
                                       x[,3]%*%t(x[,3]) )
                          return(retval)
                        }))
  xmat <- array(w,dim=c(n,27))*xmat
  amat <- array(apply(xmat,2,sum),dim=c(3,3,3))
  ymat <- apply(amat,3,function(x) {
                           ee <- eigen(x)
                           idx <- rev(order(ee$values))[1]
                           return(ee$vectors[,idx])
                         })
  rmat <- polar.decomposition(ymat)$polar
  if(det(rmat)<0) rmat[,3] <- -rmat[,3]
  extract.phiv.rotmat(rmat)
}
mean.sdr <- function(sdr,w=NULL) {
  # mean of strike, dip, rake values
  # properly taking into account the fourfold ambiguities of sdr
  # (optional weights w)
  n <- length(sdr)/3
  dim(sdr) <- c(n,3)
  phiv <- convert.sdr.phiv(sdr)
  smat <- rotmat.phiv.asrow(phiv)
  
  # convert rmat=[u a n] matrices to qmat=[t a p]
  # postmultiply rmat by rqmat to make qmat:   qmat = rmat%*%rqmat
  rqmat <- rotmat.b2mat(-pi/4)
  # convert to [t a p] matrices, and extract the Euler angles
  phiv.tap <- t(apply(smat,1,function(rmat) {
                               dim(rmat) <- c(3,3)
                               qmat <- rmat%*%rqmat
                               phiv.tap <- extract.phiv.rotmat(qmat)
                               return(phiv.tap)
                             }))
  # take the mean of these Euler angles as axes
  mphiv.tap <- mean.phiv.axis(phiv.tap, w=w)
  # convert back to rmat=[u a n]
  qmat <- rotmat.phiv(mphiv.tap)
  rmat <- qmat%*%t(rqmat)
  
  # extract Euler angles
  mphiv <- extract.phiv.rotmat(rmat)
  # convert to sdr
  msdr <- convert.phiv.sdr(mphiv)
  return(msdr)
}

#############################################################################
locate.in.array <- function(condition) {
  # find the locations of the elements of a condition array
  # condition is an array of TRUE/FALSE values
  dimx <- dim(condition)
  if(is.null(dimx)) dimx <- length(condition)
  p <- length(dimx)
  m <- c(1,sapply(1:p, function(i) prod(dimx[1:i])))
  n <- prod(dimx)
  ii <- (1:n)[as.vector(condition)]
  ni <- length(ii)
  if(ni==0) {
    idx <- NULL
  } else {
    idx <- t(sapply(ii, function(j,m,p) {
                     ret <- 1 + (j-1)%/%m[-(p+1)] - dimx*((j-1)%/%m[-1])
                     return(ret)
                  },m=m,p=p))
  }
  return(idx)
}
getq.csum <- function(x,csum,probs=c(0,0.025,0.05,0.1,0.5,0.9,0.95,0.975,1)) {
  # get quantiles from a cumulative sum
  n <- length(x)
  x <- c(0.5*(x[-1]+x[-n]),x[n]) # adjust half a bin upwards
  ii <- pmin(n-1,sapply(probs,
                        function(p,n,csum) {
                          if(p==0 || p<csum[1]) {
                            return(1)
                          } else if(p==1) {
                            return(n)
                          } else {
                            max((1:n)[csum<=p])
                          }
                        }, n=n, csum=csum))
  q <- ifelse((csum[ii+1]-csum[ii])==0 | probs==0,x[ii],
          x[ii] + (probs-csum[ii])*(x[ii+1]-x[ii])/(csum[ii+1]-csum[ii]))
  names(q) <- paste(100*probs,"%",sep="")
  return(q)
}

#############################################################################
# Averaging angular data
mean.theta <- function(theta,weights=1,k=1) {
  # Average angles around a circle
  # k=1 for directions, k=2 for axes
  theta <- theta*k
  weights <- weights/sum(weights)
  ctheta <- sum(cos(theta)*weights)
  stheta <- sum(sin(theta)*weights)
  r <- sqrt(ctheta^2+stheta^2)
  return(list(r=r,                            # mean resultant
              cvar=1-r,                       # circular variance
              theta=atan2(stheta,ctheta)/k))  # mean direction
}

median.theta <- function(theta,weights=1) {
  circular.wquantile(theta,weights,probs=0.5)
}
#############################################################################
# Quantiles
circular.wquantile <- function(theta,weights=1,probs=seq(0,1,0.25)) {
  # Circular quantiles using weights
  # NB -- if weights=a pdf, make sure it sums to 1
  n <- length(theta)
  if(length(weights)==1) weights <- rep(weights,n)
  # sort into order
  odx <- order(theta)
  theta <- theta[odx]
  weights <- weights[odx]

  # Calculate the circular median (Refer MJ, p30, eq. 3.4.18)
  # Minimise E[ pi-|pi-|(theta-phi)|| ] over phi: then phi=median
  ww <- sapply(theta,
               function(phi) {
                 sum(weights*(pi-abs(pi-abs((theta-phi)%%(2*pi)))))
               })
  idx <- which.min(ww)
  if(all(probs==0.5)) {
    # only the median is requested
    qvals <- rep(theta[idx],length=length(probs))
  } else {
    offset <- theta[idx]-pi
    theta <- (theta-offset)%%(2*pi)
    odx <- order(theta)
    theta <- theta[odx]
    weights <- weights[odx]

    # calculate the quantiles
    qvals <- wquantile(theta,weights=weights,probs=probs)
    qvals <- (qvals+offset)%%(2*pi)
  }
  names(qvals) <- paste(100*probs,"%",sep="")
  return(qvals)
}
wquantile <- function(x,weights=1,probs=seq(0,1,0.25)) {
  # Linear quantiles using weights
  n <- length(x)
  if(length(weights)==1) weights <- rep(weights,n)
  fval <- weights
  sfval <- sum(weights)
  odx <- order(x)
  x <- x[odx]
  fval <- fval[odx]
  cfval <- cumsum(fval)
  rfval <- sfval-cfval+fval[n]

  dd1 <- outer(cfval,sfval*probs,"-")
  dd2 <- outer(rfval,sfval*(1-probs),"-")
  dd1[dd1<0] <- Inf
  dd2[dd2<0] <- Inf
  jdx1 <- apply(dd1,2,which.min)
  jdx2 <- apply(dd2,2,which.min)
  qvals <- ifelse(fval[jdx1]+fval[jdx2]>0,
          (fval[jdx1]*x[jdx1]+fval[jdx2]*x[jdx2])/(fval[jdx1]+fval[jdx2]),
                  0.5*(fval[jdx1]+fval[jdx2]))
  #print(jdx1); print(x[jdx1]); print(fval[jdx1])
  #print(jdx2); print(x[jdx2]); print(fval[jdx2])
  names(qvals) <- paste(100*probs,"%",sep="")
  return(qvals)
}

# Testing for angular differences between two axial distributions
# on the unit circle
ddiff.axial <- function(f1,f2) {
  # f1, f2 densities -- tabulated from 0 to pi
  # should have matching values at 0 and pi -- ensure this
  m <- length(f1)-1
  f1[m] <- 0.5*(f1[1]+f1[m+1]); f1[1] <- f1[m]
  f2[m] <- 0.5*(f2[1]+f2[m+1]); f2[1] <- f2[m]
  f1 <- f1/sum(f1[-1])
  f2 <- f2/sum(f2[-1])
  fu <- sapply(0:m, function(k,m) sum( f1[1+(1:m)%%m] * f2[1+(1:m-k)%%m] ), m=m)
  fu <- fu/sum(fu[-1])
  return(fu)
}
ddiff.axial.quantile <- function(theta,f1,f2,probs=c(0.025,0.500,0.975)) {
  # Find quantiles
  # theta = vector from 0 to pi
  fu <- ddiff.axial(f1,f2)
  hqvals <- circular.wquantile(2*theta[-1],fu[-1],probs=c(0.5,probs))
  offset <- hqvals[1]-pi # hqvals[1]=median
  qvals <- 0.5*((hqvals-offset)%%(2*pi))+offset/2
  qvals <- qvals[-1]
  names(qvals) <- paste(100*probs,"%",sep="")
  return(qvals)
}


#############################################################################
# Euler angles and matrices

#-----------------------------------------------------------------------------
# unit vectors
unitvec.phi.theta <- function(phi.theta) {
  # unit vector with colatitude azimth phi=phi.theta[1] and theta=phi.theta[2]
  dim0 <- dim(phi.theta)
  n <- length(phi.theta)/2
  dim(phi.theta) <- c(n,2)
  st <- sin(phi.theta[,2]) # sin(theta)
  ct <- cos(phi.theta[,2]) # cos(theta)
  retval <- cbind(st*cos(phi.theta[,1]),
                  st*sin(phi.theta[,1]),
                  ct)
  if(n==1) dim(retval) <- NULL
  return(retval)
}
extract.phi.theta.unitvec <- function(u) {
  # return (phi,theta) from a set of unit vectors
  dim0 <- dim(u)
  n <- length(u)/3
  if(n==1) dim(u) <- c(n,3)
  retval <- cbind( atan2(u[,2],u[,1]), acos(u[,3]) )
  if(is.null(dim0)) dim(retval) <- NULL
  return(retval)
}

#-----------------------------------------------------------------------------
# elementary transformation matrices
rotmat.b1mat <- function(theta) {
  # Euler matrix B1
  ct <- cos(theta)
  st <- sin(theta)
  return(matrix(c(  1,  0,  0,
                    0, ct,-st,
                    0, st, ct),byrow=T,nrow=3))
}
rotmat.b2mat <- function(theta) {
  # Euler matrix B2
  ct <- cos(theta)
  st <- sin(theta)
  return(matrix(c( ct,  0, st,
                    0,  1,  0,
                  -st,  0, ct),byrow=T,nrow=3))
}
rotmat.b3mat <- function(theta) {
  # Euler matrix B3
  ct <- cos(theta)
  st <- sin(theta)
  return(matrix(c( ct,-st,  0,
                   st, ct,  0,
                    0,  0,  1),byrow=T,nrow=3))
}

# Permutation Matrices
exchaxis.1 <- function() {
  # Permutation matrix C1: exchange axes 2,3 reverse 1
  return( matrix( c(-1, 0, 0,
                     0, 0, 1,
                     0, 1, 0),byrow=T,nrow=3))
}
exchaxis.2 <- function() {
  # Permutation matrix C2: exchange axes 1,3 reverse 2
  return( matrix( c( 0, 0, 1,
                     0,-1, 0,
                     1, 0, 0),byrow=T,nrow=3))
}
exchaxis.3 <- function() {
  # Permutation matrix C3: exchange axes 1,2 reverse 3
  return( matrix( c( 0, 1, 0,
                     1, 0, 0,
                     0, 0,-1),byrow=T,nrow=3))
}
exchange.axis <- function(x,j) {
   if(j==1) return(x%*%exchaxis.1())
   if(j==2) return(x%*%exchaxis.2())
   if(j==3) return(x%*%exchaxis.3())
   return()
}
f77.exchange.axis <- function(x,j) {
   flist <- .Fortran("exchange_axis_", as.double(x), as.integer(j))
   x <- array(flist[[1]], dim=c(3,3))
   return(x)
}   

# Reverse two axes Matrices
invaxis.1 <- function() {
  # Axis Inversion matrix A1 (reverse directions of axes 2 and 3)
  return( matrix( c( 1, 0, 0,
                     0,-1, 0,
                     0, 0,-1),byrow=T,nrow=3))
}
invaxis.2 <- function() {
  # Axis Inversion matrix A2 (reverse directions of axes 1 and 3)
  return( matrix( c(-1, 0, 0,
                     0, 1, 0,
                     0, 0,-1),byrow=T,nrow=3))
}
invaxis.3 <- function() {
  # Axis Inversion matrix A3 (reverse directions of axes 1 and 2)
  return( matrix( c(-1, 0, 0,
                     0,-1, 0,
                     0, 0, 1),byrow=T,nrow=3))
}


# Cycle axes Matrices
cycleaxis.1 <- function() {
  # Cycle Axis matrix T1 (postmultiply: [x,y,z] becomes [x,y,z])(no effect)
  return( matrix( c( 1, 0, 0,
                     0, 1, 0,
                     0, 0, 1),byrow=T,nrow=3))
}
cycleaxis.2 <- function() {
  # Cycle Axis matrix T2 (postmultiply: [x,y,z] becomes [z,x,y])
  return( matrix( c( 0, 1, 0,
                     0, 0, 1,
                     1, 0, 0),byrow=T,nrow=3))
}
cycleaxis.3 <- function() {
  # Cycle Axis matrix T3 (postmultiply: [x,y,z] becomes [y,z,x])
  return( matrix( c( 0, 0, 1,
                     1, 0, 0,
                     0, 1, 0),byrow=T,nrow=3))
}

#-----------------------------------------------------------------------------
# Euler angles and rotation matrices
rotmat.phiv.slow <- function(phiv) {
  # Euler rotation matrix
  rotmat.b3mat(phiv[1]) %*% rotmat.b2mat(phiv[2]) %*% rotmat.b3mat(phiv[3])
}
rotmat.phiv <- function(phiv) {
  # Euler rotation matrix
  retval <- array( rotmat.phiv.asrow(phiv), dim=c(3,3))
  return(retval)
}
f77.rotmat.phiv <- function(phiv) {
  # Euler rotation matrix (f77 version)
  rmat <- array(0,dim=c(3,3))
  flist <- .Fortran("rotmat_phiv_", as.double(phiv), as.double(rmat))
  rmat <- array(flist[[2]],dim=c(3,3))
  return(rmat)
}
rotmat.phiv.asrow <- function(phiv) {
  # Euler rotation matrices stored as rows
  n <- length(phiv)/3
  dim(phiv) <- c(n,3)
  cp <- cos(phiv[,1]); sp <- sin(phiv[,1]) # cos/sin (phi)
  ct <- cos(phiv[,2]); st <- sin(phiv[,2]) # cos/sin (theta)
  cs <- cos(phiv[,3]); ss <- sin(phiv[,3]) # cos/sin (psi)
  retval <- cbind( cp*ct*cs-sp*ss,  sp*ct*cs+cp*ss, -st*cs,
                  -cp*ct*ss-sp*cs, -sp*ct*ss+cp*cs,  st*ss,
                   cp*st,           sp*st,           ct)
  dim(retval) <- c(n,9)
  return(retval)
}

extract.phiv.rotmat <- function(rmat) {
  # Extract Euler angles from a rotation matrix
  phiv <- rep(NA,3)
  phiv[2] <- acos( rmat[3,3] ) # theta (0<theta<pi)
  if(abs(rmat[3,3]==1)) {
    phiv[1] <- atan2(-rmat[1,2], rmat[2,2] ) # phi (-pi<phi<pi)
    phiv[3] <- 0                             # psi
  } else {
    phiv[1] <- atan2( rmat[2,3], rmat[1,3] ) # phi (-pi<phi<pi)
    phiv[3] <- atan2( rmat[3,2],-rmat[3,1] ) # psi (-pi<psi<pi)
  }
  phiv <- phiv%%(2*pi)
  return(phiv)
}
f77.extract.phiv.rotmat <- function(rmat) {
  phiv <- rep(0,3)
  flist <- .Fortran("extract_phiv_rotmat_", as.double(rmat), as.double(phiv))
  phiv <- flist[[2]]%%(2*pi)
  return(phiv)
}
extract.phiv.rotmat.asrow <- function(rmat.asrow) {
  # Extract Euler angles from a rotation matrix
  phiv <- t(apply(rmat.asrow, 1,
                  function(rmat) {
                    extract.phiv.rotmat(array(rmat,dim=c(3,3)))
                  }))
  return(phiv)
}

#-----------------------------------------------------------------------------
# Ambiguity of observations, and alternative representations
flip.rmat <- function(rmat) {
  # swap the x and z columns of matrix rmat, and reverse the y column
  return(cbind(rmat[,3], -rmat[,2], rmat[,1]))
}
flip.phiv <- function(phiv) {
  # swap the x and z axes, and reverse the y axis
  dim0 <- dim(phiv)
  dim(phiv) <- c(length(phiv)/3,3)
  phiv <- t(apply(phiv,1, function(phiv) {
                             rmat <- rotmat.phiv(phiv)
                             rmat <- cbind(rmat[,3],-rmat[,2],rmat[,1])
                             phiv <- extract.phiv.rotmat(rmat)
                             return(phiv)
                           }))
  dim(phiv) <- dim0
  return(phiv)
}
f77.flip.phiv <- function(phiv) {
  dim0 <- dim(phiv)
  n <- length(phiv)/3
  dim(phiv) <- c(n,3)
  phiv.alt <- array(0,dim=c(n,3))
  flist <- .Fortran("flip_phiv_",
                    as.integer(n), 
                    as.double(t(phiv)),
                    as.double(t(phiv.alt)))
  phiv.alt <- t(array(flist[[3]], dim=c(3,n)))
  dim(phiv.alt) <- dim0
  return(phiv.alt)
}

symm.phiv <- function(phiv) {
  # symmetries of Euler angles in a symmetric tensor
  phiv <- rbind(  phiv,                                  #+x, +y, +z
                c(phiv[1]+pi, pi-phiv[2], 2*pi-phiv[3]), #+x, -y, -z
                c(phiv[1]+pi, pi-phiv[2], pi-phiv[3]),   #-x, +y, -z
                c(phiv[1],    phiv[2],    phiv[3]+pi)    #-x, -y, +z
               )
  return(restrict.phiv(phiv))
}
restrict.phiv.theta <- function(phiv) {
  # restrict an (phiv)=(phi,theta,psi)
  # vector so that theta is in the range [0,pi/2]
  phiv <- restrict.phiv(phiv)
  dim0 <- dim(phiv)
  dim(phiv) <- c(length(phiv)/3,3)
  idx <- phiv[,2]>pi/2
  phiv[idx,1] <- phiv[idx,1]+pi
  phiv[idx,2] <- pi-phiv[idx,2]
  phiv[idx,3] <- 2*pi-phiv[idx,3]
  phiv <- restrict.phiv(phiv)
  #idx <- phiv[,3]>pi
  #phiv[idx,3] <- phiv[idx,3] - pi
  dim(phiv) <- dim0
  return(phiv)
}
restrict.phiv <- function(phiv) {
  # put phiv in the ranges ([0,2*pi], [0,pi], [0,2*pi])
  dim0 <- dim(phiv)
  dim(phiv) <- c(length(phiv)/3,3)
  # first put all angles in the range [0,2*pi]
  phiv <- phiv %% (2*pi) 
  # check for any theta values in the range [pi,2*pi]
  idx <- phiv[,2]>pi
  # add pi to phi and psi
  if(length(idx[idx])>0) {
    phiv[idx,1] <- pi+phiv[idx,1]
    phiv[idx,1] <- phiv[idx,1] %% (2*pi)
    phiv[idx,3] <- pi+phiv[idx,3]
    phiv[idx,3] <- phiv[idx,3] %% (2*pi)
    # reverse theta
    phiv[idx,2] <- 2*pi-phiv[idx,2]
  }
  dim(phiv) <- dim0
  return(phiv)
}
restrict.phiv2 <- function(phiv2) {
  # keep theta<pi/2: otherwise reverse the vector
  dim(phiv2) <- c(length(phiv2)/2,2)
  idx <- phiv2[,2]>pi/2
  phiv2[idx,1] <- phiv2[idx,1] + pi
  phiv2[idx,2] <- pi - phiv2[idx,2]
  phiv2 <- phiv2 %% (2*pi)
  if(length(phiv2)==2) dim(phiv2) <- NULL
  return(phiv2)
}

#-----------------------------------------------------------------------------
# random directions
# 2-angles (vector)
rphiv2 <- function(n=1) {
  # generate random directions distributed uniformly on the sphere
  phiv2 <- cbind( runif(n,0,2*pi), acos( runif(n,-1,1) ))
  if(length(phiv2)==2) dim(phiv2) <- NULL
  return(phiv2)
}
nvec.phiv2 <- function(phiv2) {
  dim(phiv2) <- c(length(phiv2)/2,2)
  st <- sin(phiv2[,2]) # sin(theta)
  retval <- cbind( st*cos(phiv2[,1]), st*sin(phiv2[,1]), cos(phiv2[,2]) )
  if(length(retval)==3) dim(retval) <- NULL
  return(retval)
}
ranvec <- function(n=1) {
  # Calculate random unit vector(s) uniformly distributed on the sphere
  phiv2 <- rphiv2(n)   # generate random angles uniformly on the sphere
  retval <- nvec.phiv2(phiv2) # compute the unit vectors for these directions
  return(retval)
}
extract.phiv2.nvec <- function(nvec) {
  # get the azimuth (phi) and colatitude=polar distance (theta)
  # of a unit vector
  nvec <- pmax(-1,pmin(1,as.vector(nvec)))
  dim(nvec) <- c(length(nvec)/3, 3)
  theta <- acos(nvec[,3])
  phi <- atan2( nvec[,2], nvec[,1] )
  phiv2 <- cbind(phi,theta)
  if(length(phiv2)==2) dim(phiv2) <- NULL
  return(phiv2)
}
  
# 3-angles (matrix)
rphiv <- function(n=1) {
  phiv <- cbind( runif(n,0,2*pi), acos( runif(n,-1,1) ), runif(n,0,2*pi) )
  if(length(phiv)==3) dim(phiv) <- NULL
  return(phiv)
}

#-----------------------------------------------------------------------------
# sdr convention
# Conversion between  sdr = strike/dip/rake
#   and Euler angles phiv = phi/theta/psi
convert.sdr.phiv <- function(sdr) {
  dim0 <- dim(sdr)
  dim(sdr) <- c(length(sdr)/3,3)
  phiv <- cbind( sdr[,1]+pi/2, pi-sdr[,2], sdr[,3]-pi/2 )
  phiv <- restrict.phiv(phiv)
  dim(phiv) <- dim0
  return(phiv)
}
f77.convert.sdr.phiv <- function(sdr) {
  dim0 <- dim(sdr)
  n <- length(sdr)/3
  dim(sdr) <- c(n,3)
  phiv <- array(0, dim=c(n,3))
  flist <- .Fortran("convert_sdr_phiv_", as.integer(n), 
                    as.double(t(sdr)), as.double(t(phiv)))
  phiv <- t(array(flist[[3]], dim=c(3,n)))
  dim(phiv) <- dim0
  return(phiv)
}
convert.phiv.sdr <- function(phiv, restrict.dip=F) {
  dim0 <- dim(phiv)
  dim(phiv) <- c(length(phiv)/3,3)
  sdr <- cbind( phiv[,1]-pi/2, pi-phiv[,2], phiv[,3]+pi/2 )
  sdr <- restrict.sdr(sdr)
  if(restrict.dip) sdr <- restrict.sdr.dip(sdr)
  dim(sdr) <- dim0
  return(sdr)
}
f77.convert.phiv.sdr <- function(phiv,restrict.dip=F) {
  dim0 <- dim(phiv)
  n <- length(phiv)/3
  dim(phiv) <- c(n,3)
  sdr <- array(0, dim=c(n,3))
  flist <- .Fortran("convert_phiv_sdr_", as.integer(n), 
                    as.double(t(phiv)), as.double(t(sdr)),
                    as.integer(restrict.dip))
  sdr <- t(array(flist[[3]], dim=c(3,n)))
  dim(sdr) <- dim0
  return(sdr)
}

# rotation matrix in sdr
rotmat.sdr <- function(sdr) {
  phiv <- convert.sdr.phiv(sdr)
  rotmat.phiv(phiv)
}
f77.rotmat.sdr <- function(sdr) {
  rmat <- array(0,dim=c(3,3))
  flist <- .Fortran("rotmat_sdr_", as.double(sdr), as.double(rmat))
  rmat <- array(flist[[2]],dim=c(3,3))
  return(rmat)
}
extract.sdr.rotmat <- function(rmat) {
  phiv <- extract.phiv.rotmat(rmat)
  convert.phiv.sdr(phiv)
}
restrict.sdr.dip <- function(sdr) {
  # restrict an (sdr) vector so that the dip is in the range [0,pi/2]
  sdr <- restrict.sdr(sdr)
  dim0 <- dim(sdr)
  dim(sdr) <- c(length(sdr)/3,3)
  # put everything in the range [0,2pi]
  sdr <- sdr %% (2*pi)
  # see if any dips are too large and fix these
  idx <- sdr[,2]>pi/2
  sdr[idx,] <- cbind(sdr[idx,1]+pi, pi-sdr[idx,2], 2*pi-sdr[idx,3])
  # put everything in the range [0,2pi]
  sdr <- sdr %% (2*pi)
  dim(sdr) <- dim0
  return(sdr)
}
restrict.sdr <- function(sdr) {
  # put sdr in the ranges ([0,2*pi], [0,pi], [0,2*pi])
  dim0 <- dim(sdr)
  dim(sdr) <- c(length(sdr)/3,3)
  sdr <- sdr %% (2*pi)
  idx <- sdr[,2]>pi
  sdr[idx,2] <- 2*pi-sdr[idx,2]
  sdr[idx,1] <- sdr[idx,1] + ifelse(sdr[idx,1]>=pi,-pi,pi)
  sdr[idx,3] <- sdr[idx,3] + ifelse(sdr[idx,3]>=pi,-pi,pi)
  dim(sdr) <- dim0
  return(sdr)
}

rotangle <- function(rmat,zero.tol=1.0e-15) {
  # find the rotation angle and null vector of a rotation matrix
  if(!is.orthogonal(rmat)) return(list(angle=NA,nullvec=c(NA,NA,NA)))
  ee <- eigen(rmat)
  idx <- abs(Im(ee$values))<zero.tol
  if(length(idx[idx])!=1)  return(list(angle=NA,nullvec=c(NA,NA,NA)))
  alpha <- acos((sum(diag(rmat))-1)/2)
  return(list(angle=alpha,nullvec=Re(ee$vectors[,idx])))  
}
#############################################################################
# ambiguity of observations

flip.sdr <- function(sdr) {
  # swap the x and z axes, and reverse the y axis
  # i.e. generate the alternative representation of an sdr observation
  phiv <- convert.sdr.phiv(sdr)
  phiv <- flip.phiv(phiv)
  sdr <- convert.phiv.sdr(phiv, restrict.dip=T)
  return(sdr)
}
f77.flip.sdr <- function(sdr) {
  dim0 <- dim(sdr)
  n <- length(sdr)/3
  dim(sdr) <- c(n,3)
  sdr.alt <- array(0,dim=c(n,3))
  flist <- .Fortran("flip_sdr_",
                    as.integer(n), 
                    as.double(t(sdr)),
                    as.double(t(sdr.alt)))
  sdr.alt <- t(array(flist[[3]], dim=c(3,n)))
  dim(sdr.alt) <- dim0
  return(sdr.alt)
}

#############################################################################
# angle of rotation of an orthogonal matrix

rotangle.rmat <- function(rmat) {
  # get the rotation angle of rotation matrix rmat
  tt <- sum(diag(rmat)) # trace
  t1 <- max(-2, min(2, tt-1))  # -1 <= tt <= 3
  angle <- atan2( sqrt(4-t1^2), t1)
  return(angle)
}
rotangle.rmat.slow <- function(rmat) {
  # get the rotation angle of rotation matrix rmat
  # use eigenvalues to do the calculation
  angle <- max(Arg(eigen(rmat)$values))
  return(angle)
}
f77.rotangle.rmat <- function(rmat) {
  # get the rotation angle of rotation matrix rmat
  rangle <- 0
  flist <- .Fortran("rotangle_rmat_", as.double(rmat), as.double(rangle))
  return(flist[[2]])
}

misfits.rmat <- function(rmat1, rmat2) {
  # the matrix of angles of misfit between two orthogonal matrices
  # NB - this is NOT the same as the single misfit angle between the
  #      two matrices
  acos( t(rmat1)%*%rmat2 )
}
misfits.phiv <- function(phiv1, phiv2) {
  # calculate the misfit angles between
  # parallel sets of angles phiv1 and phiv2
  rmat1.asrow <- rotmat.phiv.asrow(phiv1)
  rmat2.asrow <- rotmat.phiv.asrow(phiv2)
  t(apply( cbind(rmat1.asrow, rmat2.asrow), 1, 
           function(rmat12) {
             rmat1 <- array(rmat12[1:9], dim=c(3,3))
             rmat2 <- array(rmat12[10:18], dim=c(3,3))
             misfits.rmat(rmat1, rmat2)
           }))
}
misfits.sdr <- function(sdr1, sdr2) {
  # calculate the misfit angles between
  # parallel sets of angles sdr1 and sdr2
  phiv1 <- convert.sdr.phiv(sdr1)
  phiv2 <- convert.sdr.phiv(sdr2)
  misfits.phiv(phiv1, phiv2)
}


#############################################################################
# Circular simulation

# (1) Circular vonMises
dcircvm <- function(phi,phi0,kvm,akvm=NULL) {
  # density of circular vonMises distribution
  if(is.null(akvm)) akvm <- 1/(2*pi*besselI(kvm,0))
  exp(kvm*cos(phi-phi0)) * akvm
}
pcircvm <- function(q,phi0,kvm,akvm=NULL,phi.start=0,offset.val=NULL) {
  # cdf of circular vonMises distribution
  if(is.null(akvm)) akvm <- 1/(2*pi*besselI(kvm,0))
  retval <- apply(cbind(q), 1,
        function(q,phi0,kvm,akvm,phi.start) {
           integrate(dcircvm, phi.start, q, phi0=phi0, kvm=kvm, akvm=akvm)$value
        }, phi0=phi0, kvm=kvm, akvm=akvm, phi.start=phi.start)
  if(length(retval)==1) dim(retval) <- NULL
  if(!is.null(offset.val)) retval <- retval-offset.val
  return(retval)
}
qcircvm <- function(p,phi0,kvm,akvm=NULL,phi.start=0) {
  # quantiles of circular vonMises distribution
  if(!is.null(akvm)) akvm <- 1/(2*pi*besselI(kvm,0))
  retval <- apply(cbind(p), 1,
        function(p,phi0,kvm,akvm,phi.start) {
           if(p==0) return(phi.start)
           if(p==1) return(phi.start+2*pi)
           uniroot(pcircvm, phi.start+c(0,2*pi), 
                   phi0=phi0, kvm=kvm, akvm=akvm,
                   phi.start=phi.start, offset.val=p)$root
        }, phi0=phi0, kvm=kvm, akvm=akvm, phi.start=phi.start)
  if(length(retval)==1) dim(retval) <- NULL
  return(retval)
}
rcircvm <- function(n,phi0,kvm,akvm=NULL) {
  # draw from Circular vonMises distribution
  if(!is.null(akvm)) akvm <- 1/(2*pi*besselI(kvm,0))
  p <- runif(n)
  retval <- qcircvm(p,phi0=phi0,kvm=kvm,akvm=akvm,phi.start=0)
  return(retval)
}

# (2) Spherical vonMises = Fisher von Mises
dsphericalvm <- function(phiv2,phiv2.0,kvm,omega=F) {
  # density of spherical vonMises distribution
  # if omega=T then omit the sin(theta)/(4*pi) factor
  nvec.0 <- nvec.phiv2(phiv2.0)
  dim(phiv2) <- c(length(phiv2)/2,2)
  dd <- apply(phiv2, 1,
              function(phiv2,nvec.0) {
                 nvec  <- nvec.phiv2(phiv2)
                 return(dot(nvec,nvec.0))
              }, nvec.0=nvec.0)
  retval <- exp(kvm*dd) * kvm/sinh(kvm)
  if(!omega) retval <- retval * sin(phiv2[,2])/(4*pi)
  return(retval)
}
psphericalvm.0 <- function(q,kvm) {
  # cdf (in theta) of spherical vonMises distribution centered on the pole
  retval <- (exp(kvm)-exp(kvm*cos(q)))/(2*sinh(kvm))
  return(retval)
}
qsphericalvm.0 <- function(p,kvm) {
  # quantiles (in theta) of spherical vonMises distribution centred on the pole
  retval <- acos((log((1-p)*exp(kvm)+p*exp(-kvm)))/kvm)
  return(retval)
}
rsphericalvm <- function(n,phiv2.0,kvm) {
  # draw from Spherical vonMises distribution
  phiv.0 <- c(phiv2.0,0)
  rmat <- rotmat.phiv(phiv.0)
  phi <- runif(n,0,2*pi)
  p <- runif(n)
  theta <- qsphericalvm.0(p,kvm)
  phiv2.1 <- cbind(phi,theta)
  nvecs <- apply(phiv2.1, 1, nvec.phiv2)
  xvecs <- rmat%*%nvecs
  phiv2 <- extract.phiv2.nvec(t(xvecs))
  return(phiv2)
}

# (3) Matrix vonMises
# normalisation constant
f.matvm.nkvm.integrand <- function(u, kvm, exp.n4kvm=NULL) {
  # integrand for normalisation constant
  if(is.null(exp.n4kvm)) exp.n4kvm <- exp(-4*kvm)
  x <- 1 + cos(u)
  ifelse(x==0, 2*kvm*exp.n4kvm, (exp(2*kvm*(x-2))-exp.n4kvm)/x)
}
f.matvm.nkvm <- function(q, kvm, exp.n4kvm=NULL, nkvm=NULL, verbose=F) {
  # integral of f.matvm.nkvm.integrand
  # integrates to normalisation constant for cdf of matrix vonMises distribution
  #if(verbose) cat(sprintf("%10.4f", q))
  if(is.null(exp.n4kvm)) exp.n4kvm <- exp(-4*kvm)
  if(q>pi && is.null(nkvm)) {
     nkvm <- f.matvm.nkvm(pi, kvm, exp.n4kvm=exp.n4kvm)
     retval <- integrate(f.matvm.nkvm.integrand, 0, 2*pi-q,
                         kvm=kvm, exp.n4kvm=exp.n4kvm)$value
     retval <- 2*nkvm - retval
  } else {
     retval <- integrate(f.matvm.nkvm.integrand, 0, q,
                         kvm=kvm, exp.n4kvm=exp.n4kvm)$value
  }
  #if(verbose) cat(sprintf("%10.4f\n", retval))
  return(retval)
}
q.matvm.nkvm <- function(q, kvm, exp.n4kvm=NULL, nkvm=NULL, offsetval=0) {
  # alternative call to f.matvm.nkvm with an offset subtracted
  if(is.null(exp.n4kvm)) exp.n4kvm <- exp(-4*kvm)
  if(is.null(nkvm)) nkvm <- f.matvm.nkvm(pi, kvm, exp.n4kvm=exp.n4kvm)
  apply(cbind(q), 1,
        function(q,offsetval,...) {
           f.matvm.nkvm(q, ...) - offsetval
        }, offsetval=offsetval,
           kvm=kvm, exp.n4kvm=exp.n4kvm, nkvm=nkvm, verbose=F)
}

dmatvm <- function(phiv,kvm,phiv.0=c(0,0,0),
                   exp.n4kvm=NULL,nkvm=NULL,omega=F) {
  # density of matrix vonMises distribution
  # if omega=T then omit the sin(theta)/(8*pi^2) factor
  if(is.null(exp.n4kvm)) exp.n4kvm <- exp(-4*kvm)
  if(is.null(nkvm)) nkvm <- f.matvm.nkvm(pi, kvm, exp.n4kvm=exp.n4kvm)
  dim(phiv) <- c(length(phiv)/3,3)
  rmat.0 <- rotmat.phiv(phiv.0)
  rmat.asrow <- rotmat.phiv.asrow(phiv)
  dd <- apply(rmat.asrow, 1,
              function(rmat,rmat.0,kvm) {
                 exp( kvm*sum(diag(array(rmat,dim=c(3,3))%*%t(rmat.0))) ) 
              }, rmat.0=rmat.0, kvm=kvm)
  # normalising constant
  akvm <- 2*pi*kvm*exp(-3*kvm)/nkvm
  if(!omega) akvm <- sin(phiv[,2])/(8*pi^2)
  return(akvm*dd)
}
rmatvm <- function(n,kvm,phiv.0=c(0,0,0),exp.n4kvm=NULL,nkvm=NULL) {
  # simulate from a matrix vonMises distribution
  if(is.null(exp.n4kvm)) exp.n4kvm <- exp(-4*kvm)
  if(is.null(nkvm)) nkvm <- f.matvm.nkvm(pi, kvm, exp.n4kvm=exp.n4kvm)
  rmat.0 <- rotmat.phiv(phiv.0)
  # first draw psi
  psi <- runif(n,0,2*pi)
  # then draw phi|psi
  # first draw u from (0,pi)
  p <- runif(n)/2 # i.e. only select u in the first half of its range
  u <- apply(cbind(p), 1,
             function(p, kvm, exp.n4kvm, nkvm) {
                uniroot(q.matvm.nkvm, c(0,pi),
                        kvm=kvm, exp.n4kvm=exp.n4kvm, nkvm=nkvm,
                        offsetval=2*p*nkvm)$root
             }, kvm=kvm, exp.n4kvm=exp.n4kvm, nkvm=nkvm)
  # extend to (0,2*pi)
  u <- u*ifelse(runif(n)>0.5, 1, -1)
  # and compute phi
  phi <- u-psi
  # finally draw theta|phi,psi
  p <- runif(n)
  beta <- kvm*(1+cos(phi+psi))
  theta <- acos(1 + (log(1-p + p*exp(-2*beta)))/beta)
  phiv.1 <- cbind(phi,theta,psi)
  # now convert this draw into one from a distribution with mode rmat.0
  rmat.1 <- array( rotmat.phiv.asrow(phiv.1), dim=c(n,3,3))
  rmat.1 <- aperm(rmat.1, c(1,3,2)) # this transposes each rmat.1
  dim(rmat.1) <- c(n*3,3)
  rmat.2 <- rmat.1 %*% rmat.0
  dim(rmat.2) <- c(n,3,3)
  phiv.2 <- t(apply(rmat.2, 1, extract.phiv.rotmat))
  return(phiv.2)
}

add.vm.err.phiv <- function(phiv, kvm) {
  # Add matrix von Mises error to phiv
  dim0 <- dim(phiv)
  dim(phiv) <- c(length(phiv)/3,3)
  phiv.err <- t(apply(phiv, 1,
                      function(phiv, kvm) rmatvm(1, kvm, phiv), kvm=kvm))
  dim(phiv.err) <- dim0
  return(phiv.err)
}
add.vm.err.sdr <- function(sdr, kvm) {
  # Add matrix von Mises error to sdr
  phiv <- convert.sdr.phiv(sdr)
  phiv.err <- add.vm.err.phiv(phiv, kvm)
  sdr.err <- convert.phiv.sdr(phiv.err)
  return(sdr.err)
}

# More von Mises distrbutions
# Alternative method of simulation
dvm <- function(theta,mu=0,kappa=1) {
  # density of vonmises(mu,kappa)
  knorm <- besselI(kappa,0)
  exp(kappa*cos(theta-mu))/(2*pi*knorm)
}
rvm <- function(n,mu=0,kappa=1) {
  # random sample from vonmises(mu,kappa)
  # Use Best+Fisher (1979) method - See MJ p43 Sec 3.5
  a <- 1 + sqrt(1+4*kappa^2)
  b <- (a-sqrt(2*a))/(2*kappa)
  r <- (1+b^2)/(2*b)
  
  theta <- NULL
  while(n>0) {
    u <- array(runif(3*n),dim=c(n,3))
    z <- cos(pi*u[,1])
    f <- (1+r*z)/(r+z)
    c <- kappa*(r-f)
    tval <- mu + sign(u[,3]-0.5)*acos(f)
    accept <- (c*(2-c)-u[,2]>0) | (log(c/u[,2])+1-c>0)
    tval <- tval[accept]
    theta <- c(theta,tval)
    n <- n-length(tval)
  }
  return(theta)
}

# Axial von Mises distribution
davm <- function(theta,mu=0,kappa=1) {
  # density of AxialVonMises(mu,kappa)
  knorm <- besselI(kappa,0)
  0.5*(exp(kappa*cos(theta-mu))+exp(kappa*cos(theta-mu+pi)))/(2*pi*knorm)
}
ravm <- function(n,mu=0,kappa=1) {
  # random sample from AxialVonMises(mu,kappa)
  # Use Best+Fisher (1979) method
  theta <- rvm(n,mu=mu,kappa=kappa)
  theta <- theta + pi*(runif(n)>0.5)
  return(theta)
}

#############################################################################
# Adapative integration - test routines

ftest <- function(x) {
  # test R function
  sin(x[1])*x[2] - x[1]*x[1]
}
ftest.f77 <- function(x) {
  # wrapper for test f77 function
  ndim <- 2; nfunc <- 1; fval <- 0
  flist <- .Fortran("ftest",
                    as.double(ndim),
                    as.double(x),
                    as.double(nfunc),
                    as.double(fval)
                    )
  return(flist[[4]])
}
integrate.ftest <- function(a,b) {
   # wrapper for f77 integrator of f77 function
   ndim <- 2
   minpts <- 100; maxpts <- 10000
   epsabs <- 0.01; epsrel <- 1e-6
   result <- 0; abserr <- 0;
   neval <- 0; ifail <- 0;
   nwmax <- 1
   numfun <- 1
   key <- 1
   if(ndim==1) {
       num <- 15
   } else if(key==0) {
       if(ndim<12) {
             num <- 1 + 2*ndim*(ndim+3) + 2**ndim
       } else {
             num <- 1 + 2*ndim*(ndim+4)
       }
   } else if(key==1) {
       num <- 1 + 2*ndim*(ndim+3) + 2**ndim
   } else if(key==2) {
       num <- 1 + 4*ndim + 6*ndim*ndim + 4*ndim*(ndim-1)*(ndim-2)/3 + 2**ndim
   } else {
       num <- 1 + 2*ndim*(ndim+4)
   }
   maxsub <- (maxpts-num)/(2*num) + 1
   nwmax <- maxsub*( 2*ndim + 2*numfun + 2 ) + 7*numfun + ndim
   work <- rep(0,nwmax)
   flist <- .Fortran("integrate_ftest_",
                     as.integer(ndim),
                     as.double(a), as.double(b),
                     as.integer(minpts), as.integer(maxpts), 
                     as.double(epsabs),  as.double(epsrel), 
                     as.double(result),  as.double(abserr),
                     as.integer(neval),  as.integer(ifail), 
                     as.integer(nwmax),  as.double(work))
   return(flist[[8]])
}
int.ftest <- function(a,b,...) {
  # integrate the function ftest over (a,b) 3 different ways:
  # 1 - use R integrand, R integrator (adapt)
  # 2 - use f77 integrand, R integrator (adapt)
  # 3 - use f77 integrand, f77 integrator (adbays)
  require(adapt)
  ndim <- 2
  r1 <- adapt(ndim,a,b,minpts=100,maxpts=NULL,ftest,eps=0.01,...)
  r2 <- adapt(ndim,a,b,minpts=100,maxpts=NULL,ftest.f77,eps=0.01,...)
  r3 <- integrate.ftest(a,b)
  return(c(r1$value,r2$value,r3))
}

#############################################################################
boxcar <- function(x, w=c(1,3,1), wrap=F,idend=F, f77=T) {
  # boxcar filter the vector (x) by averaging using a set of weights w
  if(f77) {
     n <- length(x)
     nw <- length(w)
     xout <- rep(0,n)
     flist <- .Fortran("boxcar",as.integer(n),as.double(x),
                                as.integer(nw),as.double(w),
                                as.integer(wrap),as.integer(idend),
                                as.double(xout))
     return(flist[[7]])
  } else {
     if(wrap && idend) x <- x[-length(x)]
     is.even <- function(x) x==2*trunc(x/2)
     is.odd <- function(x) !is.even(x)
     n <- length(x)
     nw <- length(w)
     if(nw==1) return(x)              # no smoothing
     if(nw>trunc(n/2)) stop("w is too wide")
     if(is.even(nw)) stop("w must have an odd length")
     if(any(w<0) || sum(w)==0) stop("invalid weights w")
     w <- w/sum(w)
     nw1 <- (nw-1)/2
     if(wrap) {
       x <- c(x[n-(nw1:1)+1],x,x[1:nw1])
     } else {
       x <- c(rep(x[1],nw1),x,rep(x[n],nw1))
     }
     retval <- sapply(1:n, function(i) sum(w*x[nw1+i+(-nw1:nw1)]))
     if(wrap && idend) retval <- c(retval,retval[1])
     return(retval)
  }
}
array.column.boxcar <- function(x,s,wrap=F,idend=F) {
  # Apply a triangular boxcar filter to each column of matrix X(N1,N2)
  # Construct weights using the vector S(N2) - width 2*S+1 in each case
  n1 <- nrow(x)
  n2 <- ncol(x)
  s <- pmin(s,trunc(n1/2))
  if(length(s)!=n2) stop("s incorrect length")
  mmax <- 2*max(s)+1
  w <- rep(0,mmax)
  xtmp <- rep(0,n1)
  flist <- .Fortran("array_column_boxcar",
                    as.integer(n1),as.integer(n2),
                    as.double(x),as.integer(mmax),
                    as.integer(s),
                    as.integer(wrap),as.integer(idend),
                    as.double(w),as.double(xtmp))
  x <- array(flist[[3]],dim=dim(x))
  return(x)
}
#############################################################################

