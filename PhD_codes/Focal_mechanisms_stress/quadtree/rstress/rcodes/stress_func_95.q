#############################################################################
# Functions for earthquake models
#############################################################################

#############################################################################
# Weighting of earthquakes by distance
dist.weight.f <- function(r,rmin=0,rmax=100,pp=1) {
  # Weight earthquakes by distance r
  #   = 1 for r<rmin
  #   = 0 for r>rmax
  #   = power pp for rmin<r<rmax
  u <- ((r-rmin)/(rmax-rmin))
  w <- 1-u^pp
  w[r<=rmin] <- 1
  w[r>=rmax] <- 0
  return(w)
}
#############################################################################
# Earth model
make.earth <- function(fix=NULL) {
  # Earth model
  if(is.null(fix$nu)) {
    nu <- runif(1,0,1)
  } else {
    nu <- fix$nu
  }
  if(is.null(fix$phiv.gs)) {
    phiv.gs <- rphiv(1)
  } else {
    phiv.gs <- fix$phiv.gs
  }
  rmat.gs <- rotmat.phiv(phiv.gs)
  return(list(nu=nu, phiv.gs=phiv.gs, rmat.gs=rmat.gs))
}
make.earthquake <- function(n, earth, nvecs.g=NULL) {
  # generation of sdr of n earthquakes in model earth
  # generate random fault normals in stress coordinates
  # if nvecs is not null then use the specified values (fault normals)
  if(is.null(nvecs.g)) {
    nvecs <- ranvec(n)
    dim(nvecs) <- c(n,3)
  } else {
    n <- length(nvecs.g)/3
    dim(nvecs.g) <- c(n,3)
    nvecs <- nvecs.g %*% earth$rmat.gs
  }
  # calculate Euler angles
  phiv <- t(apply(nvecs, 1,
                function(nvec,earth) {
                   rmat.sf <- wallace.bott.mat.nvec(nvec, earth)
                   rmat.gf <- earth$rmat.gs %*% rmat.sf
                   extract.phiv.rotmat(rmat.gf)
                 }, earth=earth))
  # convert to strike/dip/rake
  sdr <- convert.phiv.sdr(phiv, restrict.dip=T)
  return(sdr)
}
as.epar <- function(earth, restrict=F) {
  # express earth parameters as a vector
  # epar = c(nu, phi, theta, psi) = c(nu, phiv.gs)
  if(restrict) earth <- restrict.earth(earth)
  c(earth$nu, earth$phiv.gs)
}
as.earth <- function(epar) {
  # express earth parameters in a list
  earth <- make.earth(fix=list(nu=epar[1], phiv.gs=epar[2:4]))
  earth <- restrict.earth(earth)
  return(earth)
}
restrict.earth <- function(earth) {
  # if nu gets out of range, or if the angles do
  # then this routine puts them back into sensible ranges
  if(is.na(earth$nu)) return(earth) # don't do anything if it's NA
  if(earth$nu<0) {
    #earth$nu <- 1/( 1-earth$nu )
    #earth$rmat.gs <- earth$rmat.gs %*% cycleaxis.3()
    #earth$phiv.gs <- extract.phiv.rotmat(earth$rmat.gs)
    earth$nu <- 1/( 1-1/earth$nu )
    earth$rmat.gs <- earth$rmat.gs %*% exchaxis.1()
    earth$phiv.gs <- extract.phiv.rotmat(earth$rmat.gs)
  } else if(earth$nu>1) {
    #earth$nu <- 1 - 1/earth$nu
    #earth$rmat.gs <- earth$rmat.gs %*% cycleaxis.2()
    #earth$phiv.gs <- extract.phiv.rotmat(earth$rmat.gs)
    earth$nu <- 1/earth$nu
    earth$rmat.gs <- earth$rmat.gs %*% exchaxis.3()
    earth$phiv.gs <- extract.phiv.rotmat(earth$rmat.gs)
  }
  earth$phiv.gs <- restrict.phiv.theta(earth$phiv.gs)
  return(earth)
}


#############################################################################
# direction of maximum horizontal compressive stress shmax
phi.shmax.f <- function(earth) {
  # azimuth of the direction of maximum horizontal compressive stress
  nu <- earth$nu
  rmat.gs <- earth$rmat.gs
  a <- rmat.gs[1,1]*rmat.gs[2,1] + nu*rmat.gs[1,2]*rmat.gs[2,2]
  b <- rmat.gs[1,1]^2-rmat.gs[2,1]^2 + nu*(rmat.gs[1,2]^2-rmat.gs[2,2]^2)
  phi.shmax <- 0.5*atan2(2*a,b)
  phi.shmax <- phi.shmax + c(0,pi/2)
  kk <- kappa.phi.f(phi.shmax,earth)
  #print(cbind(phi.shmax,kk))
  phi.shmax <- phi.shmax[order(kk)[1]] # choose the one that minimises kappa
  phi.shmax <- phi.shmax %% pi
  if(!is.na(phi.shmax) && phi.shmax<0) phi.shmax <- phi.shmax+pi
  return(phi.shmax)
}

kappa.f <- function(nvec.s,nu) {
  # kappa function for a fault normal nvec.s (in S coordinates)
  dim(nvec.s) <- c(length(nvec.s)/3,3)
  nvec.s[,1]^2 + nu*nvec.s[,2]^2
}
kappa.phi.f <- function(phi,earth) {
  # kappa function for a horizontal fault normal at azimuth phi
  nvec <- cbind(cos(phi),sin(phi),0) # nvec in G coordinates
  dim(nvec) <- c(length(nvec)/3,3)
  nvec <- (nvec)%*%(earth$rmat.gs)   # convert to S coordinates
  return(kappa.f(nvec,earth$nu))
}

#############################################################################
# Wallace Bott matrix
# supply a unit vector (fault normal) and a nu value, and get back a matrix:
# [slip, aux, normal]
wallace.bott.mat.nvec <- function(nvec, earth) {
  # supply normal vector nvec and earth parameters
  # and calculate the wallace.bott matrix for that normal vector
  nu <- earth$nu
  kappa <- nvec[1]^2 + nu*nvec[2]^2
  lambda <- 1/sqrt(kappa*(1-kappa)-nu*(1-nu)*nvec[2]^2  )
  ##!!print(c(nvec,kappa,kappa*(1-kappa)-nu*(1-nu)*nvec[2]^2,lambda))##!!==
  avec <- lambda*c(-nu*nvec[2]*nvec[3],
                       nvec[1]*nvec[3],
               -(1-nu)*nvec[1]*nvec[2] )
  uvec <- lambda*c(  (1-kappa)*nvec[1],
                   -(kappa-nu)*nvec[2],
                        -kappa*nvec[3] )
  retval <- cbind(uvec,avec,nvec)
  return(retval)
}
f77.wallace.bott.mat.nvec <- function(nvec, earth) {
  wbmat <- array(0,dim=c(3,3))
  flist <- .Fortran("wallace_bott_mat_nvec_",
                    as.double(nvec),
                    as.double(earth$nu),
                    as.double(earth$phiv.gs),
                    as.double(earth$rmat.gs),
                    as.double(wbmat)
                    )
  wbmat <- array(flist[[5]], dim=c(3,3))
  return(wbmat)
}  
wallace.bott.mat.phiv2 <- function(phiv2, earth) {
  # supply angles (phi,theta) of the normal vector,
  # and earth parameters
  # and calculate the wallace.bott matrix for that normal vector
  nvec <- nvec.phiv2(phiv2)
  wallace.bott.mat.nvec(nvec, earth=earth)
}

#############################################################################
# Angelier representation of the stress tensor
angelier.stress.tensor <- function(angelier.par) {
  psi   <- angelier.par[1]
  alpha <- angelier.par[2]
  beta  <- angelier.par[3]
  gamma <- angelier.par[4]
  # stress tensor in Angelier geographic coordinates [east,north,up]
  #amat1 <- matrix(c( cos(psi), alpha,           gamma,
  #                   alpha,    cos(psi+2*pi/3), beta,
  #                   gamma,    beta,            cos(psi+4*pi/3)),
  #                byrow=T, nrow=3)
  # stress tensor in standard geographic coordinates [north,east,down]
  amat1 <- matrix(c( cos(psi+2*pi/3), alpha,    -beta,
                     alpha,           cos(psi), -gamma,
                    -beta,           -gamma,     cos(psi+4*pi/3)),
                  byrow=T, nrow=3)
  return(amat1)
}
earth.angelier <- function(angelier.par) {
   # calculate the properties of the model of the earth
   # using Angelier et al. (1982) parameters (psi,alpha,beta,gamma)
   amat <- angelier.stress.tensor(angelier.par)
   e1 <- eigen(amat)
   # sort the eigenvalues into ascending order
   ord <- order(e1$values)
   dd <- e1$values[ord] # the principal components of stress
   dmat <- diag(dd)     # the stress tensor in S coordinates
   rmat <- e1$vectors[,ord] # axes of stress coordinates in G
   if(det(rmat)<0) rmat[,3] <- -rmat[,3]
   # test: compare amat with rmat%*%dmat%*%t(rmat)
   phi <- (dd[2]-dd[1])/(dd[3]-dd[1]) # (mid-min)/(max-min)
   nu <- (dd[3]-dd[2])/(dd[3]-dd[1]) # (max-mid)/(max-min)
   phiv.gs <- extract.phiv.rotmat(rmat)
   earth.a <- as.earth(c(nu,phiv.gs))
   return(earth.a)
}
#############################################################################
# Posterior under Model 1: one observation and no errors

nu.max.f <- function(phi,theta,nullvalue=0) {
  # maximum allowed value of nu
  c <- -tan(phi)*cos(theta)
  c2 <- c^2
  nu <- 1 + 2*c2 - 2*sqrt(c2*(1+c2))
  nu[cos(phi)*cos(theta)>=0] <- nullvalue
  return(nu)
}

post.model1.2d.phig0 <- function(phi.s, theta.s,
                                 sin.theta.mult=T,
                                 maxval=Inf, both=T, flip=F) {
  # (vector function)
  # 2D (phi,theta) posterior under model 1 with phig = c(0,0,0)

  # sin.theta.mult=T include the sin(theta) factor
  
  # both=T means do both possible interpretations of the observation
  # setting maxval to be finite prevents infinities occurring

  retval <- 0
  
  cos.theta.s <- cos(theta.s)
  sin.theta.s <- sin(theta.s)
  cos.phi.s <- cos(phi.s)
  sin.phi.s <- sin(phi.s)

  if(both || !flip) {
     cc2 <- ( sin.phi.s/cos.phi.s * cos.theta.s )^2
     temp <- cos.theta.s/(pi*cos.phi.s^2)
     temp <- abs( temp * (log((1+cc2)/cc2) - 1/(1+cc2)) )
     temp <- temp*( cos.theta.s*cos.phi.s<=0 )
     retval <- retval + temp
  }

  if(both ||  flip) {
     phi.sa <- atan2(-sin.phi.s*sin.theta.s, cos.theta.s)
     theta.sa <- acos( sin.theta.s*cos.phi.s )
     retval <- retval + post.model1.2d.phig0(phi.sa, theta.sa,
                                            sin.theta.mult=F,
                                            maxval=Inf, both=F, flip=F)
     retval <- retval/2.0
  }
  if(sin.theta.mult) retval <- sin.theta.s*retval
  
  if(maxval<Inf) retval <- pmin(retval,maxval)
  return(retval)
}
post.model1.2d.sdr <- function(phi.s, theta.s, sdr,
                               sin.theta.mult=T,
                               maxval=Inf, both=T, flip=F) {
  # (vector function)
  # 2D (phi,theta) posterior under model 1 for focal mechanism 
  
  # both=T means do both possible interpretations of the observation
  # setting maxval to be finite prevents infinities occurring

  n <- length(phi.s)
  phiv.s <- array(cbind(phi.s,theta.s,0),dim=c(n,3))
  phiv.g <- convert.sdr.phiv(sdr)
  
  # calculate the angles (phiv.s0,0) equivalent to (phiv.s,phiv.g)
  rmat.s <- array( rotmat.phiv.asrow(phiv.s), dim=c(n,3,3) )
  rmat.s <- aperm( rmat.s, c(2,1,3) )
  dim(rmat.s) <- c(3,n*3)
  
  rmat.g <- rotmat.phiv(phiv.g)
  
  rmat.s0 <- t(rmat.g) %*% rmat.s
  dim(rmat.s0) <- c(3,n,3)
  rmat.s0 <- aperm( rmat.s0, c(2,1,3) )
  phiv.s0 <- t(apply(rmat.s0, 1, extract.phiv.rotmat))

  retval <- post.model1.2d.phig0(phiv.s0[,1], phiv.s0[,2],
                                 sin.theta.mult=sin.theta.mult,
                                 maxval=maxval, both=both, flip=flip)
  
  return(retval)
}
  
post.model1.2d <- function(phi.s, theta.s, sdr,
                           sin.theta.mult=T,
                           maxval=Inf, both=T, flip=F) {
  # posterior (phi,theta) for the orientation of the direction
  # of maximum prinrcipal stress
  # for a single earthquake observation with no errors
  # evaluated on a set of (phi,theta) values (parallel vectors)
  
  # both=T means do both possible interpretations of the observation
  # setting maxval to be finite prevents infinities occurring

  if(is.list(sdr)) sdr <- unlist(sdr)
  rmat.gf <- rotmat.sdr(sdr)
  retval <- apply(cbind(phi.s, theta.s), 1,
                  post.model1.2d.sub,
                  rmat.gf=rmat.gf,
                  sin.theta.mult=sin.theta.mult,
                  maxval=maxval, both=both, flip=flip)
  return(retval)
}
post.model1.2d.sub <- function(phiv2.s, rmat.gf,
                               sin.theta.mult=T,
                               maxval=Inf, both=T, flip=F) {
  # posterior (phi,theta) for the orientation of the direction
  # of maximum principal stress
  # for a single earthquake observation with no errors
  # evaluated at a single (phi,theta) pair

  # both=T means do both possible interpretations of the observation
  # setting maxval to be finite prevents infinities occurring

  amat0 <- t(rotmat.phiv(c(phiv2.s[1],phiv2.s[2],0))) %*% rmat.gf
  retval <- 0

  if(both) {
    ivals <- 1:2
  } else if(flip) {
    ivals <- 2
  } else {
    ivals <- 1
  }
  cfact <- ifelse(both,1,2)
  for(i in ivals) {
         if(i==1) amat <- amat0
    else if(i==2) amat <- amat0%*%exchaxis.2()

    if(amat[3,1]*amat[3,3]<=0) {
       cc2 <- (amat[3,2]*amat[3,3]/amat[3,1])^2
       tempval <- amat[3,3]/(1-amat[3,2]^2)*((1+cc2)*log(1+1/cc2)-1)
       retval <- retval + abs(tempval)
    }
  }
  retval <- retval*cfact/(2*pi)
  if(sin.theta.mult) retval <- retval * sin(phiv2.s[2])
  retval <- min(retval,maxval)

  return(retval)
}

#############################################################################
# Model 2 - single observation with errors
# Q function 

old.qk.integrand <- function(phict, rmat.s, nu, kappa, k=1,
                         verbose=F, do.plot=F) {
  # phict = cbind(phi,cos(theta))
  dim(phict) <- c(length(phict)/2,2)
  phiv <- cbind(phict[,1], acos(phict[,2]), NA)
  sp <- sin(phiv[,1]) # sin(phi)
  cp <- cos(phiv[,1]) # cos(phi)
  ct <- phict[,2]     # cos(theta)

  # compute phiv[,3] = psi
  ###!!!!== swap signs: phiv[,3] <- atan2( (1-nu)*sp*cp, -(cp^2+nu*sp^2)*ct )
  phiv[,3] <- atan2( -(1-nu)*sp*cp, +(cp^2+nu*sp^2)*ct )
  rmat.f <- rotmat.phiv.asrow(phiv)
  # test
  # testval <- apply(rmat.f, 1, function(h) {
  #                                h <- array(h,c(3,3))
  #                                -h[1,3]*h[3,3]/h[2,2]
  #                             })
  # if(!all(testval>=0)) {
  #   warning("OY!!")
  #   #print(testval)
  # }
  retval <- apply(rmat.f, 1,
                  function(rmat.f, rmat.s, k) {
                     tmat <- rmat.s%*%array(rmat.f,dim=c(3,3))
                     if(k==2) tmat <- exchaxis.2()%*%tmat
                     retval <- sum(diag(tmat))
                     return(retval)
                  }, rmat.s=rmat.s, k=k)
  retval <- exp( kappa*(retval-3) )

  if(verbose) print(cbind(phict,retval))
  if(do.plot) points(phict[,1], phict[,2])
  return(retval)
}


qintegrand <- function(phict, tau, epar, rmat.gs=NULL, flip=F) {
  # integrand for the Q function
  # phict = (phi, cos(theta))
  # epar provides (nu,phiv.s)
  # except that phiv.s is ignored if rmat.gs is given
  nu <- epar[1]
  if(is.null(rmat.gs)) rmat.gs <- rotmat.phiv(epar[2:4])
  phi <- phict[1]
  sin.phi <- sin(phi); cos.phi <- cos(phi)
  cos.theta <- phict[2]
  
  num <- -(1-nu)*sin.phi*cos.phi
  den <- (cos.phi^2 + nu*sin.phi^2) * cos.theta
  psi <- atan2(num,den)
  phiv.f <- c(phi, acos(cos.theta), psi)
  rmat.sf <- rotmat.phiv(phiv.f)

  tmat <- rmat.gs %*% rmat.sf
  if(flip) {
     tr <- tmat[3,1] - tmat[2,2] + tmat[1,3]
  } else {
     tr <- sum(diag(tmat))
  }
  retval <- exp(tau*(tr-3))
  return(retval)
}
qintegrand.f77 <- function(phict, tau, epar, rmat.gs=NULL, flip=F) {
  # f77 wrapper for the integrand for the Q function
  # phict = (phi, cos(theta))
  # epar provides (nu,phiv.s)
  # except that phiv.s is ignored if rmat.gs is given
  nu <- epar[1]
  has.rmat.gs <- !is.null(rmat.gs)
  if(is.null(rmat.gs)) rmat.gs <- array(0,dim=c(3,3))
  ndim <- 2; numfun <- 1
  fval <- 0
  flist <- .Fortran("qintegrand",
                    as.integer(ndim),
                    as.double(phict),
                    as.integer(numfun),
                    as.double(fval),
                    as.integer(flip),
                    as.double(tau),
                    as.double(epar),
                    as.integer(has.rmat.gs),
                    as.double(rmat.gs))
  return(flist[[4]])
}
fnwmax <- function(ndim=2,minpts=100,maxpts=10000) {
   # compute nwmax for the adbays integrator
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
   return(nwmax)
}
integrate.qintegrand <- function(a=c(0,-1),b=c(2*pi,1),
                                 tau, epar, rmat.gs=NULL, flip=F,
                                 epsabs=0.0, epsrel=1.e-6) {
   # wrapper for f77 integrator of f77 function
   if(is.null(rmat.gs)) rmat.gs <- rotmat.phiv(epar[2:4]) 
   has.rmat.gs <- 1
   
   ndim <- 2
   minpts <- 100; maxpts <- 10000
   result <- 0; abserr <- 0;
   neval <- 0; ifail <- 0;
   nwmax <- fnwmax(ndim,minpts,maxpts)
   work <- rep(0,nwmax)
   flist <- .Fortran("integrate_qintegrand_",
                     as.integer(ndim),
                     as.double(a), as.double(b),
                     as.integer(minpts), as.integer(maxpts), 
                     as.double(epsabs),  as.double(epsrel), 
                     as.double(result),  as.double(abserr),
                     as.integer(neval),  as.integer(ifail), 
                     as.integer(nwmax),  as.double(work),
                     as.integer(flip),
                     as.double(tau), as.double(epar),
                     as.integer(has.rmat.gs),
                     as.double(rmat.gs))
    return(flist[[8]])
}
int.qintegrand <- function(a=c(0,-1),b=c(2*pi,1), tau, epar, rmat.gs=NULL,
                           flip=F, eps=1e-6) {
  # integrate the function ftest over (a,b) 3 different ways:
  # 1 - use R integrand, R integrator (adapt)
  # 2 - use f77 integrand, R integrator (adapt)
  # 3 - use f77 integrand, f77 integrator (adbays)
  require(adapt)
  if(is.null(rmat.gs)) rmat.gs <- rotmat.phiv(epar[2:4])
  
  ndim <- 2
  time.start <- Sys.time()
  r1 <- adapt(ndim,a,b,minpts=100,maxpts=NULL,qintegrand,eps=eps,
              tau=tau,epar=epar,rmat.gs=rmat.gs,flip=flip)
  time.now <- Sys.time(); r1t <- time.now-time.start
  time.start <- Sys.time()
  r2 <- adapt(ndim,a,b,minpts=100,maxpts=NULL,qintegrand.f77,eps=eps,
              tau=tau,epar=epar,rmat.gs=rmat.gs,flip=flip)
  time.now <- Sys.time(); r2t <- time.now-time.start
  time.start <- Sys.time()
  r3 <- integrate.qintegrand(a,b,
             tau=tau,epar=epar,rmat.gs=rmat.gs,flip=flip,epsabs=0,epsrel=eps)
  time.now <- Sys.time(); r3t <- time.now-time.start
  return(cbind(c(r1$value,r2$value,r3),c(r1t,r2t,r3t)))
}
qfunc <- function(epar,tau,use.f77=T,epsabs=0,epsrel=1.e-6,flip=F) {
  n <- length(epar)/4
  dim(epar) <- c(n,4)
  fval <- rep(0,n); abserr <- fval; neval <- fval; ifail <- fval
  epabs <- 0; epsrel <- 1.e-6
  nwmax <- fnwmax(); work <- rep(0,nwmax)
  if(use.f77) {
     flist <- .Fortran("qfunc",
                       as.integer(n),
                       as.double( t(epar) ),
                       as.double(tau),
                       as.integer(flip),
                       as.double(fval),
                       as.double(abserr),
                       as.integer(neval),
                       as.integer(ifail),
                       as.double(epsabs),
                       as.double(epsrel),
                       as.integer(nwmax),
                       as.double(work)
                       )
     retval <- flist[[5]]
  } else {
     retval <- apply(epar, 1,
                     function(epar,tau,epsrel,flip) {
                        rmat.gs <- rotmat.phiv(epar[2:4])
                        r1 <- adapt(ndim=2,lower=c(0,-1),upper=c(2*pi,1),
                                    minpts=100,maxpts=NULL,
                                    qintegrand,eps=epsrel,
                                    tau=tau,epar=epar,rmat.gs=rmat.gs)
                        return(r1$value)
                     }, tau=tau, epsrel=epsrel, flip=flip)
  }
  return(retval)
}

#-----------------------------------------------------------------------------
# functions to build an array of posterior values
extend.3dgrid <- function(x) {
  # extend an array x[phi,cos(theta),psi]
  # dimensions nphi = 2*mphi-1
  #            ntheta
  #            npsi
  dd <- dim(x)
  nphi <- dd[1]; ntheta <- dd[2]; npsi <- dd[3]
  nn.phi <- nphi; nn.theta <- 2*ntheta-1; nn.psi <- 2*npsi-1
  mphi <- (nphi+1)/2
  mphi1 <- mphi-1
  idx.phi <- c(mphi1+1:mphi1 , 1:(nphi-mphi1) )
  xnew <- array(NA, dim=c(nn.phi, nn.theta, nn.psi))
  xnew[1:nphi, 1:ntheta, 1:npsi] <- x              #0
  xnew[idx.phi, nn.theta:ntheta, nn.psi:npsi] <- x #1
  xnew[idx.phi, nn.theta:ntheta, npsi:1] <- x      #2
  xnew[1:nphi, 1:ntheta, npsi-1+(1:npsi)] <- x     #3
  xnew[nn.phi,,] <- xnew[1,,]
  return(xnew)
}
extend.4dgrid <- function(x) {
  # extend an array x[nu,phi,cos(theta),psi]
  # dimensions nphi = 2*mphi-1
  #            ntheta
  #            npsi
  dd <- dim(x)
  nnu <- dd[1]; nphi <- dd[2]; ntheta <- dd[3]; npsi <- dd[4]
  nn.nu <- nnu; nn.phi <- nphi; nn.theta <- 2*ntheta-1; nn.psi <- 2*npsi-1
  mphi <- (nphi+1)/2
  mphi1 <- mphi-1
  idx.phi <- c(mphi1+1:mphi1 , 1:(nphi-mphi1) )
  
  # ensure symmetry ##!!==
  x[,nphi,,] <- x[,1,,] 
  x[,,,npsi] <- x[,,,1]
  x[,idx.phi,ntheta,npsi:1] <- x[,,ntheta,]
  
  # now build the expanded array
  xnew <- array(NA, dim=c(nn.nu, nn.phi, nn.theta, nn.psi))
  xnew[, 1:nphi, 1:ntheta, 1:npsi] <- x              #0
  xnew[, idx.phi, nn.theta:ntheta, nn.psi:npsi] <- x #1
  xnew[, idx.phi, nn.theta:ntheta, npsi:1] <- x      #2
  xnew[, 1:nphi, 1:ntheta, npsi-1+(1:npsi)] <- x     #3
  xnew[,nn.phi,,] <- xnew[,1,,]
  return(xnew)
}

index.extend.grid <- function(ix, nn) {
  # index transformations to extend an array x[phi,cos(theta),psi]
  # dimensions nphi = 2*mphi-1
  #            ntheta
  #            npsi
  # nn = c(nphi,ntheta,npsi)
  nphi <- nn[1]; ntheta <- nn[2]; npsi <- nn[3]
  mphi <- (nphi+1)/2
  mphi1 <- mphi-1
  nn.phi <- nphi; nn.theta <- 2*ntheta-1; nn.psi <- 2*npsi-1
  # which of the four possible places are we in?
  iphi <- ix[1]; itheta <- ix[2]; ipsi <- ix[3]
  if(ipsi<=npsi) {
    if(itheta<=ntheta) {
      jphi <- iphi; jtheta <- itheta; jpsi <- ipsi                        # 0
    } else {
      jphi <- iphi+mphi1; jtheta <- 2*ntheta-itheta; jpsi <- npsi-ipsi+1  # 2
    }
  } else {
    if(itheta<=ntheta) {
      jphi <- iphi; jtheta <- itheta; jpsi <- ipsi-npsi+1                 # 3
    } else {
      jphi <- iphi+mphi1; jtheta <- 2*ntheta-itheta; jpsi <- 2*npsi-ipsi  # 1
    }
  }
  if(jphi>nn.phi) jphi <- jphi-nn.phi
  jx <- c(jphi,jtheta,jpsi)
  return(jx)
}
retabulate.qflist <- function(phiv.g=c(0,0,0),qflist,flip=F,
                              exchange=0,qv0.only=T) {
  # create a retabulation of an array of qvals
  # appropriate for the observation phiv.g
  # where qflist contains the null tabulation for phiv.g=c(0,0,0)
  # use exchange to swap axes: 0(no change),1(swap 2-3),2(swap 1-3)
  if(is.null(qflist$qvals1) || is.null(qflist$qvals2)) qv0.only <- T
  if(qv0.only) {
     qvals0 <- qflist$qvals0
     nn <- dim(qvals0)
  } else {
     qvals0 <- qflist$qvals0
     qvals1 <- qflist$qvals1
     qvals2 <- qflist$qvals2
     nn <- dim(qvals1)
  }

  nnu <- nn[1]; nphi <- nn[2]; ntheta <- nn[3]; npsi <- nn[4]
  numin <- qflist$epar.vecs$nu[1]
  phimin <- qflist$epar.vecs$phi[1]
  thetamin <- qflist$epar.vecs$theta[1]
  psimin <- qflist$epar.vecs$psi[1]
  dnu <- (qflist$epar.vecs$nu[nnu]-numin)/(nnu-1)
  dphi <- (qflist$epar.vecs$phi[nphi]-phimin)/(nphi-1)
  dtheta <- (qflist$epar.vecs$theta[ntheta]-thetamin)/(ntheta-1)
  dpsi <- (qflist$epar.vecs$psi[npsi]-psimin)/(npsi-1)

  if(any(is.na(qvals0))) {
    qflist$qvals0 <- array(NA,dim=nn)
  } else {
    qvals.new <- array(0,dim=nn)
    flist <- .Fortran("retabulate_qfunc",
                      as.double(phiv.g),
                      as.integer(flip), as.integer(exchange),
                      as.double(qvals.new), as.double(qvals0),
                      as.integer(nnu), as.double(numin), as.double(dnu),
                      as.integer(nphi), as.double(phimin), as.double(dphi),
                      as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                      as.integer(npsi), as.double(psimin), as.double(dpsi))
    qflist$qvals0 <- array(flist[[4]],dim=nn)
  }
  if(!qv0.only) {
     qvals.new <- array(0,dim=nn)
     if(any(is.na(qvals1))) {
       qflist$qvals1 <- array(NA,dim=nn)
     } else {
       flist <- .Fortran("retabulate_qfunc",
                         as.double(phiv.g),
                         as.integer(flip), as.integer(exchange),
                         as.double(qvals.new), as.double(qvals1),
                         as.integer(nnu), as.double(numin), as.double(dnu),
                         as.integer(nphi), as.double(phimin), as.double(dphi),
                         as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                         as.integer(npsi), as.double(psimin), as.double(dpsi))
       qflist$qvals1 <- array(flist[[4]],dim=nn)
     }
     if(any(is.na(qvals2))) {
       qflist$qvals2 <- array(NA,dim=nn)
     } else {
       flist <- .Fortran("retabulate_qfunc",
                         as.double(phiv.g), 
                         as.integer(flip), as.integer(exchange),
                         as.double(qvals.new),as.double(qvals2),
                         as.integer(nnu), as.double(numin), as.double(dnu),
                         as.integer(nphi), as.double(phimin), as.double(dphi),
                         as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                         as.integer(npsi), as.double(psimin), as.double(dpsi))
       qflist$qvals2 <- array(flist[[4]],dim=nn)
     }
  }
  return(qflist)
}

retabulate.qfphiv.list <- function(phiv.g=c(0,0,0),qfphiv.list,flip=F,
                                   exchange=0,qv0.only=T) {
  # create a retabulation of an array of qphiv.vals
  # (same as qfvals, but integrated over nu)
  # appropriate for the observation phiv.g
  # where qflist contains the null tabulation for phiv.g=c(0,0,0)
  if(is.null(qfphiv.list$qphivvals1) || is.null(qfphiv.list$qphivvals2)) qv0.only <- T
  if(qv0.only) {
     qphivvals0 <- qfphiv.list$qphivvals0
     nn <- dim(qphivvals0)
  } else {
     qphivvals0 <- qfphiv.list$qphivvals0
     qphivvals1 <- qfphiv.list$qphivvals1
     qphivvals2 <- qfphiv.list$qphivvals2
     nn <- dim(qphivvals1)
  }

  nphi <- nn[1]; ntheta <- nn[2]; npsi <- nn[3]
  phimin <- qfphiv.list$epar.vecs$phi[1]
  thetamin <- qfphiv.list$epar.vecs$theta[1]
  psimin <- qfphiv.list$epar.vecs$psi[1]
  dphi <- (qfphiv.list$epar.vecs$phi[nphi]-phimin)/(nphi-1)
  dtheta <- (qfphiv.list$epar.vecs$theta[ntheta]-thetamin)/(ntheta-1)
  dpsi <- (qfphiv.list$epar.vecs$psi[npsi]-psimin)/(npsi-1)

  if(any(is.na(qphivvals0))) {
    qfphiv.list$qphivvals0 <- array(NA,dim=nn)
  } else {
    qphivvals.new <- array(0,dim=nn)
    flist <- .Fortran("retabulate_qphiv_func",
                    as.double(phiv.g),
                    as.integer(flip), as.integer(exchange),
                    as.double(qphivvals.new), as.double(qphivvals0),
                    as.integer(nphi), as.double(phimin), as.double(dphi),
                    as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                    as.integer(npsi), as.double(psimin), as.double(dpsi))
    qfphiv.list$qphivvals0 <- array(flist[[4]],dim=nn)
  }
  if(!qv0.only) {
     if(any(is.na(qphivvals1))) {
       qfphiv.list$qphivvals1 <- array(NA,dim=nn)
     } else {
       qphivvals.new <- array(0,dim=nn)
       flist <- .Fortran("retabulate_qphiv_func",
                         as.double(phiv.g),
                         as.integer(flip), as.integer(exchange),
                         as.double(qphivvals.new), as.double(qphivvals1),
                         as.integer(nphi), as.double(phimin), as.double(dphi),
                         as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                         as.integer(npsi), as.double(psimin), as.double(dpsi))
       qfphiv.list$qphivvals1 <- array(flist[[4]],dim=nn)
     }
     if(any(is.na(qphivvals2))) {
       qfphiv.list$qphivvals2 <- array(NA,dim=nn)
     } else {
       flist <- .Fortran("retabulate_qphiv_func",
                         as.double(phiv.g), 
                         as.integer(flip), as.integer(exchange),
                         as.double(qphivvals.new),as.double(qphivvals2),
                         as.integer(nphi), as.double(phimin), as.double(dphi),
                         as.integer(ntheta), as.double(thetamin), as.double(dtheta),
                         as.integer(npsi), as.double(psimin), as.double(dpsi))
       qfphiv.list$qphivvals2 <- array(flist[[4]],dim=nn)
     }
  }
  return(qfphiv.list)
}

invert.tau <- function(sdr, weights=1, tau.qfunc=tau5.qfunc, verbose=T) {
  # invert a set of earthquakes - supplying the appropriate tabulation
  # of the posterior as tau.qfunc
  # options are tau.qfunc = tau5.qfunc (low precision)
  #                      or tau20.qfunc (high precision)
  n <- nrow(sdr)
  if(length(weights)==1) weights <- rep(weights,n)
  for(i in 1:nrow(sdr)) {
     if(verbose) cat(paste("doing",i,"..."))
  
     sdr1 <- sdr[i,]
     wgt1 <- weights[i]
     if(verbose) cat(sprintf("weight=%10.8f;",wgt1))
     phiv.g <- convert.sdr.phiv(sdr1)
     # retabulate for qfunc for this observation
     tau.qfunc.alt <- retabulate.qflist(phiv.g,tau.qfunc,qv0.only=F)
  
     if(i==1) {
       tau.post <- tau.qfunc.alt
       tau.post$qvals0 <- tau.post$qvals0^wgt1
       tau.post$qvals1 <- tau.post$qvals1^wgt1
       tau.post$qvals2 <- tau.post$qvals2^wgt1
     } else {
       tau.post$qvals0 <- tau.post$qvals0 * tau.qfunc.alt$qvals0^wgt1
       tau.post$qvals1 <- tau.post$qvals1 * tau.qfunc.alt$qvals1^wgt1
       tau.post$qvals2 <- tau.post$qvals2 * tau.qfunc.alt$qvals2^wgt1
     }
     if(verbose) cat("done\n")
  }
  return(tau.post)
}
plot.tau.post <- function(tau.post, epar=NULL, func=0,
                          pch=16, col=NULL, cex=1,
                          mfrow=c(2,3), degrees=F) {
  # plot a posterior distribution returned by invert.tau()
  # as 6 2D contour plots
  spar <- par(no.readonly=T); on.exit(par(spar))
  par(mfrow=mfrow)
  mar <- par()$mar
  mar[3:4] <- c(1,1)
  par(mar=mar)
  epar.vecs <- tau.post$epar.vecs
  
  sin.thetavec <- sin(epar.vecs$theta) # multiply by this
  if(degrees) {
    for(i in 2:4) epar.vecs[[i]] <- 180/pi*epar.vecs[[i]]
    if(!is.null(epar)) {
      epar <- array(epar,dim=c(length(epar)/4,4))
      epar[,2:4] <- epar[,2:4]*180/pi
    }
  }
  labs <- list(expression(nu),expression(phi),
               expression(theta),expression(psi))

  for(i in 1:3) {
  for(j in (i+1):4) {
     if(func==1) {
        ff <- tau.post$qvals1
     } else if(func==2) {
        ff <- tau.post$qvals2
     } else {
        ff <- tau.post$qvals0
     }
     aa <- aperm(array(sin.thetavec, dim=dim(ff)[c(3,1,2,4)]),c(2,3,1,4))
     ff <- ff*aa
     zz <- sum.mat.trap(ff,c(i,j))
     contour(epar.vecs[[i]],
             epar.vecs[[j]],zz,
             xlab=labs[[i]],ylab=labs[[j]],
             cex.lab=1.5)
     if(!is.null(epar)) {
       dim(epar) <- c(length(epar)/4,4)
       if(is.null(col)) col <- c("red","blue","green")
       points(epar[,i],epar[,j],pch=pch,cex=cex,col=col)
     }
  }}

  invisible()
}
contour.tau.post.stereogram <- function(tau.post,
           axes=c(T,T,T),
           cols=c("orange","green","blue"),
           labs=c(expression(S[min]),expression(S[mid]),expression(S[max])),
           func=0, add=F, ddz=NULL) {
  # plot a posterior distribution returned by invert.tau()
  # on a stereogram
  # if axes[i] is TRUE then plot that axis
  # 1=min 2=mid 3=max
  # if func=0 then plot q0, if ==1 then q1 and if ==2 then q2
  if(!add) start.stereogram(plot.centre=F)
  nnu <- tau.post$nepar[1]
  nphi <- tau.post$nepar[2]
  ntheta <- tau.post$nepar[3]
  npsi <- tau.post$nepar[4]
  nuvec <- tau.post$epar.vecs$nu
  phivec <- tau.post$epar.vecs$phi
  thetavec <- tau.post$epar.vecs$theta
  psivec <- tau.post$epar.vecs$psi

  if(is.null(ddz)) {
     tt.post <- list(nepar=tau.post$nepar,epar.vecs=tau.post$epar.vecs)
     if(!is.null(tau.post$qvals0)) {
       tt.post$qphivvals0 <- sum.mat.trap(tau.post$qvals0,c(2,3,4))
     }
     if(!is.null(tau.post$qvals1)) {
       tt.post$qphivvals1 <- sum.mat.trap(tau.post$qvals1,c(2,3,4))
     }
     if(!is.null(tau.post$qvals2)) {
       tt.post$qphivvals2 <- sum.mat.trap(tau.post$qvals2,c(2,3,4))
     }
     
     dd <- list(NULL,NULL,NULL)
     if(axes[1]) dd[[1]] <- retabulate.qfphiv.list(c(0,0,0),tt.post,exchange=2)
     if(axes[2]) dd[[2]] <- retabulate.qfphiv.list(c(0,0,0),tt.post,exchange=1)
     if(axes[3]) dd[[3]] <- tt.post
     
     # density of the z-axis
     ddz <- list(NULL,NULL,NULL)
     if(func==1) {
        for(i in (1:3)[axes]) ddz[[i]] <- sum.mat.trap(dd[[i]]$qphivvals1,c(1,2))
     } else if(func==2) {
        for(i in (1:3)[axes]) ddz[[i]] <- sum.mat.trap(dd[[i]]$qphivvals2,c(1,2))
     } else {
        for(i in (1:3)[axes]) ddz[[i]] <- sum.mat.trap(dd[[i]]$qphivvals0,c(1,2))
     }
     for(i in (1:3)[axes]) {
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
        ddz[[i]] <- array.column.boxcar(ddz[[i]],s,wrap=T,idend=T)
     }
  }
     
  # find the posterior mode
  pp <- list(NULL,NULL,NULL)
  for(i in (1:3)[axes]) {
     pp[[i]] <- unlist(
        expand.grid(phivec,thetavec)[as.vector(ddz[[i]])==max(ddz[[i]]),])
     # flip so that both directions of the axes are identified
     pp[[i]] <- rbind(pp[[i]], c(pp[[i]][1]+pi,pi-pp[[i]][2]))
  }

  # make the plots
  for(i in (1:3)[axes]) {
     contour.stereogram(phivec,thetavec,ddz[[i]],add=T,col=cols[i])
  }
  for(i in (1:3)[axes]) {
     points.phiv2.stereogram(pp[[i]],pch=21,cex=2,bg=cols[i])
     text.phiv2.stereogram(pp[[i]],label=labs[i],cex=2.5,pos=3)
  }
  invisible()
}
post.mean.epar.f <- function(tau.post) { 
   # calculate posterior mean of epar
   nepar <- tau.post$nepar
   epar.vecs <- tau.post$epar.vecs
   postpdf <- tau.post$qvals0
   
   postpdf <- postpdf/sum(postpdf)
   wnu <- sum.mat.trap(postpdf,1)
   wphiv <- sum.mat.trap(postpdf,c(2,3,4))

   # posterior mean
   vec.grid <- expand.grid(epar.vecs$phi,epar.vecs$theta,epar.vecs$psi)
   bar.epar <- c( sum(wnu*epar.vecs$nu),
                  mean.phiv.axis(as.matrix(vec.grid),w=as.vector(wphiv)) )
   return(list(bar.epar=bar.epar))
}
tau.posterior.summary <- function(tau.post,tabdir="tab",
                                  phi.shmax.grid=NULL,
                                  verbose=F) {
   # compute posterior summaries
  
   # first check whether the grid for shmax exists - if not then create it
   # this assumes we're using the standard grid
   have.grid <- FALSE
   if(is.null(phi.shmax.grid)) {
     if(!is.null(tabdir)) {
        ifile <- paste(tabdir,"phi.shmax.grid.Rbin",sep="/")
        if(file.exists(ifile)) {
          if(verbose) cat("Reading phi(shmax) grid from file\n")
          load(ifile)
          if(all(phi.shmax.grid.list$nepar==tau.post$nepar)) {
            have.grid <- TRUE
            phi.shmax.grid <- phi.shmax.grid.list$phi.shmax.grid
          } else {
            if(verbose) cat("Grid does not match supplied result object\n")
          }
        }
     }
   } else {
     have.grid <- TRUE
   }
   if(!have.grid) {
      if(verbose) cat("Calculating phi(shmax) grid... (this is slow) ...")
      # compute phi(shmax) on the epar grid
      phi.shmax.grid <- apply(expand.grid(tau.post$epar.vecs),1,
                           function(epar) phi.shmax.f(as.earth(epar)))
      phi.shmax.grid <- array(phi.shmax.grid,dim=tau.post$nepar)
      if(verbose) cat("Done.\n")
   }
 
   # get the relevant parts of tau.post
   nepar <- tau.post$nepar
   epar.vecs <- tau.post$epar.vecs
   qpdf <- tau.post$qvals0
   
   # scale qvals0 by sin(theta) and normalise to give the posterior density
   sin.theta <- sin(epar.vecs$theta)
   a.sin.theta <- array( rep(rep(sin.theta,
                times=rep(nepar[1]*nepar[2],nepar[3])),nepar[4]), dim=nepar)
   postpdf <- qpdf*a.sin.theta
   postpdf <- postpdf/sum(postpdf)
   
   # marginals in nu and phiv
   wnu <- sum.mat.trap(postpdf,1)
   wphiv <- sum.mat.trap(postpdf,2:4)
   
   # grids of angles
   phiv.grid <- as.matrix(expand.grid(epar.vecs[2:4]))
   rmat.grid <- rotmat.phiv.asrow(phiv.grid)
   sin.azphi.grid <- rmat.grid[,c(2,5,8)] # second row of rmat: rmat[2,]
   cos.azphi.grid <- rmat.grid[,c(1,4,7)] # first row of rmat:  rmat[1,]
   azphi.grid <- atan2( sin.azphi.grid, cos.azphi.grid ) # azimuth of each column of rmat
   
   # locate the MAP estimate
   map.idx <- locate.in.array(postpdf==max(postpdf))
   map.epar <- apply(cbind(1:4,map.idx[1,]),1,
                     function(i.idx) tau.post$epar.vecs[[i.idx[1]]][i.idx[2]])
   
   # locate the posterior mean
   vec.grid <- expand.grid(epar.vecs$phi,epar.vecs$theta,epar.vecs$psi)
   if(verbose) cat("Computing posterior mean...")
   mean.epar <- c( sum(wnu*epar.vecs$nu)/sum(wnu),
                   mean.phiv.axis(as.matrix(vec.grid),w=as.vector(wphiv)) )
   if(verbose) cat("Done\n")
   
   # (1) nu
   # posterior distribution of nu
   ff <- cumsum(wnu)/sum(wnu)
   qq.nu <- getq.csum(epar.vecs$nu,ff)
   nu.density <- wnu
   mean.nu <- sum(wnu*epar.vecs$nu)/sum(wnu)
   map.nu <- map.epar[1]
   
   # (2) Phi_S
   # look at map.epar and mean.epar
   map.phiv <- map.epar[2:4]
   mean.phiv <- mean.epar[2:4]

   # integrate over nu
   tt.post <- list(nepar=tau.post$nepar,epar.vecs=tau.post$epar.vecs)
   if(!is.null(tau.post$qvals0)) {
     tt.post$qphivvals0 <- sum.mat.trap(tau.post$qvals0,c(2,3,4))
   }
   if(!is.null(tau.post$qvals1)) {
     tt.post$qphivvals1 <- sum.mat.trap(tau.post$qvals1,c(2,3,4))
   }
   if(!is.null(tau.post$qvals2)) {
     tt.post$qphivvals2 <- sum.mat.trap(tau.post$qvals2,c(2,3,4))
   }
   
   dd <- list(NULL,NULL,NULL)
   dd[[1]] <- retabulate.qfphiv.list(c(0,0,0),tt.post,exchange=2) # exch. 1-3
   dd[[2]] <- retabulate.qfphiv.list(c(0,0,0),tt.post,exchange=1) # exch. 2-3
   dd[[3]] <- tt.post                                             # no exch.

   # density of the z-axis
   ddz <- list(NULL,NULL,NULL)
   for(i in 1:3) {
     ddz[[i]] <- sum.mat.trap(dd[[i]]$qphivvals0,c(1,2))
     # smooth each function in phi
     # grid size in phi and theta
     nphi <- dd[[i]]$nepar[2]
     ntheta <- dd[[i]]$nepar[3]
     # grid spacing in phi
     dphi <- dd[[i]]$epar.vecs$phi[2]-dd[[i]]$epar.vecs$phi[1]
     # for each row scale according to the scale phi.smooth.scale
     phi.smooth.scale <- dphi/2 #= half of grid spacing at the equator
     s <- trunc( phi.smooth.scale/dphi/sin(dd[[i]]$epar.vecs$theta) )
     s <- pmin(s,(nphi/2)) # s=half-width of triangular boxcar filter
     # for the first and last columns take the mean to enforce continuity
     ddz[[i]][,1] <- mean(ddz[[i]][,1])
     ddz[[i]][,ntheta] <- mean(ddz[[i]][,ntheta])
     # now smooth each column
     ddz[[i]] <- array.column.boxcar(ddz[[i]],s,wrap=T,idend=T)
   }
   
   # (3) the azimuths of the three principal axes

   ## old code: computes mean of azimuths - rather than the azimuth of the
   ## mean axis direction
   #rmat.gs.map <- rotmat.phiv(map.epar[2:4])
   #azphi.map <- atan2( rmat.gs.map[2,], rmat.gs.map[1,] )
   #angle.offset <- rep(azphi.map,each=nrow(azphi.grid))+pi/2
   #azphi.density <- apply( (azphi.grid-angle.offset)%%pi+angle.offset,
   #                        2, density, weights=wphiv)
   #azphi.density <- lapply(azphi.density,
   #                        function(dd) {
   #                           dd$y <- pmax(0,dd$y)
   #                           return(dd)
   #                         })
   #qq.azphi <- sapply(azphi.density,
   #                   function(dd) {
   #                     w <- dd$y
   #                     ff <- cumsum(w)/sum(w)
   #                     qq <- getq.csum(dd$x,ff)
   #                     return(qq)
   #                   })
   #mean.azphi <- sapply( azphi.density,
   #                      function(dd) {
   #                        sum(dd$x*dd$y)/sum(dd$y)
   #                      })

   # new code: compute the azimuth of the mean in each case
   # use the retabulations
   azphi.vals <- tau.post$epar.vecs$phi
   nazphi <- length(azphi.vals)
   nstore <- (nazphi+1)/2
   sin.theta.vec <- sin(tau.post$epar.vecs$theta)
   w.azphi <- lapply(ddz, function(ddzz) {
                             w <- apply(ddzz,1,function(x) sum(x*sin.theta.vec))
                             w <- w/sum(w)
                             return(w)
                          })
   azphi.density <- lapply(w.azphi,
                           function(w,x) {
                             w <- w[1:nstore]
                             x <- x[1:nstore]

                             w <- w[-nstore]
                             x <- x[-nstore]
                             xx <- c(x-pi,x,x+pi)
                             ww <- c(w,w,w)
                             if(any(is.na(ww)) || sum(ww)==0) {
                               den <- list(x=rep(NA,2*nstore+1),
                                           y=rep(NA,2*nstore+1))
                             } else {
                               ww <- ww/sum(ww)
                               bw <- x[2]-x[1]
                               den <- density(xx, bw=bw, weights=ww,
                                              from=0, to=pi, n=2*nstore+1)
                             }
                             #y <- boxcar(x=w,w=c(1,3,1),wrap=T,idend=T)
                             #y <- y/sum(y)
                             #den <- list(x=x,y=y,bw=2,n=nstore,
                             #            call=density,
                             #            data.name="w",has.na=F)
                             #attr(den,"class") <- "density"
                             return(den)
                           }, x=azphi.vals)
   qq.azphi <- sapply(azphi.density,
                      function(den) {
                        # locate the maximum - and centre around it
                        n <- length(den$x)
                        x <- den$x[-n]
                        y <- den$y[-n] # remove last values before wrapping
                        n <- n-1
                        imax <- which.max(den$y)
                        nh <- trunc(n/2)
                        ii <- 1+(imax + nh + 1)%%n
                        x <- x - pi + x[ii]
                        w <- y[1+(imax + nh + 1:n)%%n]
                        # compute mean
                        xmean <- sum(w*x)/sum(w)
                        x <- c(x,x[1]+pi)
                        w <- c(w,w[1])
                        ff <- cumsum(w)/sum(w)
                        qq <- c(getq.csum(x,ff),xmean)
                        return(qq)
                      })
   mean.azphi <- (qq.azphi[nrow(qq.azphi),] %% pi) # store the mean
   qq.azphi <- qq.azphi[-nrow(qq.azphi),] # remove mean from the quantile list
   map.azphi <- sapply(azphi.density,
                       function(den) {
                         n <- length(den$x)
                         idx <- min((1:n)[den$y==max(den$y)])
                         return(den$x[idx])
                       })
   
   
   # (4) SHmax
   map.phi.shmax <- phi.shmax.grid[map.idx][1]
   offset <- map.phi.shmax-pi/2 + pi*c(-1,0,1)
   nshmax <- 512
   if(any(is.na(postpdf))) {
     dd <- list(x=rep(NA,nshmax), y=rep(NA,nshmax))
   } else {
     dd <- density(c( (phi.shmax.grid-offset[1])%%pi+offset[1],
                      (phi.shmax.grid-offset[2])%%pi+offset[2],
                      (phi.shmax.grid-offset[3])%%pi+offset[3]),
                   weights=rep(postpdf/3,3), from=0,to=pi,n=nshmax)
   }
   x <- dd$x
   w <- pmax(0,dd$y)
   qq <- 0.5*circular.wquantile(2*x[-1],weights=w[-1],
                            probs=c(0,0.025,0.05,0.1,0.5,0.9,0.95,0.975,1))
   offset <- qq["50%"]-pi/2
   x1 <- (x-offset)%%pi
   mean.phi.shmax <- offset + sum(w[-1]*x1[-1])/sum(w[-1])
   phi.shmax.density <- dd
   ##qq.phi.shmax <- (qq+pi/2)%%pi - pi/2
   qq.phi.shmax <- (qq-offset)%%pi + offset

   retval <- list(map.idx=map.idx,
                  map.epar=map.epar,
                  mean.epar=mean.epar,
                  nu.density=nu.density,
                  qq.nu=qq.nu,
                  mean.nu=mean.nu,
                  map.nu=map.nu,
                  mean.phiv=mean.phiv,
                  map.phiv=map.phiv,
                  ddz=ddz,
                  azphi.density=azphi.density,
                  qq.azphi=qq.azphi,
                  mean.azphi=mean.azphi,
                  map.azphi=map.azphi,
                  phi.shmax.density=phi.shmax.density,
                  qq.phi.shmax=qq.phi.shmax,
                  mean.phi.shmax=mean.phi.shmax,
                  map.phi.shmax=map.phi.shmax)
   return(retval)
}
tau.posterior.plots <- function(tau.post,tau.post.summary=NULL,
                                plot.list=c(
                                  "nu.density",
                                  "phiv",
                                  "tau.post",
                                  "contour","radii",
                                  "azphi",
                                  "phi.shmax.density")
                                ) {
   # make posterior plots
   spar <- par(no.readonly=T)
  
   # get the relevant parts of tau.post
   nepar <- tau.post$nepar
   epar.vecs <- tau.post$epar.vecs
   # get the relevant parts of tau.post.summary
   if(!is.null(tau.post.summary)) {
      map.idx <- tau.post.summary$map.idx
      map.epar <- tau.post.summary$map.epar
      mean.epar <- tau.post.summary$mean.epar
      nu.density <- tau.post.summary$nu.density
      qq.nu <- tau.post.summary$qq.nu
      mean.nu <- tau.post.summary$mean.nu
      map.nu <- tau.post.summary$map.nu
      mean.phiv <- tau.post.summary$mean.phiv
      map.phiv <- tau.post.summary$map.phiv
      azphi.density <- tau.post.summary$azphi.density
      qq.azphi <- tau.post.summary$qq.azphi
      mean.azphi <- tau.post.summary$mean.azphi
      phi.shmax.density <- tau.post.summary$phi.shmax.density
      qq.phi.shmax <- tau.post.summary$qq.phi.shmax
      mean.phi.shmax <- tau.post.summary$mean.phi.shmax
      map.phi.shmax <- tau.post.summary$map.phi.shmax
   }
   par(spar)

   # (1) nu
   # posterior distribution of nu
   if(!is.na(match("nu.density",plot.list))) {
      cat(paste("Plotting nu density",sep=""))
      plot(epar.vecs$nu, tau.post.summary$nu.density,
          type="l", xlab=expression(nu), ylab="")
      title(expression(paste("Density of ",nu)))
      abline(v=mean.epar[1],col="red")
      par(spar)
   }

   # (2) Phi_S
   if(!is.na(match("phiv",plot.list))) {
      cat(paste("Plotting phiv density",sep=""))
      start.stereogram()
      points.epar.stereogram(map.epar)
      points.epar.stereogram(mean.epar,col="blue")
      title("MAP (red) and mean (blue) angles")
      par(spar)
   }
   
   if(!is.na(match("tau.post",plot.list))) {
      cat(paste("Plotting tau.post",sep=""))
      #plot.tau.post(tau.post, epar=rbind(map.epar,mean.epar))
      plot.tau.post(tau.post)
      par(spar)
   }
   
   if(!is.na(match("contour",plot.list))) {
      cat(paste("Plotting contour map",sep=""))
      contour.tau.post.stereogram(tau.post,ddz=tau.post.summary$ddz)
      #points.epar.stereogram(map.epar)
      #points.epar.stereogram(mean.epar,col="purple")
      if(!is.na(match("radii",plot.list))) {
         radius.stereogram(mean.azphi[1],col="orange",lwd=3)
         radius.stereogram(mean.azphi[2],col="green",lwd=3)
         radius.stereogram(mean.azphi[3],col="blue",lwd=3)
         radius.stereogram(mean.phi.shmax,col="red",lwd=3)
      }
      par(spar)
   }

   # (3) azimuths of the three principal axes
   if(!is.na(match("azphi",plot.list))) {
      cat(paste("Plotting azimuth 1 (min)",sep=""))
      plot(azphi.density[[1]],
           main=expression(paste("Posterior density of ",phi,"(Smin)")))
      abline(v=mean.azphi[1],col="orange",lwd=3)
      cat(paste("Plotting azimuth 2 (mid)",sep=""))
      plot(azphi.density[[2]],
           main=expression(paste("Posterior density of ",phi,"(Smid)")))
      abline(v=mean.azphi[2],col="green",lwd=3)
      cat(paste("Plotting azimuth 3 (max)",sep=""))
      plot(azphi.density[[3]],
           main=expression(paste("Posterior density of ",phi,"(Smax)")))
      abline(v=mean.azphi[3],col="blue",lwd=3)
      
      #radius.stereogram(azphi.map[1],col="orange",lwd=3,lty=3)
      #radius.stereogram(azphi.map[2],col="green",lwd=3,lty=3)
      #radius.stereogram(azphi.map[3],col="blue",lwd=3,lty=3)
      par(spar)
   }

   # (4) SHmax
   if(!is.na(match("phi.shmax.density",plot.list))) {
      cat(paste("Plotting azimuth 4 phi(shmax)",sep=""))
      plot(phi.shmax.density,
           main=expression(paste("Posterior density of ",phi[SHmax])))
      abline(v=map.phi.shmax,col="red",lwd=3)
      par(spar)
   }

   invisible()
}
compare.inversion <- function(stress.out1,stress.out2,
                              draw.grid=F,draw.map.rmat=F) {
  # Compare the outputs of two stress inversions

  par(opar)
  par(mfcol=c(2,2))
  # Draw the two stereonets
  stress.stereonet(stress.out1$invresult,stress.out1$psummary,
                   draw.grid=draw.grid,draw.map.rmat=draw.map.rmat)
  #text(par()$usr[1],par()$usr[4],adj=c(-1,1),lab="1",xpd=T)
  text(-1,1,adj=c(-1,1),lab="1",xpd=T,cex=1.2)
  stress.stereonet(stress.out2$invresult,stress.out2$psummary,
                   draw.grid=draw.grid,draw.map.rmat=draw.map.rmat)
  text(-1,1,adj=c(-1,1),lab="2",xpd=T,cex=1.2)
  par(mar=opar$mar)
  # Draw the SHmax distributions
  shmax.density1 <- stress.out1$psummary$phi.shmax.density
  shmax.density2 <- stress.out2$psummary$phi.shmax.density
  x1 <- shmax.density1$x
  y1 <- shmax.density1$y
  x2 <- shmax.density2$x
  y2 <- shmax.density2$y
  y1 <- y1/sum(y1[-1])
  y2 <- y2/sum(y2[-1])
  plot(180/pi*x1,y1,col="black",type="l",
       xlim=c(0,180),ylim=c(0,max(c(y1,y2))),
       ylab="Posterior density",xlab="Azimuth",
       main="SHmax azimuth")
  lines(180/pi*x2,y2,lty=2)
  legend(par()$usr[2],par()$usr[4],xjust=1,yjust=1,
         legend=c(1,2),lty=c(1,2))

  # Credible intervals for SHmax azimuth
  f1 <- stress.out1$psummary$phi.shmax.density
  qq1 <- 0.5*circular.wquantile(2*f1$x[-1],
                                weights=f1$y[-1],probs=c(0.025,0.5,0.975))
  qq1 <- ((qq1+pi/2)%%pi-pi/2)
  f2 <- stress.out2$psummary$phi.shmax.density
  qq2 <- 0.5*circular.wquantile(2*f2$x[-1],
                                weights=f2$y[-1],probs=c(0.025,0.5,0.975))
  qq2 <- ((qq2+pi/2)%%pi-pi/2)

  # Distribution of the difference
  qq.diff <- ddiff.axial.quantile(x1,y1,y2)
  if(qq.diff["50%"]>pi/2) qq.diff <- qq.diff-pi
  xp <- 180/pi*c(x1-pi,x1)
  yp <- ddiff.axial(y1,y2)
  yp <- rep(yp,2)
  plot(xp,yp,type="l",
       xlab="Angular Difference",ylab="p(diff)",
       main=paste("Difference in SHmax azimuth\n",
                  sprintf("%.1f (%.1f,%.1f)",
                          180/pi*qq.diff[2],
                          180/pi*qq.diff[1],
                          180/pi*qq.diff[3]),sep=""))
  abline(v=qq.diff*180/pi,lty=2)

  # Finished
  invisible(list(qq1=qq1,qq2=qq2,qq.diff=qq.diff))
}
stress.stereonet <- function(invresult,psummary,newpage=F,
                             draw.grid=T,draw.map.rmat=T,
                             cols=c("red","green","blue"),
                             store.stereonet.file=NULL) {
  # Draw the stereonet of the stress tensor axes
  
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
  ststore <- list(NA,NA,NA)
  if(newpage) {
     par(opar)
     par(mar=c(0.5,0.5,0.5,0.5))
     par(mfrow=c(2,2))
     par(mfg=c(1,1,2,2))
  }
  start.stereogram(plot.centre=T)               # Construct basemap
  if(draw.grid) stereogram.grid()               # Add grid
  if(!all(is.na(ddz[[1]]))) {
    ststore[[1]] <- contour.stereogram(phivec,thetavec,ddz[[1]],
                                       add=T,col=cols[1]) # S1 (Smin)
  } else {
    ststore[[1]] <- NA
  }
  if(!all(is.na(ddz[[2]]))) {
    ststore[[2]] <- contour.stereogram(phivec,thetavec,ddz[[2]],
                                       add=T,col=cols[2]) # S2 (Smid)
  } else {
    ststore[[2]] <- NA
  }
  if(!all(is.na(ddz[[3]]))) {
    ststore[[3]] <- contour.stereogram(phivec,thetavec,ddz[[3]],
                                       add=T,col=cols[3]) # S3 (Smax)
  } else {
    ststore[[3]] <- NA
  }
  # Note that in the following we must interpret the first (third) column
  # of the rotation matrix as S3 (S1) and colour it accordingly
  points.rmat.stereogram(mean.rmat[,3],pch=19,col=cols[1])
  points.rmat.stereogram(mean.rmat[,2],pch=19,col=cols[2])
  points.rmat.stereogram(mean.rmat[,1],pch=19,col=cols[3])
  if(draw.map.rmat) {
     points.rmat.stereogram(map.rmat[,3],pch=21,col=cols[1],bg="white")
     points.rmat.stereogram(map.rmat[,2],pch=21,col=cols[2],bg="white")
     points.rmat.stereogram(map.rmat[,1],pch=21,col=cols[3],bg="white")
  }
  radius.stereogram(psummary$mean.phi.shmax,col="black",lwd=1,lty="dashed")
  radius.stereogram(psummary$mean.phi.shmax+pi,col="black",lwd=1,lty="dashed")

  # Write out the stereogram contour line segments
  if(!is.null(store.stereonet.file)) {
    # store stereogram file
    for(i in 1:3) {
      # break up any line segments with NA values before writing
      ofile <- paste(store.stereonet.file,"_S",i,".dat",sep="")
      cat(paste("Writing",ofile,"\n"))
      sink(ofile)
      ststore[[i]] <- delist1(lapply(ststore[[i]],split.na))
      lapply(ststore[[i]], function(x) {
         cat("####\n")
         x <- matrix(x,ncol=2)
         write.table(x,append=T,sep=",",row.names=F,col.names=F)
         invisible()
      })
      sink()
    }
  }

  # finished
  invisible()
}

#############################################################################
# Alternative methods
#############################################################################
# 1 Michael inversion

michael.invert <- function(sdr,verbose=F,flip=F) {
   # Michael (1984) inversion
   if(is.data.frame(sdr)) sdr <- as.matrix(sdr)
   dim(sdr) <- c(length(sdr)/3,3)
   if(flip) sdr <- rbind(sdr,flip.sdr(sdr))
   amat.f <- function(n) {
     # the A matrix for a unit vector n
     amat <- array(
     c(n[1]-n[1]^3+n[1]*n[3]^2, -n[2]*n[1]^2+n[2]*n[3]^2, -n[3]*n[1]^2-n[3]+n[3]^3,
            n[2]-2*n[2]*n[1]^2,       n[1]-2*n[1]*n[2]^2,        -2*n[1]*n[2]*n[3],
            n[3]-2*n[3]*n[1]^2,        -2*n[1]*n[2]*n[3],       n[1]-2*n[1]*n[3]^2,
      -n[1]*n[2]^2+n[1]*n[3]^2,  n[2]-n[2]^3+n[2]*n[3]^2, -n[2]^2*n[3]-n[3]+n[3]^3,
             -2*n[1]*n[2]*n[3],       n[3]-2*n[3]*n[2]^2,       n[2]-2*n[2]*n[3]^2
          ), dim=c(3,5))
     return(amat)
   }
   # prepare the data
   phiv <- convert.sdr.phiv(sdr)
   nphiv <- nrow(phiv)
   rmats <- rotmat.phiv.asrow(phiv)
   nvecs <- rmats[,7:9] # fault normals (n)
   uvecs <- rmats[,1:3] # slip vectors (u)
   dim(nvecs) <- c(nphiv,3)
   amats <- array( apply(nvecs,1,amat.f), dim=c(3,5,nphiv))
   amats <- aperm(amats, c(1,3,2))
   dim(amats) <- c(3*nphiv,5)
   all.uvecs <- as.vector(t(uvecs)) # stack up the slip vectors

   # estimate the deviatoric stress tensor
   xtx.inv <- try(solve(t(amats)%*%amats),silent=T)
   if(class(xtx.inv)=="try-error") {
     # XtX is not invertible: no solution exists
     epar <- c(NA,NA,NA,NA)
     return(epar)
   } 
   sigvec <- xtx.inv %*% t(amats) %*% all.uvecs
   sigmat <- array(c( sigvec[1], sigvec[2], sigvec[3],
                      sigvec[2], sigvec[4], sigvec[5],
                      sigvec[3], sigvec[5], -(sigvec[1]+sigvec[4])), dim=c(3,3))
   fit.uvecs <- amats %*% sigvec
   #vv <- as.numeric(var(all.uvecs-fit.uvecs))
   #cov.sigvec <- xtx.inv*vv

   # properties of the fit
   fuvecs <- t(array(fit.uvecs,dim=c(3,nphiv)))
   tau <- apply(fuvecs,1,function(x) sqrt(sum(x^2))) # lengths of fitted vectors
   beta <- acos( apply(fuvecs*uvecs,1,sum)/tau ) # angles between actual and fitted
   dev <- diag(nvecs %*% sigmat %*% t(nvecs)) # normal stress from deviatoric component of the stress tensor
   ival <- -tau/0.8 - dev # remainder = normal stress from isotropic component
   if(verbose) {
     cat(paste("mean tau=",round(mean(tau),2)," (sd=",
               round(sd(tau),2),")\n",sep=""))
     cat(paste("mean beta=",round(180/pi*mean(beta),2)," (sd=",
               round(180/pi*sd(beta),2),")\n",sep=""))
     cat(paste("mean ival=",round(mean(ival),2)," (sd=",
               round(sd(ival),2),")\n",sep=""))
   }
   
   # principal axes
   eigen.sigmat <- eigen(sigmat)
   idx <- rev(order(eigen.sigmat$values)) # sort into descending order
   # puts the most compressive - most negative - stress along the z-axis
   evals <- eigen.sigmat$values[idx]
   
   #print("sigvec"); print(sigvec)
   #print("sigmat"); print(sigmat)
   #print("eigen.sigmat"); print(eigen.sigmat)
   #print("evals"); print(evals)
   
   rmat.gs <- eigen.sigmat$vectors[,idx]
   # ensure that this is a right handed matrix: det=+1
   # reverse the middle vector if not
   if(det(rmat.gs)<0) rmat.gs[,2] <- -rmat.gs[,2]
   phiv.gs <- extract.phiv.rotmat(rmat.gs)
   # ordered x=1=least to z=3=most compressive
   phi <- (evals[2]-evals[1])/(evals[3]-evals[1])
   nu <- 1-phi
   epar <- c(nu, phiv.gs)
   return(epar)
}   
boot.michael.invert.statistic <- function(sdr,idx,flip=F) {
  michael.invert(sdr[idx,],flip=flip)
}
boot.epar.michael <- function(sdr1, R=50, flip=F) {
   # boottrap estimates of epar
   require(boot)
   b1 <- boot(sdr1, boot.michael.invert.statistic, R=R, flip=flip)
   # centre each bootstrap estimate
   dd <- b1$t
   dd[,1] <- (b1$t[,1]-b1$t0[1])
   dd[,2] <- (b1$t[,2]-b1$t0[2]+pi  )%%(2*pi) - pi
   dd[,3] <- (b1$t[,3]-b1$t0[3]+pi/2)%%(pi)   - pi/2
   dd[,4] <- (b1$t[,4]-b1$t0[4]+pi  )%%(2*pi) - pi
   retval <- cbind(b1$t0,
                   t(apply(dd,2,function(dd) c(mean(dd),sd(dd)))),
                   b1$t0+t(apply(dd,2,function(dd) {
                          (quantile(dd,
                   probs=c(0.025,0.05,0.10,0.50,0.90,0.95,0.975)))}))
                   )
   dimnames(retval)[[1]] <- c("nu","phi","theta","psi")
   dimnames(retval)[[2]] <- c("estimate","bias","stderr",
                      "2.5%","5%","10%","50%","90%","95%","97.5%")
   return(retval)
}
boot.shmax.michael.invert.statistic <- function(sdr,idx,flip=F) {
  epar <- michael.invert(sdr[idx,],flip=flip)
  phi.shmax <- phi.shmax.f(as.earth(epar))
  return(phi.shmax)
}
boot.shmax.michael <- function(sdr1, R=50, flip=F) {
   # boottrap estimates of SHmax
   require(boot)
   b1 <- boot(sdr1, boot.shmax.michael.invert.statistic, R=R, flip=flip)
   # centre each bootstrap estimate, and put in the range (-pi,+pi)
   dd <- (b1$t-b1$t0+pi/2)%%(pi) - pi/2
   retval <- c(b1$t0,mean(dd),sd(dd),
               (b1$t0+quantile(dd,
                   probs=c(0.025,0.05,0.10,0.50,0.90,0.95,0.975)))%%(2*pi))
   names(retval) <- c("estimate","bias","stderr",
                      "2.5%","5%","10%","50%","90%","95%","97.5%")
   return(retval)
}

make.michael.cc <- function(srcdir="michael-algorithm/src/",bindir="bin/") {
  # compile the Michael programs
  curdir <- getwd()
  on.exit(setwd(curdir))
  # go to the source directory
  setwd(paste(curdir,srcdir,sep=.Platform$file.sep))
  # create slick and slfast
  system("make -f makeslick")
  system("make -f makeslfast")
  # move these to the bin directory
  mv.cmd <- ifelse(.Platform$OS.type=="windows","ren","mv")
  system(paste(mv.cmd,"slick",paste(curdir,bindir,sep=.Platform$file.sep)))
  system(paste(mv.cmd,"slfast",paste(curdir,bindir,sep=.Platform$file.sep)))
  # return to original working directory
  setwd(curdir)
  # finished
  invisible()
}

michael.invert.cc <- function(sdr, verbose=F, flip=F,
                              method="slick",
                              bindir="bin/",tmpdir="tmp/") {
   # Michael (1984) inversion - using Michael's own program
   if(is.data.frame(sdr)) sdr <- as.matrix(sdr)
   dim(sdr) <- c(length(sdr)/3,3)
   if(flip) sdr <- rbind(sdr,flip.sdr(sdr))
   # convert sdr to dip direction/dip/rake, and convert back to degrees
   ddr.deg <- 180/pi*cbind( sdr[,1]+pi/2, sdr[,2], sdr[,3] )

   if(method=="slick") {
      # write dipdir/dip/rake out to a temporary file
      write.table(ddr.deg,paste(tmpdir,"slick.input",sep=""),row.names=F,col.names=T)
      # run slick
      #system("./michael/slick slick.input")
      system(paste(bindir,"slick ",tmpdir,"slick.input",sep=""))
      # read back the output file
      sdat <- read.table(paste(tmpdir,"slick.input.oput",sep=""),sep="\n",as.is=T)
      system(paste("\\rm ",tmpdir,"slick.input",sep=""))
      system(paste("\\rm ",tmpdir,"slick.input.oput",sep=""))
      if(verbose) print(sdat)
      # extract values
      tmp <- sdat[4:6,]
      stress.tensor <- array(as.numeric(array(unlist(strsplit(tmp,split=" ")),
                                              dim=c(6,3))[c(1,3,5),]),dim=c(3,3))
      tmp <- sdat[8:10,]
      tmp <- t(array(as.numeric(array(unlist(strsplit(tmp,split=" ")),
                                     dim=c(11,3))[c(1,3,5,7,9,11),]),dim=c(6,3)))
      stress.eigen <- tmp[,1]
      stress.vecs <- t(tmp[,2:4])
      stress.direct <- tmp[,5]
      stress.plunge <- tmp[,6]
      variance <- as.numeric(strsplit(sdat[11,],split="=")[[1]][2])
      phivalue <- as.numeric(strsplit(sdat[12,],split="=")[[1]][2])
   
      # make some checks:
      # (1) reorder the vectors and directions from least to most compressive
      #     note that Michael takes compressive stress as negative
      idx <- rev(order(stress.eigen)) # from least to most compressive
      stress.eigen <- stress.eigen[idx]
      stress.tensor <- stress.tensor[,idx]
      stress.tensor <- stress.tensor[idx,]
      stress.vecs <- stress.vecs[,idx]
      stress.direct <- stress.direct[idx]
      stress.plunge <- stress.plunge[idx]
      # (2) convert the eigen values into a phi value
      check.phi <- (stress.eigen[2]-stress.eigen[1])/(stress.eigen[3]-stress.eigen[1])
      c(phivalue,check.phi) # correct
      # (3) create eigenvectors from direction and plunge, and compare
      stress.direct <- stress.direct*pi/180
      stress.plunge <- stress.plunge*pi/180
      a <- -rbind( cos(-stress.plunge)*sin(stress.direct),
                   cos(-stress.plunge)*cos(stress.direct),
                   sin(-stress.plunge) )
      cbind(stress.vecs,a)
      
      # (4) coordinates in Michael are (E,N,UP) whereas we have (N,E,DOWN)
      #     so swap first coordinates and reverse the 3rd
      rmat.gs <- exchaxis.3()%*%stress.vecs
      rmat.gs <- rmat.gs%*%invaxis.1()
      
      # convert these parameters to epar
      nu <- 1-phivalue
      phiv.gs <- extract.phiv.rotmat(rmat.gs)
      epar <- c(nu,phiv.gs)
      
   } else if(method=="slfast") {
      # write dipdir/dip/rake out to a temporary file
      write.table(ddr.deg,paste(tmpdir,"slfast.input",sep=""),row.names=F,col.names=T)
      # run slfast
      system("\\rm slfast.input.slboot")
      #system("./michael/slfast slfast.input")
      system(paste(bindir,.Platform$OS.type,"slfast slfast.input",sep=""))
      # read back the output file
      sdat <- as.matrix(read.table(paste(tmpdir,"slfast.input.slboot",sep=""),as.is=T))
      system(paste("\\rm ",tmpdir,"slfast.input",sep=""))
      system(paste("\\rm ",tmpdir,"slfast.input.slboot",sep=""))
      if(verbose) print(sdat)
      # extract values
      stress.tensor <- array(sdat[1,c(2,3,4, 3,5,6, 4,6,7)], dim=c(3,3))
      tmp <- eigen(stress.tensor)
      stress.eigen <- tmp$values
      stress.vecs <- tmp$vectors %*% invaxis.2()
      stress.direct <- sdat[2,c(6,4,2)]
      stress.plunge <- sdat[2,c(7,5,3)]
      variance <- sdat[1,1]
      phivalue <- sdat[2,1]
   
      # make some checks:
      # (1) reorder the vectors and directions from least to most compressive
      #     note that Michael takes compressive stress as negative
      idx <- rev(order(stress.eigen)) # from least to most compressive
      stress.eigen <- stress.eigen[idx]
      stress.tensor <- stress.tensor[,idx]
      stress.tensor <- stress.tensor[idx,]
      stress.vecs <- stress.vecs[,idx]
      stress.direct <- stress.direct[idx]
      stress.plunge <- stress.plunge[idx]
      # (2) convert the eigen values into a phi value
      check.phi <- (stress.eigen[2]-stress.eigen[1])/(stress.eigen[3]-stress.eigen[1])
      c(phivalue,check.phi) # correct
      # (3) create eigenvectors from direction and plunge, and compare
      stress.direct <- stress.direct*pi/180
      stress.plunge <- stress.plunge*pi/180
      a <- -rbind( cos(-stress.plunge)*sin(stress.direct),
                   cos(-stress.plunge)*cos(stress.direct),
                   sin(-stress.plunge) )
      cbind(stress.vecs,a)
      
      # (4) coordinates in Michael are (E,N,UP) whereas we have (N,E,DOWN)
      #     so swap first coordinates and reverse the 3rd
      rmat.gs <- exchaxis.3()%*%stress.vecs
      rmat.gs <- rmat.gs%*%invaxis.1()
      
      # convert these parameters to epar
      nu <- 1-phivalue
      phiv.gs <- extract.phiv.rotmat(rmat.gs)
      epar <- c(nu,phiv.gs)
      
   } else {
     epar <- c(NA,NA,NA,NA)
     stop("method not recognised")
   }
   return(epar)
}

##!!== checked to here
#############################################################################
create.tau.qfunc <- function(tau,mepar=c(11,51,31,11)) {
   # create a new tabulation for a particular value of tau

   cat("Creating vectors...")
   epar.min <- c(0, 0,    0,  0)
   epar.max <- c(1, 2*pi, pi, 2*pi)
   epar.smin <- c(0, 0,    0,    0)
   epar.smax <- c(1, 2*pi, pi/2, pi)
   nepar <- c(mepar[1:2], 2*mepar[3:4]-1)
   epar.vecs <- list()
   for(i in 1:4) epar.vecs[[i]] <- seq(from=epar.min[i],to=epar.max[i],
                                       length=nepar[i])
   names(epar.vecs) <- c("nu","phi","theta","psi")
   epar.grid <- expand.grid(epar.vecs[[1]],epar.vecs[[2]],
                            epar.vecs[[3]],epar.vecs[[4]])
   names(epar.grid)  <- c("nu","phi.s","theta.s","psi.s")
   epar.grid <- as.matrix(epar.grid)

   time.start <- Sys.time()
   cat(paste("Doing tau=",tau,"; started at",time.start,"....",sep=""))
   sqmat1.tau <- qfunc(epar.sgrid, tau)
   sqmat1.tau.alt <- qfunc(epar.sgrid, tau, flip=T)
   dim(sqmat1.tau) <- mepar
   dim(sqmat1.tau.alt) <- mepar
   tau.qfunc <- list(nepar=nepar, epar.vecs=epar.vecs,
                      qvals1=extend.4dgrid(sqmat1.tau),
                      qvals2=extend.4dgrid(sqmat1.tau.alt))
   tau.qfunc$qvals0 <- tau.qfunc$qvals1+tau.qfunc$qvals2
   time.now <- Sys.time()
   cat(paste("Completed at",time.now,"\n"))
   cat(paste("Duration: ",time.now-time.start,"\n"))
   
   return(tau.qfunc)
}
write.tau.qfunc.bin <- function(tau.qfunc.obj,filename,...) {
   writeBin(unlist(tau.qfunc.obj),filename,...)
   invisible()
}
read.tau.qfunc.bin <- function(filename,...) {
   tt <- readBin(filename,what=0,n=4)
   ntt <- length(tt) + sum(tt) + 3*prod(tt)
   k1 <- length(tt)+sum(tt)
   k2 <- prod(tt)
   tt <- readBin(filename,what=0,n=ntt,...)
   tt1 <- list( nepar=tt[1:4],
             epar.vecs=list(nu=tt[4+1:tt[1]],
                           phi=tt[4+tt[1]+1:tt[2]],
                         theta=tt[4+tt[1]+tt[2]+1:tt[3]],
                           psi=tt[4+tt[1]+tt[2]+tt[3]+1:tt[4]]),
             qvals1=array(tt[k1+1:k2],dim=tt[1:4]),
             qvals2=array(tt[k1+k2+1:k2],dim=tt[1:4]),
             qvals0=array(tt[k1+k2+k2+1:k2],dim=tt[1:4]))
   return(tt1)
}

invert.sdr <- function(sdr,errs=10*pi/180,weights=1,
                       tauvals=c(1000,60,30,20,10,5,2,1),
                       tabdir="tab",verbose=T) {
   # combine a set of focal mechanisms specified by sdr (in radians)
   # with specified errors (also in radians)
   # requires a set of pretabulated single observation
   # posterior distributions
   # all observations are assumed to be ambiguous
   errvals <- exp(3.9155-0.5659*log(tauvals)) * pi/180
   nsdr <- length(sdr)/3
   dim(sdr) <- c(nsdr,3)
   if(length(weights)==1) weights <- rep(weights,nsdr)
   if(length(errs)==1) errs <- rep(errs,nsdr)
   
   for(i in 1:nsdr) {
      if(verbose) cat(paste("doing",i,"..."))

      wgt1 <- weights[i]
      if(verbose) cat(sprintf("weight=%10.8f;",wgt1))
      
      err1 <- errs[i]
      tau <- tauvals[which.min(abs(err1-errvals))]
      if(verbose) cat(paste("Obs. error=",err1*180/pi,"degrees, using tau=",tau))
      fname.bin <- paste(tabdir,"/tau",tau,".tab.qfunc.Rbin",sep="")
      tau.qfunc <- read.tau.qfunc.bin(fname.bin)
      
      sdr1 <- sdr[i,]
      phiv.g <- convert.sdr.phiv(sdr1)
   
      # retabulate 
      tau.qfunc.alt <- retabulate.qflist(phiv.g,tau.qfunc,qv0.only=T)
   
      if(i==1) {
        qfcombined <- tau.qfunc.alt
        qfcombined$qvals0 <- qfcombined$qvals0^wgt1
      } else {
        qfcombined$qvals0 <- qfcombined$qvals0 * (tau.qfunc.alt$qvals0)^wgt1
      }
      print(c(range(qfcombined$qvals0),
              NA,
              length(qfcombined$qvals0),
              length(qfcombined$qvals0[qfcombined$qvals0>1e-10])))##!!==
      # rescale to prevent over/underflows
      qfcombined$qvals0 <- qfcombined$qvals0/max(qfcombined$qvals0)
      
      if(verbose) cat("...done\n")
   }
   return(qfcombined)
}
#############################################################################
s1.contour.plot <- function(result,ofile,view=F) {
  # Draw a contour plot of the S1 axis
  
  # Marginalise the pdf to obtain p(phi,cos(theta)|data)
  zz <- sum.mat.trap(result$qvals0,c(2,3))
  # Plotting vectors
  vphi <- result$epar.vecs$phi
  vtheta <- result$epar.vecs$theta
  # Set graphics parameters
  par(opar)
  # Make plot
  epsfile <- open.eps(ofile)           
  contour.stereogram(vphi,vtheta,zz,plot.centre=T)
  close.eps(epsfile,view=view)
  invisible()
}
s123.contour.plot <- function(result,ofile,view=F) {
  # Make a contour plot of 3 principal axes at once

  # Plotting vectors
  vphi <- result$epar.vecs$phi
  vtheta <- result$epar.vecs$theta
  
  # Copy the results - for Smax[S1] (this is the z-axis)
  dd0 <- result				# Make a duplicate of the results
  # Swap the y and z axes - i.e. Smid[S2] is now the z-axis
  dd1 <- retabulate.qflist(c(0,0,0),dd0,exchange=1) 
  # Swap the x and z axes - i.e. Smin[S3] is now the z-axis
  dd2 <- retabulate.qflist(c(0,0,0),dd0,exchange=2) 
  # Density of the z-axis
  ddz0 <- sum.mat.trap(dd0$qvals0,c(2,3))	# Marginalise over psi and nu for (S1)
  ddz1 <- sum.mat.trap(dd1$qvals0,c(2,3))	# Ditto (S2)
  ddz2 <- sum.mat.trap(dd2$qvals0,c(2,3))	# Ditto (S3)

  # Find the posterior mode
  pp0 <- unlist(expand.grid(vphi,vtheta)[as.vector(ddz0)==max(ddz0),])
  pp1 <- unlist(expand.grid(vphi,vtheta)[as.vector(ddz1)==max(ddz1),])
  pp2 <- unlist(expand.grid(vphi,vtheta)[as.vector(ddz2)==max(ddz2),])

  # Flip so that both directions of the axes are identified
  pp0 <- rbind(pp0, c(pp0[1]+pi,pi-pp0[2]))
  pp1 <- rbind(pp1, c(pp1[1]+pi,pi-pp1[2]))
  pp2 <- rbind(pp2, c(pp2[1]+pi,pi-pp2[2]))

  # Create the graph
  par(opar)
  epsfile <- open.eps(ofile)
  par(mar=c(0,0,0,0))
  cols <- c("blue","green","orange")
  acol <- "red"
  start.stereogram(plot.centre=F)
  contour.stereogram(vphi,vtheta,ddz0,add=T,col=cols[1])
  contour.stereogram(vphi,vtheta,ddz1,add=T,col=cols[2])
  contour.stereogram(vphi,vtheta,ddz2,add=T,col=cols[3])
  points.phiv2.stereogram(pp0,pch=21,cex=2,bg=cols[1])
  points.phiv2.stereogram(pp1,pch=21,cex=2,bg=cols[2])
  points.phiv2.stereogram(pp2,pch=21,cex=2,bg=cols[3])
  text.phiv2.stereogram(pp0,label=expression(S[max]),cex=2.5,pos=3)
  text.phiv2.stereogram(pp1,label=expression(S[mid]),cex=2.5,pos=3)
  text.phiv2.stereogram(pp2,label=expression(S[min]),cex=2.5,pos=3)
  close.eps(epsfile,view=view)

  # Finished
  invisible()
}
nu.contour.plots <- function(result,ofile,view=F) {
  # Posterior distributions of nu (stress ratio)

  nepar <- result$nepar
  nuvec <- result$epar.vecs$nu
  phivec <- result$epar.vecs$phi
  thetavec <- result$epar.vecs$theta
  psivec <- result$epar.vecs$psi

  # sin(theta) array for scaling
  sin.theta <- sin(thetavec)
  a.sin.theta <- array( rep(rep(sin.theta,
           times=rep(nepar[1]*nepar[2],nepar[3])),nepar[4]), dim=nepar)
  postpdf1 <- result$qvals0         # density over cos(theta)
  postpdf2 <- postpdf1*a.sin.theta  # density over theta
     
  par(opar)
  epsfile <- open.eps(ofile,width=6,aspect=0.7)
  par(mfrow=c(1,3))

  # p(nu,phi\data)
  zz <- sum.mat.trap(postpdf2,c(1,2))
  contour(nuvec,phivec,zz,
          xlab=expression(nu),ylab=expression(phi),
          cex.lab=1.5)
  # p(nu,cos(theta)|data)
  zz <- sum.mat.trap(postpdf1,c(1,3))
  contour(nuvec,thetavec,zz,
          xlab=expression(nu),ylab=expression(theta),
          cex.lab=1.5)
  # p(nu,psi|data)
  zz <- sum.mat.trap(postpdf2,c(1,4))
  contour(nuvec,psivec,zz,
          xlab=expression(nu),ylab=expression(psi),
          cex.lab=1.5)
  
  close.eps(epsfile,view=view)  
  par(mfrow=c(1,1))

  # Finished
  invisible()
}
stress.inversion <- function(infile,fname=NULL,
                             inputdir="Input",outputdir="Output",
                             tabdir="tab",
                             skip=1,err=NULL,sep=",",verbose=F,
                             save=T,do.graphics=T,text.out=T,
                             opar=NULL) {
  # takes an input file <filename> 
  # containing strike,dip,rake,err,wgt (comma separated)
  #    with all angles in ***degrees***
  #    if wgt is missing it is set to 1 for each observation
  # if err=NULL then err is in the file
  #    otherwise err=constant value for each event
  # skip    = number of rows to skip before data start
  # sep     = separator in the data file
  # verbose = print comments during run?
  # save    = save the output object to an R dump file?
  # do.graphics = make a graphical summary?
  # text.out    = make a text summary?

  # file stem name
  if(is.null(fname)) fname <- unlist(strsplit(infile,split=".",fixed=T))[1]
  if(verbose) cat(paste("File stem name is ",fname,"\n",sep=""))
  # read the input file
  infile <- paste(inputdir,infile,sep="/")
  sdr <- read.table(infile,skip=skip,sep=sep,as.is=T)
  nc <- ncol(sdr)
  if(nc==5) {
    if(verbose) cat(paste("err and weights specified in file\n",sep=""))
    names(sdr) <- c("strike","dip","rake","err","weights")
  } else if(nc==4) {
    if(!is.null(err)) {
      names(sdr) <- c("strike","dip","rake","weights")
      sdr$err <- err
      if(verbose) cat(paste("weights specified in file\n",sep=""))
    } else {
      names(sdr) <- c("strike","dip","rake","err")
      sdr$weights <- 1
      if(verbose) cat(paste("weights not specified in file\n",sep=""))
    }
  } else if(nc==3) {
    if(is.null(err)) error("Format of input file incorrect.")
    names(sdr)[1:3] <- c("strike","dip","rake")
    sdr$err <- err
    sdr$weights <- 1
    if(verbose) cat(paste("weights not specified in file\n",sep=""))
  }
  if(verbose) {
    cat(paste("Read ",nrow(sdr)," earthquakes from file ",infile,"\n",sep=""))
    print(sdr)
  }
  weights <- sdr$weights
  sdr <- as.matrix(sdr[,c("strike","dip","rake","err")])*pi/180

  # inversion
  # do the inversion
  if(verbose) cat("Starting the inversion now:\n")
  invresult <- invert.sdr(sdr[,1:3],errs=sdr[,4],weights=weights,
                          tabdir=tabdir,verbose=verbose)

  print(names(invresult))##!!==
  print(range(invresult$qvals0))##!!==
  #print(range(invresult$qvals1))
  #print(range(invresult$qvals2))
  
  if(verbose) cat("Complete.\n")

  # Michael (1984,1987) stress inversion method
  michael.epar <- michael.invert(sdr[,1:3])
  # Shmax for this
  michael.shmax <- phi.shmax.f(as.earth(michael.epar))

  # Trap errors
  if(any(is.nan(range(invresult$qvals0)))) {
    warning(paste(fname,": Posterior pdf has numerical errors (NaN) -- possible underflow -- Exiting"))
    retval <- list(fname=fname,
                   sdr=sdr,
                   invresult=invresult,
                   psummary=list(michael.epar=michael.epar,
                                 michael.shmax=michael.shmax))
    return(retval)
  }
  
  # posterior summary
  if(verbose) cat("Computing posterior summaries:\n")
  psummary <- tau.posterior.summary(invresult,tabdir=tabdir,verbose=verbose)
  psummary$michael.epar <- michael.epar
  psummary$michael.shmax <- michael.shmax
  if(verbose) cat("Complete.\n")

  # Make graphical output
  if(do.graphics) {
     epsfile <- open.eps(paste(outputdir,paste(fname,".eps",sep=""),sep="/"),
                         width=9,height=8)
     single.page.summary(invresult,psummary,opar=opar)
     close.eps(epsfile)
  }

  # Make text output files
  if(text.out) {
    # Output strike/dip grid
    phivec <- invresult$epar.vecs$phi
    thetavec <- invresult$epar.vecs$theta
    dphi <- (max(phivec)-min(phivec))/(length(phivec)-1)
    dtheta <- (max(thetavec)-min(thetavec))/(length(thetavec)-1)
    outfile <- paste(outputdir,paste(fname,".s123grid.dat",sep=""),sep="/")
    cat("# Phi/Theta Grid\n",file=outfile)
    cat("# Phi (strike):\n",file=outfile,append=T)
    cat(paste(length(phivec),",",180/pi*min(phivec),",",
              180/pi*max(phivec),",",180/pi*dphi,
              "\n",sep=""),file=outfile,append=T)
    cat("# Theta (dip+90):\n",file=outfile,append=T)
    cat(paste(length(thetavec),",",180/pi*min(thetavec),",",
              180/pi*max(thetavec),",",180/pi*dtheta,
              "\n",sep=""),file=outfile,append=T)
    cat("# Phi:\n",file=outfile,append=T)
    cat(paste(phivec*180/pi,collapse="\n"),file=outfile,append=T)
    cat("\n",file=outfile,append=T)
    cat("# Theta:\n",file=outfile,append=T)
    cat(paste(thetavec*180/pi,collapse="\n"),file=outfile,append=T)
    cat("\n",file=outfile,append=T)
    
    # Posterior density of S1(Smax)=3, S2(Smid)=2, S3(Smin)=1
    outfile <- paste(outputdir,paste(fname,".s1density.dat",sep=""),sep="/")
    odf <- psummary$ddz[[3]] # S1(Smax)=3
    write.table(odf,outfile,col.names=F,row.names=F,sep=",",quote=F)
    outfile <- paste(outputdir,paste(fname,".s2density.dat",sep=""),sep="/")
    odf <- psummary$ddz[[2]] # S2(Smid)=2
    write.table(odf,outfile,col.names=F,row.names=F,sep=",",quote=F)
    outfile <- paste(outputdir,paste(fname,".s3density.dat",sep=""),sep="/")
    odf <- psummary$ddz[[1]] # S3(Smin)=1
    write.table(odf,outfile,col.names=F,row.names=F,sep=",",quote=F)
    
    # Output posterior density of nu as a text file
    outfile <- paste(outputdir,paste(fname,".nudensity.dat",sep=""),sep="/")
    odf <- data.frame(nu=invresult$epar.vecs$nu, density=psummary$nu.density)
    write.table(odf,outfile,col.names=T,row.names=F,sep=",",quote=F)
    # Output posterior density of phi.shmax as a text file
    outfile <- paste(outputdir,paste(fname,".shmaxdensity.dat",sep=""),sep="/")
    odf <- data.frame(phi=psummary$phi.shmax.density$x*180/pi,
                      density=psummary$phi.shmax.density$y)
    write.table(odf,outfile,col.names=T,row.names=F,sep=",",quote=F)

    # Estimates and credible intervals for 1D parameters
    # Azimuths and the stress ratio
    outfile <- paste(outputdir,paste(fname,".1dparameters.dat",sep=""),sep="/")
    odf <- array(NA,dim=c(5,5))
    dimnames(odf) <- list(c("S1","S2","S3","Shmax","nu"),
                          c("mean","map","median","2.5%","97.5%"))
    odf[1:3,"mean"] <- psummary$mean.azphi[3:1]
    odf[1:3,"map"] <- psummary$map.azphi[3:1]
    odf[1:3,c("median","2.5%","97.5%")] <- t(psummary$qq.azphi[c("50%","2.5%","97.5%"),3:1])
    odf[4,"mean"] <- psummary$mean.phi.shmax
    odf[4,"map"] <- psummary$map.phi.shmax
    odf[4,"median"] <- psummary$qq.phi.shmax[c("50%")]
    odf[4,"2.5%"] <- psummary$qq.phi.shmax[c("2.5%")]
    odf[4,"97.5%"] <- psummary$qq.phi.shmax[c("97.5%")]
    odf <- 180/pi*odf
    odf[5,"mean"] <- psummary$mean.nu
    odf[5,"map"] <- psummary$map.nu
    odf[5,"median"] <- psummary$qq.nu[c("50%")]
    odf[5,"2.5%"] <- psummary$qq.nu[c("2.5%")]
    odf[5,"97.5%"] <- psummary$qq.nu[c("97.5%")]
    odf <- data.frame(parameter=row.names(odf),odf)
    write.table(odf,outfile,col.names=T,row.names=F,sep=",",quote=F)
    # Strike and Dip of axes
    outfile <- paste(outputdir,paste(fname,".2dparameters.dat",sep=""),sep="/")
    odf <- data.frame(parameter=c("nu","phi","theta","psi"),
                      mean=psummary$mean.epar,
                      map=psummary$map.epar,
                      michael=psummary$michael.epar)
    # convert these to phi/theta for the principal axes
    rmat.mean <- rotmat.phiv(psummary$mean.epar[2:4])
    rmat.map <- rotmat.phiv(psummary$map.epar[2:4])
    rmat.michael <- rotmat.phiv(psummary$michael.epar[2:4])
    phiv2.mean <- apply(rmat.mean[,3:1],2,extract.phiv2.nvec)
    phiv2.map <- apply(rmat.map[,3:1],2,extract.phiv2.nvec)
    phiv2.michael <- apply(rmat.michael[,3:1],2,extract.phiv2.nvec)
    phiv2.mat <- data.frame(parameter=c("S1:Phi","S1:Theta",
                              "S2:Phi","S2:Theta","S3:Phi","S3:Theta"),
                            mean=as.vector(phiv2.mean),
                            map=as.vector(phiv2.map),
                            michael=as.vector(phiv2.michael))
    odf <- rbind(odf,phiv2.mat)
    odf[-1,2:4] <- 180/pi*odf[-1,2:4]
    write.table(odf,outfile,col.names=T,row.names=F,sep=",",quote=F)
  }

  # return value
  retval <- list(fname=fname,
                 sdr=sdr,
                 invresult=invresult,
                 psummary=psummary)
  # save this
  if(save) {
    if(verbose) cat("Saving output\n")
    outfile <- paste(outputdir,paste(fname,".out",sep=""),sep="/")
    save(retval,file=outfile)
  }
  
  # Finished
  if(verbose) cat("Finished.\n")
  return(retval)
}
#############################################################################
#############################################################################
