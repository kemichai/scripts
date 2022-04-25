#-----------------------------------------------------------------------------
# Post processing
#-----------------------------------------------------------------------------
reload()

# Using the output object created in R
names(stress.out)
single.page.summary(stress.out$invresult,stress.out$psummary,opar=opar)

#-----------------------------------------------------------------------------
# Using the output text files
outputdir <- "outdata"
#rname <- "dixie"
rname <- "GEON2"
#rname <- "montserrat_period1"

# 1. rname.eps = single page graphical summary of the inversion
# 2. rname.1dparameters.dat
#        5 line comma separated file:
#          name, mean, map, median, 10%, 90%
#    name = S1,S2,S3 (for max,mid,min), SHmax, nu
#    mean = posterior mean
#    map  = posterior mode
#    median = posterior median
#    10%, 90% = posterior 10% and 90% quantiles
#    Values are azimuths of the axes (for S1,S2,S3,SHmax) in DEGREES
#           or value of the stress ratio nu = (Smax-Smid)/(Smax-Smin)
# 3. rname.2dparameters.dat
#       11 line comma separated file (1 row header)
#          name, mean, map, michael
#    name = nu, Phi, Theta, Psi,
#           S1:Phi, S1:Theta,  S2:Phi, S2:Theta,  S3:Phi, S3:Theta,
#    mean, map as above
#    michael = equivalent value from Michael 1984,1987 inversion
#    all angles in degrees
# 4. rname.s123grid.dat
#    information on grid used
# 5. rname.s1density.dat, rname.s2density.dat, rname.s3density.dat
#    2D tabulation (on the defined phi/theta grid) of the
#    posterior densities of the three axes (S1,S2,S3)=(Smax,Smid,Smin)
#    NB - these are densities in (phi,cos(theta))
#         and should be scaled by sin(theta) if they are to be integrated
#         over theta
# 6. rname.nudensity.dat
#    tabulation of the posterior density of the stress ratio nu
# 7. rname.shmax.density.dat
#    tabulation of the posterior desnity of the SHmax azimuth

fname <- paste(outputdir,"/",rname,".s123grid.dat",sep="")
# Get vector of azimuths for plotting
tmp <- as.numeric(as.vector(read.table(fname,sep=",",skip=2,nrows=1)[1,]))
phivec <- seq(from=tmp[2],to=tmp[3],length=tmp[1])
# Get vector of dips for plotting
tmp <- as.numeric(as.vector(read.table(fname,sep=",",skip=4,nrows=1)[1,]))
thetavec <- seq(from=tmp[2],to=tmp[3],length=tmp[1])

# Posterior densities of S1 (Smax)
fname <- paste(outputdir,"/",rname,".s1density.dat",sep="")
z <- read.table(fname,sep=",")
z <- as.matrix(z)
# flat contour map of Smin
contour(phivec,thetavec,z,xlab=expression(phi),ylab=expression(theta))
# stereogram of Smin
contour.stereogram(phivec*pi/180,thetavec*pi/180,z)
par(opar)

# Posterior densities of S2 (Smid)
fname <- paste(outputdir,"/",rname,".s2density.dat",sep="")
z <- read.table(fname,sep=",")
z <- as.matrix(z)
# flat contour map Smid
contour(phivec,thetavec,z,xlab=expression(phi),ylab=expression(theta))
# stereogram of Smid
contour.stereogram(phivec*pi/180,thetavec*pi/180,z)
par(opar)

# Posterior densities of S3 (Smin)
fname <- paste(outputdir,"/",rname,".s3density.dat",sep="")
z <- read.table(fname,sep=",")
z <- as.matrix(z)
# flat contour map of Smax
contour(phivec,thetavec,z,xlab=expression(phi),ylab=expression(theta))
# stereogram of Smax
contour.stereogram(phivec*pi/180,thetavec*pi/180,z)
par(opar)

stress.stereonet(stress.out$invresult,
                 stress.out$psummary,
                 newpage=F,
                 draw.grid=T,draw.map.rmat=T,
                 cols=c("red","green","blue"),
                 store.stereonet.file="outdata/junk") 

names(stress.out)
names(stress.out$psummary)

phiv2 <- cbind( c(45,50,55), c(65,65,65) )*pi/180

points.phiv2.stereogram(phiv2, pch=17, col="red", cex=2)
lines.phiv2.stereogram(phiv2, lwd=2, col="blue")

lines.phiv2.stereogram(cbind(c(45,225),c(90,90))*pi/180, lwd=2, col="blue",
                       lty=2)
text.phiv2.stereogram(c(270,45)*pi/180, label="hello")
axis(1); axis(2)

text(-1,-1, label="jjj", cex=2)

# End
#-----------------------------------------------------------------------------



