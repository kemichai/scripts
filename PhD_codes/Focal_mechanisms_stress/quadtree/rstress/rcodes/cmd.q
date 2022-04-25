#############################################################################
reload()
#-----------------------------------------------------------------------------
# For Zara - 4 July 2011
sdr <- c(268,50,290)*pi/180
rmat <- rotmat.sdr(sdr)
sdr1 <- extract.sdr.rotmat(rmat)
phiv1 <- convert.sdr.phiv(sdr1)
180/pi*sdr1

sdr2 <- flip.sdr(sdr)
phiv2 <- convert.sdr.phiv(sdr2)
sdr2*180/pi
##[1]  58.47985  43.95821 247.82401

start.stereogram()
points.phiv.stereogram(phiv1, revvecs=T, clip=F)

start.stereogram()
points.phiv.stereogram(phiv2, revvecs=T, clip=F)

points.phiv.stereogram(rbind(phiv1,phiv2), revvecs=T, clip=F)

start.stereogram()
points.phiv.stereogram(rbind(phiv1,phiv2), revvecs=T, clip=T)

start.stereogram()
points.rmat.stereogram(rmat, revvecs=T, clip=T)

#############################################################################

#############################################################################
