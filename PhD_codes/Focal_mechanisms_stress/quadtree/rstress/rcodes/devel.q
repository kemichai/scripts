#############################################################################
# Routines in development
#############################################################################
source("first.q")
reload()
par(opar)

# Want to calculate Smax in various planes
# Shmax is in the horizontal plane

outputdir <- "outdata"
rname <- "dixie"
load(paste(outputdir,rname,".out",sep=""))
# loads the object "stress.out"
length(stress.out)
names(stress.out)
stress.out$fname
stress.out$sdr
names(stress.out$invresult)
stress.out$invresult$nepar
stress.out$invresult$epar.vecs
dim(stress.out$invresult$qvals0) # 11 51 61 21 = nu phi theta psi

ls()
length(phi.shmax.f)

#############################################################################
