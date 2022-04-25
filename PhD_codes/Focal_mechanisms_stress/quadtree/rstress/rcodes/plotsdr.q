#############################################################################
# Plotting sdr and means on a stereonet
#-----------------------------------------------------------------------------
# processing a set of .csv files with strike/dip/rake/error
# build and load routines
source("first.q")
reload()

inputdir <- "indata/"
fname <- "sdr_sigma.out"
sdr.df <- read.table(paste(inputdir,fname,sep=""),sep=",",header=F,as.is=T)
names(sdr.df) <- c("strike","dip","rake","err","weight")
sdr.df

sdr <- as.matrix(sdr.df[,c("strike","dip","rake")])
sdr <- pi/180*sdr
msdr <- mean.sdr(sdr)

par(mfrow=c(6,6))
apply(sdr, 1, beach.ball, col="grey")

start.stereogram()
beach.ball(msdr,col="red",add=T)
points.sdr.stereogram(sdr,
                      pch=c("u",NA,"n"),
                      col=c("red","green","blue"))
points.sdr.stereogram(msdr,
                      pch=c(16,NA,16),
                      col=c("red","green","blue"), clip=T)
par(mfrow=c(1,1))


# Note - there is a cleverer beach.ball() function SDRfoc in the
# R package RFOC (needs a couple of other packages too)
par(mfrow=c(1,2))
beach.ball(sdr[1,])
require(RFOC)
M <- SDRfoc(180/pi*sdr[1,1], 180/pi*sdr[1,2], 180/pi*sdr[1,3], 
            u = FALSE, ALIM = c(-1, -1, +1, +1), PLOT=TRUE)
par(mfrow=c(1,1))

# Average a small subset
idx <- 1:1
sdr1 <- sdr[idx,,drop=F]
reload()
msdr1 <- mean.sdr(sdr1,plot=T)

par(mfrow=c(1,1))
start.stereogram()
beach.ball(sdr1,col="grey",add=T)
beach.ball(msdr1,col="black",add=T)
points.sdr.stereogram(sdr1,
                      pch=c("u",NA,"n"),
                      col=c("red","green","blue"))
points.sdr.stereogram(msdr1,
                      pch=c(16,NA,16),
                      col=c("red","green","blue"), clip=T)



idx <- 5+1:2
sdr1 <- sdr[idx,,drop=F]
reload()
msdr1 <- mean.sdr(sdr1,plot=T)

par(mfrow=c(1,1))
start.stereogram()
apply(sdr1,1,beach.ball,col="grey",add=T)
beach.ball(msdr1,col="black",add=T)
points.sdr.stereogram(sdr1,
                      pch=c("u",NA,"n"),
                      col=c("red","green","blue"))
points.sdr.stereogram(msdr1,
                      pch=c(16,NA,16),
                      col=c("red","green","blue"), clip=T)





#############################################################################
