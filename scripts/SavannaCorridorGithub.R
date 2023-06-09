#####################################################################################################################
## Analysis of palaeoecological records across South-East Asia to determine the evidence for regime shifts between ## 
## open savannas and dense tropical forests occurred since the Last Glacial Maximum
#####################################################################################################################

## Rebecca Hamilton, Corey Bradshaw, Frédérik Saltré
## June 2023

## load libraries
library(spatstat)
library(gstat)
library(maps)
library(sp)
library(ape)
library(permute)
library(ggplot2)
library(dplyr)
library(boot)
library(tmvnsim)
library(wCorr)
library(hrbrthemes)

###########################################################################
## STEP 1 : DATA INPUT
###########################################################################

#POLLEN for analysis
dat1 <- read.csv("./data/G6-4.csv", header=T)
dat2 <- read.csv("./data/SH19014_grass.csv", header=T)
dat3 <- read.csv("./data/G5_6_149P2.csv", header=T)
dat4 <- read.csv("./data/NPK2.csv", header=T)
dat5 <- read.csv("./data/hordorli.csv", header=T)
dat6 <- read.csv("./data/18300.csv", header=T)
dat7 <- read.csv("./data/18323.csv", header=T)
dat8 <- read.csv("./data/18302.csv", header=T)
dat9 <- read.csv("./data/CB19.csv", header=T)
dat10 <- read.csv("./data/MD063075.csv", header=T)
dat11 <- read.csv("./data/NS-0725.csv", header=T)
dat12 <- read.csv("./data/17964.csv", header=T)
dat13 <- read.csv("./data/DDA.csv", header=T)
dat14 <- read.csv("./data/G4_K12P1.csv", header=T)
dat15 <- read.csv("./data/PB-A.csv", header=T)
dat16 <- read.csv("./data/GEOB10069_3.csv", header=T)
dat17 <- read.csv("./data/PSS.csv", header=T)
dat18 <- read.csv("./data/GeoB10053_7.csv", header=T)
dat19 <- read.csv("./data/G5_2_056P.csv", header=T)
dat20 <- read.csv("./data/LL2.csv", header=T)
dat21 <- read.csv("./data/RD-3.csv", header=T)
dat22 <- read.csv("./data/KUM3.csv", header=T)
dat23 <- read.csv("./data/G4_K4P3.csv", header=T)
dat24 <- read.csv("./data/BYK2.csv", header=T)

#d13C for analysis
dat25 <- read.csv("./data/TOW9.csv", header=T)
dat26 <- read.csv("./data/BJ8_03_91GGC.csv", header=T)
dat27 <- read.csv("./data/GEOB10069_3.csv", header=T)
dat28 <- read.csv("./data/MAT10_2B.csv", header=T)
dat29 <- read.csv("./data/SO189_144KL.csv", header=T)
dat30 <- read.csv("./data/MC1.csv", header=T)
dat31 <- read.csv("./data/GeoB10053_7.csv", header=T)
dat32 <- read.csv("./data/mbelen.csv", header=T)

# make variable of interest consistent
dat1$open <- dat1$grassland_drylandpc
dat2$open <- dat2$grass_drylandpc
dat3$open <- dat3$grassland_drylandpc
dat4$open <- dat4$grass_drylandpc
dat5$open <- dat5$grass_drylandpc
dat6$open <- dat6$grass_drylandpc
dat7$open <- dat7$grass_drylandpc
dat8$open <- dat8$grass_drylandpc
dat9$open <- dat9$grass_drylandpc
dat10$open <- dat10$grass_drylandpc
dat11$open <- dat11$grass_drylandpc
dat12$open <- dat12$herbs_pc
dat13$open <- dat13$grass_drylandpc
dat14$open <- dat14$grass_woodland_drylandpc
dat15$open <- dat15$grass_pc
dat16$open <- dat16$C4_pc
dat17$open <- dat17$non_forest_drylandpc
dat18$open <- dat18$grass_pc
dat19$open <- dat19$grass_woodland_pc
dat20$open <- dat20$grass_dryland_pc
dat21$open <- dat21$herb_drylandpc
dat22$open <- dat22$herbs_pc
dat23$open <- dat23$grass_woodland_drylandpc
dat24$open <- dat24$grass_drylandpc
dat25$open <- dat25$d13C_C28
dat26$open <- dat26$d13C_C30_FA
dat27$open <- dat27$d13C_C30FA
dat28$open <- dat28$d13C_C28
dat29$open <- dat29$d13C_C30_AA
dat30$open <- dat30$Guano_d13C
dat31$open <- dat31$C33_d13C
dat32$open <- dat32$Guano_d13C

## visualise data together vs. medial age-depth profile
# pollen
par(mfrow=c(5,5)) # adjust for number of records
plot(dat1$median, dat1$open, type="l", xlab="years BP")
plot(dat2$median, dat2$open, type="l", xlab="years BP")
plot(dat3$median, dat3$open, type="l", xlab="years BP")
plot(dat4$median, dat4$open, type="l", xlab="years BP")
plot(dat5$median, dat5$open, type="l", xlab="years BP")
plot(dat6$median, dat6$open, type="l", xlab="years BP")
plot(dat7$median, dat7$open, type="l", xlab="years BP")
plot(dat8$median, dat8$open, type="l", xlab="years BP")
plot(dat9$median, dat9$open, type="l", xlab="years BP")
plot(dat10$median, dat10$open, type="l", xlab="years BP")
plot(dat11$median, dat11$open, type="l", xlab="years BP")
plot(dat12$median, dat12$open, type="l", xlab="years BP")
plot(dat13$median, dat13$open, type="l", xlab="years BP")
plot(dat14$median, dat14$open, type="l", xlab="years BP")
plot(dat15$median, dat15$open, type="l", xlab="years BP")
plot(dat16$median, dat16$open, type="l", xlab="years BP")
plot(dat17$median, dat17$open, type="l", xlab="years BP")
plot(dat18$median, dat18$open, type="l", xlab="years BP")
plot(dat19$median, dat19$open, type="l", xlab="years BP")
plot(dat20$median, dat20$open, type="l", xlab="years BP")
plot(dat21$median, dat21$open, type="l", xlab="years BP")
plot(dat22$median, dat22$open, type="l", xlab="years BP")
plot(dat23$median, dat23$open, type="l", xlab="years BP")
plot(dat24$median, dat24$open, type="l", xlab="years BP")
par(mfrow=c(1,1))

# d13C
par(mfrow=c(3,4)) # adjust for number of records
plot(dat25$median, dat25$open, type="l", xlab="years BP")
plot(dat26$median, dat26$open, type="l", xlab="years BP")
plot(dat27$median, dat27$open, type="l", xlab="years BP")
plot(dat28$median, dat28$open, type="l", xlab="years BP")
plot(dat29$median, dat29$open, type="l", xlab="years BP")
plot(dat30$median, dat30$open, type="l", xlab="years BP")
plot(dat31$median, dat31$open, type="l", xlab="years BP")
plot(dat32$median, dat32$open, type="l", xlab="years BP")
par(mfrow=c(1,1))

# output file names (add relevant names for new records)
name.out <- c("G6_4_grassland.out","SH19014_grass.out","G5_6_149P2_grassland.out","NPK2_grass.out","hordorli_grass.out",
              "SO18300_grass.out","SO18323_grass.out","SO18302_grass.out","CB19_grass.out","MD063075_grass.out",
              "NS_0725_grass.out","d17964_herbs.out","DDA_grass.out","G4_K12P1_grasswood.out","PB_A_grass.out",
              "GEOB10069_3_C4.out","PSS_nonforest.out","GeoB10053_7_grass.out","G5_2_056P.out","LL2.out",
              "RD_3.out","KUM3.out","G4_K4P3.out","BYK2.out","TOW9.out",
              "BJ8_91GCC.out","GEOB69_3.out","MAT10_2B.out","SO189_144KL.out",
              "MC1.out","GeoB10053_7.out","mbelen.out")


###################################################################################
## STEP 2: DATA STANDARISATION (MAKE AGE INTERVALS, TIMESPAN & SCALING CONSISTENT)
###################################################################################

# 2A - SCALING FOR CORRELATIONS
# number of iterations
iter <- 10000
itdiv <- iter/10

# standardise to set age interval (min, max, interval)
agest <- seq(2000,34000,2000)
dat.vec <- paste(rep("dat",32),1:length(name.out),sep="") # adjust for number of files

for (d in 1:length(dat.vec)) {
  
  dat <- eval(as.name(dat.vec[d]))
  
  valapproxmat <- valscmat <- valrmat <- matrix(data=NA, nrow=iter, ncol=length(agest))
  for (i in 1:iter) {
    
    age.it <- sd <- rep(NA,dim(dat)[1])
    for (t in 1:dim(dat)[1]) {
      sd[t] <- mean(c((dat$max[t] - dat$median[t]), (dat$median[t] - dat$min[t]))) / 1.96
      age.it[t] <- rnorm(1, dat$median[t], sd[t])
    } # end t
    
    val.approx <- approx(age.it, dat$open, xout = agest) # grassland/dryland %
    valapproxmat[i,] <- val.approx$y
    valscmat[i,] <- scale(val.approx$y, center=F, scale=T) # scale #change to center = F for SPAC
    valrmat[i,2:length(val.approx$x)] <- rev(log(rev(valscmat[i,])[2:length(valscmat[i,])] / rev(valscmat[i,])[1:(length(valscmat[i,])-1)])) # instantaneous exponential rate of change; r <- log(Nt1/Nt)
    
    if (i %% itdiv==0) print(i) 
    
  } # end i
  
  # r
  r.md <- apply(valrmat, MARGIN=2, median, na.rm=T)
  r.md <- ifelse(is.infinite(r.md)==T, NA, r.md)
  r.lo <- apply(valrmat, MARGIN=2, quantile, probs=0.025, na.rm=T)
  r.lo <- ifelse(is.infinite(r.lo)==T, NA, r.lo)
  r.up <- apply(valrmat, MARGIN=2, quantile, probs=0.975, na.rm=T)
  r.up <- ifelse(is.infinite(r.up)==T, NA, r.up)
   
  approx.med <- approx(dat$median, dat$open, xout = agest) # grass/dryland %
  approx.scmed <- scale(approx.med$y, center=F, scale=T) # scale
  approx.scr <- rev(log(rev(approx.scmed)[2:length(approx.scmed)] / rev(approx.scmed)[1:(length(approx.scmed)-1)])) # instantaneous exponential rate of change; r <- log(Nt1/Nt)
  
  # scaled values
  sc.md <- apply(valscmat, MARGIN=2, median, na.rm=T)
  sc.lo <- apply(valscmat, MARGIN=2, quantile, probs=0.025, na.rm=T)
  sc.up <- apply(valscmat, MARGIN=2, quantile, probs=0.975, na.rm=T)
  
  # interpolated values
  val.md <- apply(valapproxmat, MARGIN=2, median, na.rm=T)
  val.lo <- apply(valapproxmat, MARGIN=2, quantile, probs=0.025, na.rm=T)
  val.up <- apply(valapproxmat, MARGIN=2, quantile, probs=0.975, na.rm=T)
  
  # save output data.frame
  out <- data.frame(agest, val.md, val.up, val.lo, sc.md, sc.up, sc.lo, r.md, r.up, r.lo)
  assign(name.out[d], out)
  
  print("####################")
  print(name.out[d])
  print("####################")
  
}

#save outputs to environment
save.image("SC_data_34K_at_2K.RData")

#export data files as .csv
write.csv(G6_4_grassland.out, "G6_4_grasslandOut.csv")
write.csv(SH19014_grass.out,"SH19014_grassOut.csv")
write.csv(G5_6_149P2_grassland.out,"G5_6_149P2_grasslandOut.csv")
write.csv(NPK2_grass.out, "NPK2_grassOUT.csv")
write.csv(hordorli_grass.out, "hordorli_grassOut.csv")
write.csv(SO18300_grass.out,"SO18300_grassOut.csv")
write.csv(SO18323_grass.out, "SO18323_grassOut.csv")
write.csv(SO18302_grass.out,"SO18302_grassOut.csv")
write.csv(CB19_grass.out,"CB19_grassOut.csv")
write.csv(MD063075_grass.out,"MD063075_grassOut.csv")
write.csv(NS_0725_grass.out,"NS_0725_grassOut.csv")
write.csv(d17964_herbs.out,"d17964_herbsOut.csv")
write.csv(DDA_grass.out,"DDA_grassOut.csv")
write.csv(G4_K12P1_grasswood.out,"G4_K12P1_grasswoodOut.csv")
write.csv(PB_A_grass.out,"PB_A_grassOut.csv")
write.csv(GEOB10069_3_C4.out,"GEOB10069_3_C4Out.csv")
write.csv(PSS_nonforest.out,"PSS_nonforestOut.csv")
write.csv(GeoB10053_7_grass.out,"GeoB10053_7_grassOut.csv")
write.csv(G5_2_056P.out,"G5_2_056POut.csv")
write.csv(LL2.out,"LL2Out.csv")
write.csv(RD_3.out,"RD_3.out.csv")
write.csv(KUM3.out,"KUM3Out.csv")
write.csv(G4_K4P3.out,"G4_K4P3Out.csv")
write.csv(BYK2.out,"BYK2out.csv")
write.csv(TOW9.out,"TOW9out.csv")
write.csv(BJ8_91GCC.out, "BJ8_91GCC.out.csv")
write.csv(GEOB69_3.out,"GEOB69_3out.csv")
write.csv(MAT10_2B.out,"MAT10_2Bout.csv")
write.csv(SO189_144KL.out,"SO189_144KLout.csv")
write.csv(MC1.out ,"MC1out.csv")
write.csv(GeoB10053_7.out, "GeoB10053_7out.csv")
write.csv(mbelen.out,"mbelenout.csv")


###############################################################################################
## STEP 3: SPATIAL AUTOCORRELATION
### this correlates data (step 3b) while accounting for spatial autocorrelation bias (step 3a)
###############################################################################################

###################################################################################
## STEP 3A: autocorrelation - testing for spatial autocorrelations using Moran's I
###################################################################################

# need .csv file with x,y coordinates & scaled, non centred median values. This is built from exported files above.
# Easiest to order (e.g. dat 1 to 32 as above) for ease of use

## NOTE 1. D13C are all negative ##
## NOTE 2. SC_SACor_test.csv built for tested pollen data only; SC_all_SACor.csv built with used pollen and d13C data

## MORAN's I

# create subset of data to analyse
# these are all pollen samples that include at least 1 x LGM sample, and do not have giant data gaps

spatdat <- read.table("./data/SC_all_SACor.csv", header=T,sep=",") 
spatdat.size <- dim(spatdat)
nld <- 10000

select.sub <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
spatdatupl <- spatdat[select.sub ,]

dat.use <- spatdatupl

spatdat.size <- dim(dat.use)

spatdat.dists <- as.matrix(dist(cbind(dat.use$x, dat.use$y)))
spatdat.dists.inv <- 1/spatdat.dists
diag(spatdat.dists.inv) <- 0

spatdatSpat <- rep(NA, spatdat.size[2]-3,1)
cpt <- 0

for (k in 4:spatdat.size[2]) {
  cpt <- cpt+1
  print((cpt/(spatdat.size[2]-3)*100))
  spatdatvar <- dat.use[,k]
  mor_obs <- Moran.I(spatdatvar, spatdat.dists.inv, na.rm = TRUE)
  out <- rep(NA, nld, 1)
  
  for (i in 1:nld) {
    mat <- spatdatvar[shuffle(spatdat.size[1])]
    MorI <- Moran.I(mat, spatdat.dists.inv,na.rm = TRUE)
    out[i] <- MorI$observed
  }
  
  Moran_out <- as.data.frame(out)

  nbout <- which(abs(Moran_out) >= abs(mor_obs$observed))
  spatdatSpat[cpt] <- length(nbout)/nld
  rm(nbout, Moran_out, MorI,mat)
}

plot(spatdatSpat,type = "l")
write.csv(spatdatSpat, "data_pollenall_Spat.csv")


# these are all samples (pollen and d13C) that include at least 1 x LGM sample, and do not have giant data gaps
select.sub2 <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,25,26,27,28,29,30,31,32)
spatdatupl2 <- spatdat[select.sub2, ]

dat.use <- spatdatupl2

spatdat.size <- dim(dat.use)

spatdat.dists <- as.matrix(dist(cbind(dat.use$x, dat.use$y)))
spatdat.dists.inv <- 1/spatdat.dists
diag(spatdat.dists.inv) <- 0

spatdatSpat2 <- rep(NA, spatdat.size[2]-3,1)
cpt <- 0

for (k in 4:spatdat.size[2]) {
  cpt <- cpt+1
  print((cpt/(spatdat.size[2]-3)*100))
  spatdatvar <- dat.use[,k]
  mor_obs <- Moran.I(spatdatvar, spatdat.dists.inv, na.rm = TRUE)
  out <- rep(NA, nld, 1)
  
  for (i in 1:nld) {
    mat <- spatdatvar[shuffle(spatdat.size[1])]
    MorI <- Moran.I(mat, spatdat.dists.inv,na.rm = TRUE)
    out[i] <- MorI$observed
  }
  
  Moran_out <- as.data.frame(out)
  
  nbout <- which(abs(Moran_out) >= abs(mor_obs$observed))
  spatdatSpat2[cpt] <- length(nbout)/nld
  rm(nbout, Moran_out, MorI,mat)
}

plot(spatdatSpat2,type = "l")
write.csv(spatdatSpat2, "data_all_Spat.csv")


##############################################################################
### STEP3B - correlate data
##############################################################################

## all datasets
out1 <- get(name.out[1])
out2 <- get(name.out[2])
out3 <- get(name.out[3])
out4 <- get(name.out[4])
out5 <- get(name.out[5])
out6 <- get(name.out[6])
out8 <- get(name.out[7])
out8 <- get(name.out[8])
out9 <- get(name.out[9])
out10 <- get(name.out[10])
out11 <- get(name.out[11])
out12 <- get(name.out[12])
out13 <- get(name.out[13])
out14 <- get(name.out[14])
out15 <- get(name.out[15])
out16 <- get(name.out[16])
out17 <- get(name.out[17])
out18 <- get(name.out[18])
out25 <- get(name.out[25])
out26 <- get(name.out[26])
out27 <- get(name.out[27])
out28 <- get(name.out[28])
out29 <- get(name.out[29])
out30 <- get(name.out[30])
out31 <- get(name.out[31])
out32 <- get(name.out[32])


## import spatial autocorrelation probabilities
allSpatP <- read.csv("data_all_Spat.csv")

## CORRELATE DATASETS using all Spat data
## the following code is an example of a bivariate spatially constrained correlation between
## proxies 8 and 32; the full correlation matrix values are presented in the following section

par(mfrow=c(2,1))
plot(out8$agest, out8$val.md, type="l", ylim=c(min(out8$val.lo,na.rm=T),max(out8$val.up,na.rm=T)))
lines(out8$agest, out8$val.lo, lty=2, col="red")
lines(out8$agest, out8$val.up, lty=2, col="red")
plot(out32$agest, out32$val.md, type="l", ylim=c(min(out32$val.lo,na.rm=T),max(out32$val.up,na.rm=T)))
lines(out32$agest, out32$val.lo, lty=2, col="red")
lines(out32$agest, out32$val.up, lty=2, col="red")
par(mfrow=c(1,1))

# rescale
out8scale <- scale(out8$val.md, center=T, scale=T)
out8loscale <- scale(out8$val.lo, center=attr(out8scale, "scaled:center"), scale=attr(out8scale, "scaled:scale"))
out8upscale <- scale(out8$val.up, center=attr(out8scale, "scaled:center"), scale=attr(out8scale, "scaled:scale"))
plot(out8$agest, out8scale, type="l", ylim=c(min(out8loscale,na.rm=T),max(out8upscale,na.rm=T)))
lines(out8$agest, out8loscale, lty=2, col="red")
lines(out8$agest, out8upscale, lty=2, col="red")

out32scale <- scale(out32$val.md, center=T, scale=T)
out32loscale <- scale(out32$val.lo, center=attr(out32scale, "scaled:center"), scale=attr(out32scale, "scaled:scale"))
out32upscale <- scale(out32$val.up, center=attr(out32scale, "scaled:center"), scale=attr(out32scale, "scaled:scale"))
plot(out32$agest, out32scale, type="l", ylim=c(min(out32loscale,na.rm=T),max(out32upscale,na.rm=T)))
lines(out32$agest, out32loscale, lty=2, col="red")
lines(out32$agest, out32upscale, lty=2, col="red")

# stochastic correlation (bivariate)
iter6 <- 10000
itdiv6 <- iter6/10

# change these accordingly (right-hand side refers to number for specific dataset)
input1lo <- out32loscale
input1up <- out32upscale
input2lo <- out8loscale
input2up <- out8upscale

cor.vecN <- cor.vecSA <- rep(NA,iter6)
for (c in 1:iter6) {
  out1it.val <- out2it.val <- rep(NA, length(out1$agest))
  for (t in 1:length(out1$agest)) {
    # resample for each age from each series
    out1it.val[t] <- runif(1,min=input1lo[t], max=input1up[t])
    out2it.val[t] <- runif(1,min=input2lo[t], max=input2up[t])
  } # end t
  cor.vecN[c] <- cor(out1it.val, out2it.val, use="na.or.complete", method="spearman")
  # combine into data.frame to remove NAs
  it.dat <- na.omit(data.frame(out1it.val,out2it.val,allSpatP$x))
  cor.vecSA[c] <- weightedCorr(x=it.dat[,1], y=it.dat[,2], method="Spearman",
                               weights = it.dat[,3])
  if (c %% itdiv6==0) print(c) 
  
} # end c

# correlation test with spatial autocorrelation
corSA.md <- median(cor.vecSA, na.rm=T)
corSA.lo <- quantile(cor.vecSA, probs=0.025, na.rm=T)
corSA.up <- quantile(cor.vecSA, probs=0.975, na.rm=T)
print(c(corSA.lo, corSA.up))

# correlation test without spatial autocorrelation
corN.md <- median(cor.vecN, na.rm=T)
corN.lo <- quantile(cor.vecN, probs=0.025, na.rm=T)
corN.up <- quantile(cor.vecN, probs=0.975, na.rm=T)
print(c(corN.lo, corN.up))


##############################################################################
## STEP 4 - CORRELATION HEAT MAP
##############################################################################
## this step generates the correlation matrix shown in Figure 5 in the manuscript
# import data
samp1 <- read.csv("./data/cor_matrix_allmean.csv")
dsNam.vec <- samp1$X
samp1.exp <- expand.grid(X=dsNam.vec, Y=dsNam.vec)
samp1.exp$cor <- 0

for (i in 1:dim(samp1.exp)[1]) {
  x.var <- as.character(samp1.exp[i,1])
  y.var <- as.character(samp1.exp[i,2])
  samp1.exp$cor[i] <- samp1[which(samp1$X == x.var), which(attr(samp1[1,], "name") == y.var)]
}
ggplot(samp1.exp, aes(X, Y, fill= cor)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", limits=c(-1,1)) +
  labs(x = "", y = "")



##############################################################################
## STEP 5 - VOLATILITY
##############################################################################
## this step generates the combined time series limits shown in Figure 6 in the manuscript
## Volatility using interpolated values scaled to a single series (standardisation)

# Datasets1 = all non-lowland records
val.md.comb <- c(NPK2_grass.out$val.md,
               hordorli_grass.out$val.md,
               DDA_grass.out$val.md,
               PB_A_grass.out$val.md,
               PSS_nonforest.out$val.md,
               TOW9.out$val.md,
               MAT10_2B.out$val.md,
               mbelen.out$val.md)

val.md.comb <- data.frame(NPK2_grass.out$val.md,
                          hordorli_grass.out$val.md,
                          DDA_grass.out$val.md,
                          PB_A_grass.out$val.md,
                          PSS_nonforest.out$val.md,
                          TOW9.out$val.md,
                          MAT10_2B.out$val.md,
                          mbelen.out$val.md)

# scale first
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb[,4] <- scale(val.md.comb[,4], center=T, scale=T)
valsc.md.comb[,5] <- scale(val.md.comb[,5], center=T, scale=T)
valsc.md.comb[,6] <- scale(val.md.comb[,6], center=T, scale=T)
valsc.md.comb[,7] <- scale(val.md.comb[,7], center=T, scale=T)
valsc.md.comb[,8] <- scale(val.md.comb[,8], center=T, scale=T)
valsc.md.comb

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
lines(agest,valsc.md.comb[,4],lty=4)
lines(agest,valsc.md.comb[,5],lty=5)
lines(agest,valsc.md.comb[,6],lty=6)
lines(agest,valsc.md.comb[,7],lty=7)
lines(agest,valsc.md.comb[,8],lty=8)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(NPK2_grass.out$sc.lo,
                         hordorli_grass.out$sc.lo,
                         DDA_grass.out$sc.lo,
                         PB_A_grass.out$sc.lo,
                         PSS_nonforest.out$sc.lo,
                         TOW9.out$sc.lo,
                         MAT10_2B.out$sc.lo,
                         mbelen.out$sc.lo)

sc.up.comb <- data.frame(NPK2_grass.out$sc.up,
                         hordorli_grass.out$sc.up,
                         DDA_grass.out$sc.up,
                         PB_A_grass.out$sc.up,
                         PSS_nonforest.out$sc.up,
                         TOW9.out$sc.up,
                         MAT10_2B.out$sc.up,
                         mbelen.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")


# Datasets2 = all non-lowland pollen
val.md.comb <- c(NPK2_grass.out$val.md,
                 hordorli_grass.out$val.md,
                 DDA_grass.out$val.md,
                 PB_A_grass.out$val.md,
                 PSS_nonforest.out$val.md)

val.md.comb <- data.frame(NPK2_grass.out$val.md,
                          hordorli_grass.out$val.md,
                          DDA_grass.out$val.md,
                          PB_A_grass.out$val.md,
                          PSS_nonforest.out$val.md)

# scale first
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb[,4] <- scale(val.md.comb[,4], center=T, scale=T)
valsc.md.comb[,5] <- scale(val.md.comb[,5], center=T, scale=T)
valsc.md.comb

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
lines(agest,valsc.md.comb[,4],lty=4)
lines(agest,valsc.md.comb[,5],lty=5)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=31400, lty=2, col="red")
abline(v=15500, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(NPK2_grass.out$sc.lo,
                         hordorli_grass.out$sc.lo,
                         DDA_grass.out$sc.lo,
                         PB_A_grass.out$sc.lo,
                         PSS_nonforest.out$sc.lo)

sc.up.comb <- data.frame(NPK2_grass.out$sc.up,
                         hordorli_grass.out$sc.up,
                         DDA_grass.out$sc.up,
                         PB_A_grass.out$sc.up,
                         PSS_nonforest.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")


# Datasets3 = non-lowland d13C
val.md.comb <- c(TOW9.out$val.md,
                 MAT10_2B.out$val.md,
                 mbelen.out$val.md)

val.md.comb <- data.frame(TOW9.out$val.md,
                          MAT10_2B.out$val.md,
                          mbelen.out$val.md)

# scale first
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(TOW9.out$sc.lo,
                         MAT10_2B.out$sc.lo,
                         mbelen.out$sc.lo)

sc.up.comb <- data.frame(TOW9.out$sc.up,
                         MAT10_2B.out$sc.up,
                         mbelen.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")


# Datasets4 = all lowland records
val.md.comb <- c(G6_4_grassland.out$val.md,
                SH19014_grass.out$val.md,
                G5_6_149P2_grassland.out$val.md,
                SO18300_grass.out$val.md,
                SO18323_grass.out$val.md,
                SO18302_grass.out$val.md,
                CB19_grass.out$val.md,
                MD063075_grass.out$val.md,
                NS_0725_grass.out$val.md,
                d17964_herbs.out$val.md,
                G4_K12P1_grasswood.out$val.md,
                GEOB10069_3_C4.out$val.md,
                GeoB10053_7_grass.out$val.md,
                BJ8_91GCC.out$val.md,
                GEOB69_3.out$val.md,
                SO189_144KL.out$val.md,
                MC1.out$val.md,
                GeoB10053_7.out$val.md)

val.md.comb <- data.frame(G6_4_grassland.out$val.md,
                 SH19014_grass.out$val.md,
                 G5_6_149P2_grassland.out$val.md,
                 SO18300_grass.out$val.md,
                 SO18323_grass.out$val.md,
                 SO18302_grass.out$val.md,
                 CB19_grass.out$val.md,
                 MD063075_grass.out$val.md,
                 NS_0725_grass.out$val.md,
                 d17964_herbs.out$val.md,
                 G4_K12P1_grasswood.out$val.md,
                 GEOB10069_3_C4.out$val.md,
                 GeoB10053_7_grass.out$val.md,
                 BJ8_91GCC.out$val.md,
                 GEOB69_3.out$val.md,
                 SO189_144KL.out$val.md,
                 MC1.out$val.md,
                 GeoB10053_7.out$val.md)

# scale first
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb[,4] <- scale(val.md.comb[,4], center=T, scale=T)
valsc.md.comb[,5] <- scale(val.md.comb[,5], center=T, scale=T)
valsc.md.comb[,6] <- scale(val.md.comb[,6], center=T, scale=T)
valsc.md.comb[,7] <- scale(val.md.comb[,7], center=T, scale=T)
valsc.md.comb[,8] <- scale(val.md.comb[,8], center=T, scale=T)
valsc.md.comb[,9] <- scale(val.md.comb[,9], center=T, scale=T)
valsc.md.comb[,10] <- scale(val.md.comb[,10], center=T, scale=T)
valsc.md.comb[,11] <- scale(val.md.comb[,11], center=T, scale=T)
valsc.md.comb[,12] <- scale(val.md.comb[,12], center=T, scale=T)
valsc.md.comb[,13] <- scale(val.md.comb[,13], center=T, scale=T)
valsc.md.comb[,14] <- scale(val.md.comb[,14], center=T, scale=T)
valsc.md.comb[,15] <- scale(val.md.comb[,15], center=T, scale=T)
valsc.md.comb[,16] <- scale(val.md.comb[,16], center=T, scale=T)
valsc.md.comb[,17] <- scale(val.md.comb[,17], center=T, scale=T)
valsc.md.comb[,18] <- scale(val.md.comb[,18], center=T, scale=T)

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
lines(agest,valsc.md.comb[,4],lty=2)
lines(agest,valsc.md.comb[,5],lty=3)
lines(agest,valsc.md.comb[,6],lty=2)
lines(agest,valsc.md.comb[,7],lty=3)
lines(agest,valsc.md.comb[,8],lty=2)
lines(agest,valsc.md.comb[,9],lty=3)
lines(agest,valsc.md.comb[,10],lty=2)
lines(agest,valsc.md.comb[,11],lty=3)
lines(agest,valsc.md.comb[,12],lty=2)
lines(agest,valsc.md.comb[,13],lty=3)
lines(agest,valsc.md.comb[,14],lty=2)
lines(agest,valsc.md.comb[,15],lty=3)
lines(agest,valsc.md.comb[,16],lty=2)
lines(agest,valsc.md.comb[,17],lty=3)
lines(agest,valsc.md.comb[,18],lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(G6_4_grassland.out$sc.lo,
                        SH19014_grass.out$sc.lo,
                        G5_6_149P2_grassland.out$sc.lo,
                        SO18300_grass.out$sc.lo,
                        SO18323_grass.out$sc.lo,
                        SO18302_grass.out$sc.lo,
                        CB19_grass.out$sc.lo,
                        MD063075_grass.out$sc.lo,
                        NS_0725_grass.out$sc.lo,
                        d17964_herbs.out$sc.lo,
                        G4_K12P1_grasswood.out$sc.lo,
                        GEOB10069_3_C4.out$sc.lo,
                        GeoB10053_7_grass.out$sc.lo,
                        BJ8_91GCC.out$sc.lo,
                        GEOB69_3.out$sc.lo,
                        SO189_144KL.out$sc.lo,
                        MC1.out$sc.lo,
                        GeoB10053_7.out$sc.lo)

sc.up.comb <- data.frame(G6_4_grassland.out$sc.up,
                         SH19014_grass.out$sc.up,
                         G5_6_149P2_grassland.out$sc.up,
                         SO18300_grass.out$sc.up,
                         SO18323_grass.out$sc.up,
                         SO18302_grass.out$sc.up,
                         CB19_grass.out$sc.up,
                         MD063075_grass.out$sc.up,
                         NS_0725_grass.out$sc.up,
                         d17964_herbs.out$sc.up,
                         G4_K12P1_grasswood.out$sc.up,
                         GEOB10069_3_C4.out$sc.up,
                         GeoB10053_7_grass.out$sc.up,
                         BJ8_91GCC.out$sc.up,
                         GEOB69_3.out$sc.up,
                         SO189_144KL.out$sc.up,
                         MC1.out$sc.up,
                         GeoB10053_7.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")


# Datasets5 = lowland pollen records
val.md.comb <- c(G6_4_grassland.out$val.md,
                 SH19014_grass.out$val.md,
                 G5_6_149P2_grassland.out$val.md,
                 SO18300_grass.out$val.md,
                 SO18323_grass.out$val.md,
                 SO18302_grass.out$val.md,
                 CB19_grass.out$val.md,
                 MD063075_grass.out$val.md,
                 NS_0725_grass.out$val.md,
                 d17964_herbs.out$val.md,
                 G4_K12P1_grasswood.out$val.md,
                 GEOB10069_3_C4.out$val.md,
                 GeoB10053_7_grass.out$val.md)

val.md.comb <- data.frame(G6_4_grassland.out$val.md,
                          SH19014_grass.out$val.md,
                          G5_6_149P2_grassland.out$val.md,
                          SO18300_grass.out$val.md,
                          SO18323_grass.out$val.md,
                          SO18302_grass.out$val.md,
                          CB19_grass.out$val.md,
                          MD063075_grass.out$val.md,
                          NS_0725_grass.out$val.md,
                          d17964_herbs.out$val.md,
                          G4_K12P1_grasswood.out$val.md,
                          GEOB10069_3_C4.out$val.md,
                          GeoB10053_7_grass.out$val.md)

# scale first because we are using values
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb[,4] <- scale(val.md.comb[,4], center=T, scale=T)
valsc.md.comb[,5] <- scale(val.md.comb[,5], center=T, scale=T)
valsc.md.comb[,6] <- scale(val.md.comb[,6], center=T, scale=T)
valsc.md.comb[,7] <- scale(val.md.comb[,7], center=T, scale=T)
valsc.md.comb[,8] <- scale(val.md.comb[,8], center=T, scale=T)
valsc.md.comb[,9] <- scale(val.md.comb[,9], center=T, scale=T)
valsc.md.comb[,10] <- scale(val.md.comb[,10], center=T, scale=T)
valsc.md.comb[,11] <- scale(val.md.comb[,11], center=T, scale=T)
valsc.md.comb[,12] <- scale(val.md.comb[,12], center=T, scale=T)
valsc.md.comb[,13] <- scale(val.md.comb[,13], center=T, scale=T)

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
lines(agest,valsc.md.comb[,4],lty=2)
lines(agest,valsc.md.comb[,5],lty=3)
lines(agest,valsc.md.comb[,6],lty=2)
lines(agest,valsc.md.comb[,7],lty=3)
lines(agest,valsc.md.comb[,8],lty=2)
lines(agest,valsc.md.comb[,9],lty=3)
lines(agest,valsc.md.comb[,10],lty=2)
lines(agest,valsc.md.comb[,11],lty=3)
lines(agest,valsc.md.comb[,12],lty=2)
lines(agest,valsc.md.comb[,13],lty=3)

abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(G6_4_grassland.out$sc.lo,
                         SH19014_grass.out$sc.lo,
                         G5_6_149P2_grassland.out$sc.lo,
                         SO18300_grass.out$sc.lo,
                         SO18323_grass.out$sc.lo,
                         SO18302_grass.out$sc.lo,
                         CB19_grass.out$sc.lo,
                         MD063075_grass.out$sc.lo,
                         NS_0725_grass.out$sc.lo,
                         d17964_herbs.out$sc.lo,
                         G4_K12P1_grasswood.out$sc.lo,
                         GEOB10069_3_C4.out$sc.lo,
                         GeoB10053_7_grass.out$sc.lo)

sc.up.comb <- data.frame(G6_4_grassland.out$sc.up,
                         SH19014_grass.out$sc.up,
                         G5_6_149P2_grassland.out$sc.up,
                         SO18300_grass.out$sc.up,
                         SO18323_grass.out$sc.up,
                         SO18302_grass.out$sc.up,
                         CB19_grass.out$sc.up,
                         MD063075_grass.out$sc.up,
                         NS_0725_grass.out$sc.up,
                         d17964_herbs.out$sc.up,
                         G4_K12P1_grasswood.out$sc.up,
                         GEOB10069_3_C4.out$sc.up,
                         GeoB10053_7_grass.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")


# Datasets6 = lowland d13C records
val.md.comb <- c(BJ8_91GCC.out$val.md,
                 GEOB69_3.out$val.md,
                 SO189_144KL.out$val.md,
                 MC1.out$val.md,
                 GeoB10053_7.out$val.md)

val.md.comb <- data.frame(BJ8_91GCC.out$val.md,
                          GEOB69_3.out$val.md,
                          SO189_144KL.out$val.md,
                          MC1.out$val.md,
                          GeoB10053_7.out$val.md)

# scale first
valsc.md.comb <- val.md.comb
valsc.md.comb[,1] <- scale(val.md.comb[,1], center=T, scale=T)
valsc.md.comb[,2] <- scale(val.md.comb[,2], center=T, scale=T)
valsc.md.comb[,3] <- scale(val.md.comb[,3], center=T, scale=T)
valsc.md.comb[,4] <- scale(val.md.comb[,4], center=T, scale=T)
valsc.md.comb[,5] <- scale(val.md.comb[,5], center=T, scale=T)

# visualise to check
plot(agest,valsc.md.comb[,1], type="l", xlab="YBP", ylab="scaled value", ylim=c(min(valsc.md.comb, na.rm=T),max(valsc.md.comb, na.rm=T)))
lines(agest,valsc.md.comb[,2],lty=2)
lines(agest,valsc.md.comb[,3],lty=3)
lines(agest,valsc.md.comb[,4],lty=2)
lines(agest,valsc.md.comb[,5],lty=3)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")

# upper/lower (scaled)
sc.lo.comb <- data.frame(BJ8_91GCC.out$sc.lo,
                         GEOB69_3.out$sc.lo,
                         SO189_144KL.out$sc.lo,
                         MC1.out$sc.lo,
                         GeoB10053_7.out$sc.lo)

sc.up.comb <- data.frame(BJ8_91GCC.out$sc.up,
                         GEOB69_3.out$sc.up,
                         SO189_144KL.out$sc.up,
                         MC1.out$sc.up,
                         GeoB10053_7.out$sc.up)

iter6 <- 1000
itdiv6 <- iter6/10
trend.mean.mat <- matrix(data=NA, nrow=iter6, ncol=length(agest))

for (i in 1:iter6) {
  comb.ran.it <- matrix(data=NA,nrow=length(agest), ncol=dim(sc.up.comb)[2])
  for (t in 1:length(agest)) {
    for (r in 1:dim(sc.up.comb)[2]) {
      comb.ran.it[t,r] <- runif(1, min=sc.lo.comb[t,r], max=sc.up.comb[t,r])
    } # end r
  } # end t
  trend.mean.mat[i,] <- apply(comb.ran.it, MARGIN=1, median, na.rm=T)
  if (i %% itdiv6==0) print(i) 
}
trend.mean.lo <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
trend.mean.up <- apply(trend.mean.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
plot(agest, trend.mean.lo, type="l", lty=2, ylab="mean scaled val", ylim=c(min(trend.mean.lo,na.rm=T), max(trend.mean.up, na.rm=T)))
lines(agest, trend.mean.up, lty=2)
abline(h=0, lty=3, lwd=3, col="blue")
abline(v=29000, lty=2, col="red")
abline(v=12000, lty=2, col="red")
