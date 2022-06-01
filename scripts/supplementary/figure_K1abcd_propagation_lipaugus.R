# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARY
#
#=======================================================================================================#
if (!require("signal")) install.packages("signal")
if (!require("oce")) install.packages("oce")
if (!require("tuneR")) install.packages("tuneR")
if (!require("seewave")) install.packages("seewave", repos="http://cran.at.r-project.org/")
library(oce, warn.conflicts = F, quietly = T) # image plotting functions and nice color maps
library(signal, warn.conflicts = F, quietly = T) # signal processing functions
library(tuneR, warn.conflicts = F, quietly = T) 
library(seewave) 

#=======================================================================================================#
#
#                                           SET WORKING DIRECTORIES 
#                                           LOAD MY TOOLBOX
#
#=======================================================================================================#

##### Change working directory to the current script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### CALL functions from MY TOOLBOX
source("./toolbox/toolbox_propa.R", chdir=T)

#=================================================================================================
# simulation attenuation on real data : ex Lipaugus vociferans
############## Lipaugus vociferans
filename = "../data/guiana/LIPAUGUS/S4A09154_20191120_134500.wav"

R0 = 200 # reference distance
R = c(200, 300, 400, 600, 800)  # final distances

# extract signals (lipaugus and background)
piha = readWave(filename, from=4.3, to=6.4,units = 'seconds')
bkg  = readWave(filename, from=1.8, to=3.9,units = 'seconds')

# Filtered the signal around the vocal signature of the lipaugus
piha = bwfilter(piha, n=6, from=1000, to= 5500,  output='Wave')

# get the sampling frequency and the waves
fs      = piha@samp.rate
p0.piha = piha@left - mean(piha@left)
p0.bkg  = bkg@left  - mean(bkg@left)

fig = list()
for (ii in 1:length(R))
{
  sprintf("final distance %d m", R[ii])
  # apply attenuation. Environmental values are the same as in the article
  p.piha.att= propa.apply.att(p0.piha, fs, r=R[ii], r0= R0, t=24, rh=87, pa=101325, a0=0.019)
  # mixup with original background with signal after propagation because attenuation process also attenuate background...
  p.piha.att = p.piha.att + p0.bkg*1
  N = 512
  res = specgram(x=p.piha.att,n=N,Fs=fs,window=N,overlap=N/2)
  # frequency vector
  f = res$f
  # PSD
  P = (abs(res$S))
  # convert to dB
  L = psd2dBSPL(P, gain=30, sensitivity=-35)
  # config time axis
  t = res$t
  # plot spectrogram
  image = imagep(x = t,
         y = res$f,
         z = t(L),
         col = oce.colorsViridis,
         ylim = c(0,11000),
         ylab = 'Frequency [Hz]',
         xlab = 'Time [s]',
         zlim = c(20,50),
         drawPalette = T,
         decimate = T,
         useRaster = T)
  fig[[ii]] = list(image)
}


