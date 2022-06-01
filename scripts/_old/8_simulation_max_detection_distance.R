# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARY
#
#=======================================================================================================#
if (!require("plotly")) install.packages("plotly")
library(plotly)  

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

#=======================================================================================================#
#
#                                           SET GLOBAL FUNCTIONS
#
#=======================================================================================================#

propa.plot.active_distance <- function (f, d, title='Detection distance', xtitle = 'Distances [m]', ytitle='Frequency [kHz]', name=NULL)
{
  if (is.null(name)) { SHOWLEGEND = FALSE}
  else {SHOWLEGEND = TRUE}
  
  p = plot_ly(x=d, y=f, type = 'bar', orientation = 'h', name=name, showlegend=SHOWLEGEND)
  p = layout(p, 
             title = title,
             paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
             legend = list(orientation = 'v', 
                           automargin = TRUE,
                           font = list(size = 10)),
             xaxis = list(title =xtitle,
                          gridcolor = 'rgb(255,255,255)',
                          showgrid = TRUE,
                          showline = FALSE,
                          showticklabels = TRUE,
                          tickcolor = 'rgb(127,127,127)',
                          ticks = 'outside',
                          zeroline = FALSE),
             yaxis = list(title =ytitle,
                          range = c(0, 20),
                          gridcolor = 'rgb(255,255,255)',
                          showgrid = TRUE,
                          showline = FALSE,
                          showticklabels = TRUE,
                          tickcolor = 'rgb(127,127,127)',
                          ticks = 'outside',
                          zeroline = FALSE,
                          automargin = TRUE))
  return (p)
}

propa.plot.active_distance_compiled <- function (L0, L_bkg, f, d, r0=1, t=20, rh=60, pa=101325, a0=0.02, 
                                                 showlegend=FALSE, fmin=1, fmax=6, dbmin_spl=0, dbmax_spl=120, 
                                                 dbmin_att=0,dbmax_att=100, dmin=0, dmax=1000)
{
  Ageo.dB = propa.Ageo(r=d,r0)$db
  Aatm.dB = rep(NA,length(d))
  Ahab.dB = rep(NA,length(d))
  for (ii in (1:length(d)))
  {
    Aatm.dB[ii] = propa.Aatm(f[ii], r=d[ii],r0, t, rh, pa)$db
    Ahab.dB[ii] = propa.Ahab(f[ii],r=d[ii],r0,a0)$db
  }
  
  fig0 = plot_ly(x=L0, y=f, type = 'scatter', mode = 'lines',orientation = 'h', name="source")
  fig0 = add_trace(fig0, x=L_bkg, y=f, type = 'scatter', mode = 'lines', line = list(dash="dot"),orientation = 'h', name="ambient sound")
  fig0 = layout(fig0, 
                # colorway = mypalette, 
                xaxis = list(title ="Level [dB SPL]",
                             range = c(dbmin_spl, dbmax_spl)),
                yaxis = list(title ="Frequency [kHz]",
                             range = c(fmin, fmax)))
  
  fig1 = plot_ly(x=Ageo.dB, y=f, type = 'bar', orientation = 'h', name="A<sub>geo</sub>")
  fig1 = add_trace(fig1, x=Aatm.dB, y=f, name="A<sub>atm</sub>", type = 'bar', orientation = 'h')
  fig1 = add_trace(fig1, x=Ahab.dB, y=f, name="A<sub>hab</sub>", type = 'bar', orientation = 'h')
  fig1 = layout(fig1, 
                barmode = 'stack',
                xaxis = list(title ="Attenuation [dB SPL]",
                             range = c(dbmin_att,dbmax_att)),
                yaxis = list(title ="Frequency [kHz]",
                             range = c(fmin, fmax)))
  
  fig2 = plot_ly(x=d, y=f, type = 'bar', orientation = 'h', name="Detection distance")
  fig2 = layout(fig2, 
                xaxis = list(title='Distance [m]', 
                             range = c(dmin, dmax)),
                yaxis = list(title ="Frequency [kHz]"))
  
  fig_recap = subplot(fig0, fig1, fig2, shareY = TRUE)
  fig_recap = layout(fig_recap,
                     paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
                     showlegend = showlegend,
                     margin=list(t=0, b=0),
                     xaxis = list(title='Sound Pressure Level [dB]',
                                  gridcolor = 'rgb(255,255,255)',
                                  showgrid = TRUE,
                                  showline = FALSE,
                                  showticklabels = TRUE,
                                  tickcolor = 'rgb(127,127,127)',
                                  ticks = 'outside',
                                  zeroline = FALSE),
                     xaxis2 = list(title='Attenuation [dB]',
                                   gridcolor = 'rgb(255,255,255)',
                                   showgrid = TRUE,
                                   showline = FALSE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = FALSE),
                     xaxis3 = list(title='Distance [m]',
                                   gridcolor = 'rgb(255,255,255)',
                                   showgrid = TRUE,
                                   showline = FALSE,
                                   showticklabels = TRUE,
                                   tickcolor = 'rgb(127,127,127)',
                                   ticks = 'outside',
                                   zeroline = FALSE),
                     yaxis = list(gridcolor = 'rgb(255,255,255)',
                                  showgrid = TRUE,
                                  showline = FALSE,
                                  showticklabels = TRUE,
                                  tickcolor = 'rgb(127,127,127)',
                                  ticks = 'outside',
                                  zeroline = FALSE,
                                  automargin = TRUE))
  return (fig_recap)
}

#=======================================================================================================#
#
#                                           LOAD DATA
#
#=======================================================================================================#

#=================================================================================================
# simulation 
#=================================================================================================
# # ########## howler monkey
# up to 90 dB sound pressure level (SPL) at 5 m of distance (Whitehead 1995 )
# L0 = 104 # howler monkey @5m => 104dB @1m / Production of Loud and Quiet Calls in Howler Monkeys Da cunha 2015
# F0 = 0.3
# F1 = 1
# DELTA_F = 0.1
# R0 = 1
# # data for sonometer in Guyane
# INDEX = 4
# A0 = df$A0[INDEX]
# DATA_FILENAME = FILENAME
# # for display
# dBMAX_SPL=100
# dBMAX_ATT=90
# D_MAX=1000

######### Lipaugus vociferans
# FILENAME_ROOT = "guiana_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 24   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 87     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 1.01340e5 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.019 # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 111.5     # initial souund level (Lipaugus vociferans / MEASURING THE SOUNDPRESSURE LEVEL OFTHE SONG OF THESCREAMING PIHA LIPAUGUSVOCIFERANS: ONE OF THELOUDEST BIRDS IN THE WORLD? Nemeth 2004
# F0 = 1      # min frequency of the call/song/signal 
# F1 = 5     # max frequency of the call/song/signal 
# DELTA_FBIN = 0.1 # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# 
# # for display
# dBMAX_SPL=120
# dBMAX_ATT=90
# D_MAX=750

##########  Allobates Femoralis
# L0 = 92 # Allobates Femoralis/ Ringler & al, Acoustic ranging in poison frogs—it is not about signal amplitude alone, behav Ecol Sociobiol (2017)
# L_thresh = 56 # phonotactic threshold 56-68dB/ Sound radiation pattern of the advertisement call of the highly territorial poison frog Allobates femoralis. Set as background level
# F0 = 2.9
# F1 = 3.9
# DELTA_FBIN = 0.1
# R0=0.5
# dBMAX_SPL=100
# dBMAX_ATT=70
# D_MAX=300

# ########## Song Thrush (Grive musicienne)
# FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 80     # initial souund level (Mathevon & Aubin book 2020 chap 2)
# F0 = 1.95      # min frequency of the call/song/signal (no ref, from sounds found on XenoCanto close to Jura)
# F1 = 4.60     # max frequency of the call/song/signal (no ref, from sounds found on XenoCanto close to Jura)
# DELTA_FBIN = 0.1 # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# 
# # for display
# dBMAX_SPL=80
# dBMAX_ATT=80
# D_MAX=300

# ########## Wood Pigeon
# FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 80     # initial sound level (no ref)
# F0 = 0.3      # min frequency of the call/song/signal (no ref, from sounds found on XenoCanto close to Jura)
# F1 = 0.6      # max frequency of the call/song/signal (no ref, from sounds found on XenoCanto close to Jura)
# DELTA_FBIN = 0.05  # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# 
# # for display
# dBMAX_SPL=90
# dBMAX_ATT=90
# D_MAX=400

##########  pure tone @ 90dB SPL (Maclaren ECOEVOL 2019)
# FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 90     # pure tone @ 90dB SPL avec filtre A (Maclaren ECOEVOL 2019)
# F0 = 7      # min frequency of the call/song/signal
# F1 = 7.1      # max frequency of the call/song/signal
# DELTA_FBIN = 0.1  # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# # for display
# dBMAX_SPL=100
# dBMAX_ATT=100
# D_MAX=500

# ##########  average bird 75B SPL @ 1m with call between 2-8kHz
# filename root of the data
# FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 75     # sound level in dB SPL at R0
# F0 = 2      # min frequency of the call/song/signal
# F1 = 8     # max frequency of the call/song/signal
# DELTA_FBIN = 0.5  # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# # for display
# dBMAX_SPL=100
# dBMAX_ATT=90
# D_MAX=500

##########  average bird 100 SPL @ 1m with call between 2-8kHz
# filename root of the data
# FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
# FILE_DIR = "../data/psd/" 
# FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
# CORRECTION_RECORDER = TRUE # SM4 gain correction
# TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
# RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
# PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
# A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019) 
# L0 = 100     # sound level in dB SPL at R0
# F0 = 2      # min frequency of the call/song/signal
# F1 = 8     # max frequency of the call/song/signal
# DELTA_FBIN = 0.5  # frequency resolution in kHz
# R0= 1       # distance at which the initial sound level L0 was measured (in m)
# # for display
# dBMAX_SPL=100
# dBMAX_ATT=90
# D_MAX=500

##########  white noise @ 80dB SPL between 0-20kHz
# filename root of the data
FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/"
FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
CORRECTION_RECORDER = TRUE # SM4 gain correction
TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019)
L0 = 80     # sound level in dB SPL at R0
F0 = 0      # min frequency of the call/song/signal
F1 = 20     # max frequency of the call/song/signal
DELTA_FBIN = 0.5  # frequency resolution in kHz
R0= 1       # distance at which the initial sound level L0 was measured (in m)
# for display
dBMAX_SPL=80
dBMAX_ATT=80
D_MAX=150

# load data : p.bkg, p.exp, TEMP, RH, PS0, A0
load(FILENAME)

# convert spectrum into specbins
f     = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$f
P_bkg = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$s

# create a vector with the index of the selected distances
DISTANCE_SELECT = DISTANCES>=R0 
# create a vector with the index of the selected frequencies
FREQUENCY_SELECT = f>= F0 & f<=F1

# select the frequencies and the distances
f     = f        [FREQUENCY_SELECT]
P_bkg = P_bkg    [FREQUENCY_SELECT, DISTANCE_SELECT]

# Experimental value of the background noise P_bkg (mean values (average of all distances) for each frequency bins)
L_bkg  =  psd2dBSPL(apply(P_bkg,1,mean), gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
# sound level pressure of the source corresponding to the frequency band
L0_per_bin = propa.dBSPL_per_bin(L=L0, f=f)$db

# Adjust the SPL value to add the frequency response of the SM4 in order to simulate what is recorded/detected by the SM4 (which is bad for high frequencies)
if ((CORRECTION_RECORDER == TRUE) && (grepl("sm4", FILENAME_ROOT) == TRUE))
{
  # load the GENERIC frequency response of the SM4. 
  # The gain to correct the frequency response of the SM4 was measured on a single SM4 recorder, 
  # assuming that all SM4 have the exact same frequency response
  # It would be better to adapt the correction for each SM4 if possible.
  load('./toolbox/SM4_gain.Rdata')
  FREQ_GAIN_CORR = SM4.G[,1]
  GAIN_CORR = SM4.G[,2] 
  # linear interpolation of SM4.G in order to match the frequency bin
  # SM4.gain contains the gain to add to the result obtained with a SM4 
  # in order to correct as much as possible the frequency response of the SM4.
  SM4.gain = -power2dB(approx(FREQ_GAIN_CORR, dB2power(GAIN_CORR), f, rule=2)$y)
  # Set L0 as it is "seen" by the SM4 (higher sensitivity around (1-10kHz), low sensitivity for >10kHz)
  L0_per_bin  = L0_per_bin + SM4.gain
} 

# get the maximum listening distance 
dmax = propa.detection_distance(L_bkg=L_bkg, L0=L0_per_bin, f=f, r0= R0, delta_r=1, t=TEMP, rh=RH, pa=PS0, a0=A0)

# distance min, max, avg
dmax.min = min(dmax[,2])
dmax.max = max(dmax[,2])
dmax.avg = mean(dmax[,2])

fig1 = propa.plot.active_distance_compiled(L_bkg=L_bkg, L0=L0_per_bin, f=dmax[,1], d=dmax[,2], r0= R0, t=TEMP, rh=RH, pa=PS0, a0=A0, 
                                           showlegend=TRUE, fmin=F0-DELTA_FBIN, fmax=F1+DELTA_FBIN, dbmin_spl=0, dbmax_spl=dBMAX_SPL, dbmin_att=0, dbmax_att=dBMAX_ATT, dmin=0, dmax=D_MAX)

layout(fig1,
       bargap = 0,
       autosize = F, 
       width = 790, 
       height = 230)

#=================================================================================================
# Sound pressure level at r0 for the frequency f1
#=================================================================================================
F1 = 5
L = 13
R = 100
R0 = 1
TEMP = 17 # Temperature
RH = 67 # relative humidity
PS0 = 87999 # atmospheric pressure in Pa
A0 = 0.020 # coef attenuation of the habitat (JURA)
Ageo = propa.Ageo(r=R, r0=R0)$db
Aatm = propa.Aatm(f=F1, r=R, r0=R0, t=TEMP, rh=RH, pa= PS0)$db
Ahab = propa.Ahab(f=F1, r=R, r0=R0, a0=A0)$db

# Sound pressure level at r0 for the frequency f1
L0 = L + Ageo + Aatm + Ahab

F1 = 12
L = 8
R = 52
R0 = 1
TEMP = 17 # Temperature
RH = 67 # relative humidity
PS0 = 87999 # atmospheric pressure in Pa
A0 = 0.020 # coef attenuation of the habitat (JURA)
Ageo = propa.Ageo(r=R, r0=R0)$db
Aatm = propa.Aatm(f=F1, r=R, r0=R0, t=TEMP, rh=RH, pa= PS0)$db
Ahab = propa.Ahab(f=F1, r=R, r0=R0, a0=A0)$db

# Sound pressure level at r0 for the frequency f1
L0 = L + Ageo + Aatm + Ahab

#=================================================================================================
# simulation attenuation on real data : ex Lipaugus vociferans

# extract signals

############## Lipaugus vociferans
filename = "XC442394-GUY2014_013_Lipaugus_vociferans.wav" # from Juan Ulloa
FILE_DIR = "../data/guiana/LIPAUGUS/"
R0 = 10   # reference distance
R = 100 # final distance

library(signal, warn.conflicts = F, quietly = T) # signal processing functions
library(oce, warn.conflicts = F, quietly = T) # image plotting functions and nice color maps
library(tuneR) 
library(seewave) 

# extract signals

piha.wav = readWave(paste(FILE_DIR, filename,sep=""), from=50.5, to=52,units = 'seconds')
bkg.wav = readWave(paste(FILE_DIR, filename,sep=""), from=54.5, to=56,units = 'seconds')

fs = piha.wav@samp.rate
p0.piha.wav = piha.wav@left - mean(piha.wav@left)
p0.bkg.wav = bkg.wav@left - mean(bkg.wav@left)

# apply attenuation. Environmental values are the same as in the article
p.piha.wav= propa.apply.att(p0.piha.wav, fs, r=R, r0= R0, t=24, rh=87, pa=101325, a0=0.019)
# add original background to the data after propagation because attenuation process also attenuate background...
p.piha.wav = p.piha.wav + p0.bkg.wav*0.5
p0.piha.wav = p0.piha.wav + p0.bkg.wav*0.5

# Leq
wav2leq(p0.bkg.wav,fs,gain=30, sensitivity = -35)
wav2leq(p0.piha.wav,fs,gain=30, sensitivity = -35)

# create spectrogram with initial song at 50m
p = p0.piha.wav
N = 512
res = specgram(x=p,n=N,Fs=fs,window=N,overlap=N/2)
# frequency vector
f = res$f
# PSD
P = (abs(res$S))
# convert to dB
L = psd2dBSPL(P, gain=30, sensitivity=-35)
# config time axis
t = res$t
# plot spectrogram
imagep(x = t,
       y = res$f,
       z = t(L),
       col = oce.colorsViridis,
       ylim = c(0,11000),
       ylab = 'Frequency [Hz]',
       xlab = 'Time [s]',
       zlim = c(30,65),
       drawPalette = T,
       decimate = T)

# create spectrogram with initial song at 50m
p = p.piha.wav
N = 512
res = specgram(x=p,n=N,Fs=fs,window=N,overlap=N/2)
# frequency vector
f = res$f
# PSD
P = (abs(res$S))
# convert to dB
L = psd2dBSPL(P, gain=26+4, sensitivity=-35)
Leq = psd2leq(P, gain=26+4, sensitivity=-35)
# config time axis
t = res$t
# plot spectrogram
imagep(x = t,
       y = res$f,
       z = t(L),
       col = oce.colorsViridis,
       ylim = c(0,11000),
       ylab = 'Frequency [Hz]',
       xlab = 'Time [s]',
       zlim = c(30,60),
       drawPalette = T,
       decimate = T)

#=================================================================================================
# simulation attenuation on real data : ex Lipaugus vociferans

# extract signals

############## Lipaugus vociferans
#=================================================================================================
# simulation attenuation on real data : ex Lipaugus vociferans
############## Lipaugus vociferans
filename = "S4A09154_20191120_134500.wav"
FILE_DIR = "../data/guiana/LIPAUGUS/"
R0 = 150 # reference distance
R = c(150,200, 400, 600, 800)  # final distances

# extract signals (lipaugus and background)
piha = readWave(paste(FILE_DIR, filename,sep=""), from=4.3, to=6.4,units = 'seconds')
bkg  = readWave(paste(FILE_DIR, filename,sep=""), from=1.8, to=3.9,units = 'seconds')

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
  p.piha.att = p.piha.att + p0.bkg*0.5
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
                 zlim = c(25,47),
                 drawPalette = T,
                 decimate = T,
                 useRaster = T)
  fig[[ii]] = list(image)
}