# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARIES
#
#=======================================================================================================#

library(seewave)
library(tuneR)
library(signal)

#=======================================================================================================#
#
#                                           SET WORKING DIRECTORIES 
#                                           LOAD MY TOOLBOX
#
#=======================================================================================================#

##### Change working directory to the current script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### CALL functions from MY TOOLBOX
source("./toolbox/toolbox_propagation_wav_extraction.R", chdir=T)

#=======================================================================================================#
#
#                                           SET GLOBAL PARAMETERS
#
#=======================================================================================================#
# save results ?
SAVE = TRUE
# Which channel ?
CHANNEL ='right' 
# length of the window for the spectrogram or use for averaging the spectra
NFFT = 2048
# number of repetition per distance
Nrep = 3

# SM4 PARAMETERS
# Sensbility microphone
S = -35   #dBV
# Amplification gain
G = 26+16 # total amplification gain : preamplifier gain = 26dB and gain = 16dB
# Analog to digital convertion (Vpp [peak to peak])
VADC = 2
# Reference pressure (in the air : 20µPa)
P_REF = 20e-6

######## TIME LINE
silbetween_dur <- 0.500 
am_dur <- 0.500
silbetween22s_dur <- 22
sw_dur <- 1.1
silbetweenSweep_dur <- 4.4 
wn_dur <- 22
snkHz_dur <- 0.100
silbetweenSnkHz_dur <- 0.4
snkHz_comp_dur <- 44*(snkHz_dur + silbetweenSnkHz_dur)
snkHz_comp_dur_without_sil <- 44*(snkHz_dur)

one_dur <- am_dur+silbetween22s_dur+4*(sw_dur+silbetweenSweep_dur)+wn_dur+silbetween_dur+ snkHz_comp_dur + silbetween_dur + am_dur
start_ref = am_dur
start_sw =  start_ref+ silbetween22s_dur 
start_wn =  start_sw + 4*(sw_dur+silbetweenSweep_dur) 
start_pt =  start_wn + wn_dur + silbetween_dur

total_dur <- one_dur*3 + 2*silbetween_dur # 4min 30s

# # ########### FILE PARAMETERS
filename = "../data/jura/SM4/site04_S4A09600_NE_long_signals.wav"
savefile = "../data/psd/jura_sm4_wn" 
start_1m = 16.569
start_2m = 5*60 + 0.179
start_5m = 9*60 + 44.250
start_10m = 14*60 + 33.420
start_20m = 19*60 + 19.834
start_30m = 24*60 + 1.911
start_40m = 28*60 + 42.421
start_50m = 33*60 + 25.451
start_60m = 38*60 + 10.651
start_70m = 42*60 + 56.976
start_80m = 47*60 + 43.510
start_90m = 52*60 + 25.136
start_100m = 57*60 + 13.219
DISTANCES <- c(1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100) # in meters
am_timeline <- c(start_1m, start_2m, start_5m, start_10m, start_20m, start_30m, start_40m, start_50m, start_60m, start_70m, start_80m, start_90m, start_100m)

##### Sampling Frequency (Recording)
# read a portion of the wave file in order to know the sampling frequency (fs_rec)
wavSM4 = readWave(filename, from = 0, to = 0.001, units = 'seconds')
fs_rec = wavSM4@samp.rate
bit = wavSM4@bit

# data frame with distances and timeline for the different tests
df_timeline = data.frame(DISTANCES, am_timeline, am_timeline+one_dur+silbetween_dur, am_timeline+2*(one_dur+silbetween_dur))
colnames(df_timeline) <- c("distance", "am_timeline1", "am_timeline2", "am_timeline3")

#====================================================================================================#
#
#                             EXTRACT SOUND 
#
#====================================================================================================#

for (ii in 1:(length(df_timeline)-1))
{
  print(paste('timeline:',ii))
  
  timeline = df_timeline[[ii+1]]
  distance = df_timeline$distance
  
  #================= Extract background noise signal (wave)
  print('Extracting Background noise...')
  bkg.wav  = extractNoise(filename, channel =CHANNEL, am_timeline=timeline, total_dur=total_dur, t0=start_ref*2, dur=silbetween22s_dur-start_ref)
  # convert sounds (wave) into SPL and subtract DC offset
  bkg.wav = bkg.wav - mean(bkg.wav)
  # Leq
  #bkg.wav.Leq = apply(bkg.wav, 2, wav2leq, fs_rec, G, deltaT, S, bit, VADC, P_REF)
  # average Leq from wav
  #bkg.wav.Leq = power2dB(colMeans(dB2power(bkg.wav.Leq)))
  # Power Density Spetrum : PSD
  res = meanPSDfromNoise(bkg.wav, fs_rec, NFFT, timeline)
  bkg.freq = res$freq
  bkg.PSD  = res$PSD
  # average Leq from PSD
  #bkg.PSD.Leq  = apply(bkg.PSD, 2, psd2leq, G, S, bit, VADC, P_REF)
  
  #================= Extract white noise signal (wave)
  print('Extracting White noise...')
  sig.wav= extractNoise(filename, channel =CHANNEL, am_timeline=timeline, total_dur=total_dur, t0=start_wn, dur=wn_dur)
  # convert sounds (wave) into SPL and subtract DC offset
  sig.wav  = sig.wav  - mean(sig.wav)
  # Leq
  #sig.wav.Leq.mat = apply(sig.wav, 2, wav2leq, fs_rec, G, deltaT, S, bit, VADC, P_REF)
  # average Leq from wav
  #sig.wav.Leq = power2dB(colMeans(dB2power(sig.wav.Leq.mat)))
  # Power Density Spetrum : PSD
  res = meanPSDfromNoise(sig.wav, fs_rec, NFFT, timeline)
  sig.freq = res$freq
  sig.PSD  = res$PSD
  # average Leq from PSD
  #sig.PSD.Leq  = apply(sig.PSD, 2, psd2leq, G, S, bit, VADC, P_REF)
  
  ########## Save data
  if (SAVE==TRUE)
  {
    dir.create(file.path('../data/', 'psd'))
    FREQUENCY = bkg.freq
    save(filename, 
         DISTANCES, FREQUENCY,
         CHANNEL, Nrep,
         NFFT, fs_rec,
         S, G, bit, VADC, P_REF, 
         bkg.PSD, sig.PSD,
         #bkg.PSD.Leq, sig.PSD.Leq,
         #bkg.wav.Leq, sig.wav.Leq,
         file=paste(savefile, '_timeline',ii,'.Rdata', sep=""))
  }
}