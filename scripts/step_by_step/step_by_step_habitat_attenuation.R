# clear workspace
rm(list = ls())

#=======================================================================================================#
#                                           LOAD LIBRARIES
#=======================================================================================================#
library(seewave)
library(tuneR)
library(signal)
if (!require("robustbase")) install.packages("robustbase")
library(robustbase)

#=======================================================================================================#
#                                           SET WORKING DIRECTORIES 
#                                           LOAD MY TOOLBOX
#=======================================================================================================#
##### Change working directory to the current script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##### CALL functions from MY TOOLBOX
source("../toolbox/toolbox_propagation_wav_extraction.R", chdir=T)
source("../toolbox/toolbox_propa.R", chdir=T)

#=======================================================================================================#
#                                           SET GLOBAL PARAMETERS
#=======================================================================================================#
# save results ?
SAVE = TRUE
# plot results ?
PLOT = TRUE
# length of the window for the spectrogram or use for averaging the spectra
NFFT = 2048
# BIN size when transforming PSD into histogram
DELTA_FBIN = 0.5 # in kHz
# Select the frequency range
F_RANGE = c(0.5,15)
# reference distance used to estimate a0
LIST_R0 = c(10,20,30,40,50)
# correction of the background noise
CORRECTION_BKG = TRUE

# ENVIRONMENTAL FACTORS
# Temperature
TEMP = 17.0 # 
# relative humidity
RH = 67# 
# atmospheric pressure in Pa
PS0 = 87999

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
start_marker_duration <- 0.5
ambient_sound_duration <- 22
wn_duration <- 22
total_dur <- (start_marker_duration + ambient_sound_duration + wn_duration)*10
# # ########### FILE PARAMETERS
filename = "../../data/jura/FIELD_RECORDINGS/field_recordings_SM4_jura.wav"
savefile = "../../data/psd/jura_sm4_wn" 
start_10m = 0.721
start_20m = 45.276
start_30m = 1*60 + 29.952
start_40m = 2*60 + 15.105
start_50m = 2*60 + 59.620
start_60m = 3*60 + 44.075
start_70m = 4*60 + 28.913
start_80m = 5*60 + 13.216
start_90m = 5*60 + 57.999
start_100m = 6*60 + 42.922
DISTANCES <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100) # in meters
am_timeline <- c(start_10m, start_20m, start_30m, start_40m, start_50m, start_60m, start_70m, start_80m, start_90m, start_100m)

##### Sampling Frequency (Recording)
# read a portion of the wave file in order to know the sampling frequency (fs_rec)
wavSM4 = readWave(filename, from = 0, to = 0.001, units = 'seconds')
fs_rec = wavSM4@samp.rate
bit = wavSM4@bit

# data frame with distances and timeline for the different tests
df_timeline = data.frame(DISTANCES, am_timeline)
colnames(df_timeline) <- c("distance", "am_timeline")

#====================================================================================================#
#                             EXTRACT SOUND 
#====================================================================================================#

for (ii in 1:(length(df_timeline)-1))
{
  print(paste('timeline:',ii))
  
  timeline = df_timeline[[ii+1]]
  distance = df_timeline$distance
  
  #================= Extract background noise signal (wave)
  print('Extracting Background noise...')
  bkg.wav  = extractNoise(filename, 
                          am_timeline=timeline, 
                          total_dur=total_dur, 
                          t0=start_10m+start_marker_duration+1, 
                          dur=ambient_sound_duration-1)
  # convert sounds (wave) into SPL and subtract DC offset
  bkg.wav = bkg.wav - mean(bkg.wav)
  # Power Density Spetrum : PSD
  res = meanPSDfromNoise(bkg.wav, fs_rec, NFFT, timeline)
  bkg.freq = res$freq
  bkg.PSD  = res$PSD
  
  #================= Extract white noise signal (wave)
  print('Extracting White noise...')
  sig.wav= extractNoise(filename, 
                        am_timeline=timeline, 
                        total_dur=total_dur, 
                        t0=start_10m+start_marker_duration+ambient_sound_duration+1, 
                        dur=wn_duration-1)
  # convert sounds (wave) into SPL and subtract DC offset
  sig.wav  = sig.wav  - mean(sig.wav)
  # Power Density Spetrum : PSD
  res = meanPSDfromNoise(sig.wav, fs_rec, NFFT, timeline)
  sig.freq = res$freq
  sig.PSD  = res$PSD
  
  ########## Save data
  if (SAVE==TRUE)
  {
    FREQUENCY = bkg.freq
    save(filename, 
         DISTANCES, FREQUENCY,
         NFFT, fs_rec,
         bkg.PSD, sig.PSD,
         #bkg.PSD.Leq, sig.PSD.Leq,
         #bkg.wav.Leq, sig.wav.Leq,
         file=paste(savefile, '_step_by_step','.Rdata', sep=""))
  }
}

#====================================================================================================#
#               convert spectrum into specbins
#====================================================================================================#

f         = specbin(bkg.PSD, FREQUENCY, DELTA_FBIN)$f
P_bkg     = specbin(bkg.PSD, FREQUENCY, DELTA_FBIN)$s
P         = specbin(sig.PSD, FREQUENCY, DELTA_FBIN)$s

# keep only the frequencies corresponding to the frequency range (F_RANGE)
index_select = f>=F_RANGE[1] & f<=F_RANGE[2]
f         = f[index_select]
P_bkg     = P_bkg[index_select,]
P         = P[index_select,]

#====================================================================================================#
#               Attenuation and Excess Attenuation (EA)
#====================================================================================================#
first = 1
for (dd in 1:length(LIST_R0))
{ 
  DISTANCE_SELECT = DISTANCES>=LIST_R0[dd] 
  
  # select distances use for the calculation : VECTOR of size M                                   
  r = DISTANCES[DISTANCE_SELECT]         
  # select the initial distance
  r0 = r[1]
  
  #### ENERGY (no unit)
  # P.select:  energy of the propagated signal : MATRIX [N*M]
  P.select  = P[,DISTANCE_SELECT]           
  # P_bkg.select : energy of the ambient sound 
  P_bkg.select = P_bkg[,DISTANCE_SELECT]  
  
  if (CORRECTION_BKG == TRUE)
  {
    #### CORRECTED SIGNAL AFTER AMBIENT SOUND SUBRACTION
    # Set to NA the data below the ambient sound in order to avoid using them for the calculation
    index = P.select < 4*P_bkg.select         
    P.select[index] = NA      
    # subtract the noise level to the original sound level
    P.select = P.select-P_bkg.select 
  }
  
  # geometrical or spheric or spreading loss attenuation : MATRIX [N*M]
  Ageo.dB = propa.Ageo(r=r, r0=r0)$db
  Ageo.dB = matrix(rep(Ageo.dB,each=length(f)),nrow=length(f))
  # atmospheric attenuation 
  Aatm.dB = propa.Aatm(f=f, r=r, r0=r0, t=TEMP, rh=RH, pa=PS0)$db
  
  # Convert psd into dB scale
  L = power2dB(P.select)
  # create a matrix with L0 (for distance r0)
  L0 = t(matrix(rep(L[,1],each=length(r)),nrow=length(r)))
  
  # excess attenuation : MATRIX [N*M]
  EA =  L0 - L - Ageo.dB  - Aatm.dB
  
  # create matrix of r and f in order to be able to divide the matrix EA
  r_mat = matrix(rep(r,each=length(f)),nrow=length(f))
  f_mat = t(matrix(rep(f,each=length(r)),nrow=length(r)))*1000
  
  if (first == 1)
  {
    EA_per_d = power2dB(rowMeans(dB2power(EA/(r_mat-r0)), na.rm=TRUE))
    EA_per_d.mean = matrix(EA_per_d)
    first = 0
  }
  else
  {
    EA_per_d = power2dB(rowMeans(dB2power(EA/(r_mat-r0)), na.rm=TRUE))
    EA_per_d.mean = cbind(EA_per_d.mean, EA_per_d)
  }
}

########################################################################""

# robust linear regression with lmrob (library(robustbase))
y = EA_per_d.mean
x = f
y[y==0] = NA
y[which(is.nan(y))] = NA
y[which(y==Inf)] = NA
y = as.vector(y)
x = as.vector(x)
data = data.frame(f=x, EA=y)
model = lmrob(EA~f+0, data, setting = "KS2014", na.rm=TRUE)
summary(model)

# a0
if (summary(model)$coefficients[4]>0.001)
{
  sprintf('test de linéarité rejeté')
  a0 = NA
  r.squared= NA
  pvalue = "  p>0.001"
}else
{
  a0 = coef(model)['f']
  r.squared = summary(model)$r.squared
  pvalue = "  p<0.001"
  cat(sprintf('a0 = %2.3f \nr² = %1.2f \np = %1.4f',a0, r.squared, summary(model)$coefficients[4]))
}



#====================================================================================================#
#               PLOT
#====================================================================================================#
if (PLOT == TRUE)
{
  if (!require("plotly")) install.packages("plotly")
  library(plotly)
  
  XAXIS_FREQ = list(title = "Frequency [kHz]",
                    range = c(0, 20),
                    gridcolor = 'rgb(255,255,255)',
                    showgrid = TRUE,
                    showline = FALSE,
                    showticklabels = TRUE,
                    tickcolor = 'rgb(127,127,127)',
                    ticks = 'outside',
                    zeroline = FALSE,
                    automargin = TRUE)
  
  YAXIS_SPECTRUM = list(title = "SPL [dB]",
                        range = c(-5, 70),
                        gridcolor = 'rgb(255,255,255)',
                        showgrid = TRUE,
                        showline = FALSE,
                        showticklabels = TRUE,
                        tickcolor = 'rgb(127,127,127)',
                        ticks = 'outside',
                        zeroline = FALSE)
  
  XAXIS_DIST = list(title = "Distances [m]",
                    gridcolor = 'rgb(255,255,255)',
                    showgrid = TRUE,
                    showline = FALSE,
                    showticklabels = TRUE,
                    tickcolor = 'rgb(127,127,127)',
                    ticks = 'outside',
                    zeroline = FALSE)
  
  #====================================================================================================#
  #               PLOT LEQ
  #====================================================================================================#
  sig.PSD.Leq = rep(NA,length(DISTANCES))
  bkg.PSD.Leq = rep(NA,length(DISTANCES))
  
  for (ii in 1:length(DISTANCES))
  {
    sig.PSD.Leq[ii] = psd2leq(sig.PSD[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
    bkg.PSD.Leq[ii] = psd2leq(bkg.PSD[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  }
  
  p <-
    plot_ly(
      x = DISTANCES,
      y = sig.PSD.Leq ,
      type = 'scatter',
      mode = 'lines',
      name = 'White noise',
      showlegend = TRUE
    )
  p <-
    add_trace(
      p,
      y = bkg.PSD.Leq ,
      type = 'scatter',
      mode = 'lines',
      name = 'Ambient sound',
      showlegend = TRUE
    )
  
  p3 <- layout(
    p,
    autosize = F,
    width = 450,
    height = 250,
    paper_bgcolor = 'rgb(255,255,255)',
    plot_bgcolor = 'rgb(225,225,225)',
    legend = list(x = 0.6, y = 0.05),
    xaxis = XAXIS_DIST,
    yaxis = YAXIS_SPECTRUM
  )
  
  p3
  
  #====================================================================================================#
  #               Plot raw spectrum for each distance
  #====================================================================================================#
  
  # vector r with all distances
  r = DISTANCES
  
  # convert spectrum into specbins
  f     = specbin(bkg.PSD, FREQUENCY, DELTA_FBIN)$f
  P_bkg = specbin(bkg.PSD, FREQUENCY, DELTA_FBIN)$s
  P     = specbin(sig.PSD, FREQUENCY, DELTA_FBIN)$s
  
  # Convert psd (P) into dB sound pressure level (L)
  # get the original sound level (with noise) Lexp 
  Lexp = psd2dBSPL(P, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 
  # get the noise level
  L_bkg = psd2dBSPL(P_bkg, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 
  
  #### CORRECTED SIGNAL AFTER AMBIENT SOUND SUBRACTION
  # Set to NA the data below the ambient sound in order to avoid using them for the calculation
  index = P< P_bkg         
  P[index] = NA      
  # subtract the noise level to the original sound level
  P_corr= P-P_bkg 
  L_corr = psd2dBSPL(P_corr, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF)   
  
  fig_spectrum <- plot_ly(width = 450, height = 350,
                          x = f,  y = L[,1], type = 'scatter', mode = 'lines+markers', name = paste(r[1],'m'),showlegend = TRUE, marker=list(opacity=0.75, symbol =1, color= 1))
  
  for (ii in 2:length(r))
  {
    fig_spectrum <- add_trace(fig_spectrum, x =f , y =L[,ii], type = 'scatter', mode = 'lines+markers', name = paste(r[ii],'m'),showlegend = TRUE, marker=list(opacity=0.75, symbol =I(ii), color= ii))
  }
  
  fig_spectrum =layout(fig_spectrum, 
                       paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
                       legend = list(x = 1, y = 0.9,
                                     automargin = FALSE,
                                     font = list(size = 10)),
                       xaxis = XAXIS_FREQ ,
                       yaxis =  YAXIS_SPECTRUM)
  show(fig_spectrum)
  
  #====================================================================================================#
  #               Plot mean SPL vs frequencies (frequency bins)
  #====================================================================================================#
  
  # map of dB SPL vs distance and frequency
  fig_L = plot_ly(x = r, y = f, z =L, type = "heatmap", zmin=0, zmax=60) %>% hide_colorbar()
  
  fig_L_bkg = plot_ly(x = r, y = f, z =L_bkg, type = "heatmap", zmin=0, zmax=60) %>% hide_colorbar()
  
  fig_L_corr = plot_ly(x = r, y = f, z =L_corr,type = "heatmap", zmin=0, zmax=60, colorbar=list(title='dB SPL')) 
  
  fig_LdBSPL =subplot(fig_L, fig_L_bkg, fig_L_corr, shareY = TRUE)
  
  fig_LdBSPL =layout(fig_LdBSPL, 
                     width = 450, height = 350, 
                     paper_bgcolor='rgb(255,255,255)', 
                     plot_bgcolor='rgb(229,229,229)',
                     xaxis2 = list(title = "Distances [m]"),
                     yaxis = XAXIS_FREQ)  
  
  show(fig_LdBSPL)
  
  
  
}
