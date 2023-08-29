rm(list = ls()) # clear workspace

if (!require("robustbase")) install.packages("robustbase")
library(robustbase)
if (!require("plotly")) install.packages("plotly")
library(plotly)  



source("../toolbox/toolbox_propagation_wav_extraction.R", chdir=T)
source("../toolbox/toolbox_propa.R", chdir=T)


SAVE = TRUE # save results ?
NFFT = 2048 # length of the window for the spectrogram or use for averaging the spectra
DELTA_FBIN = 0.5 # BIN size in kHz when transforming PSD into histogram
F_RANGE = c(0.5,15) # Select the frequency range in kHz
LIST_R0 = c(10,20,30,40,50) # reference distance used to estimate the habitat attenuation coefficient a0
CORRECTION_BKG = TRUE # correction of the background noise

####### ENVIRONMENTAL FACTORS
TEMP = 17.0 # Temperature
RH = 67 # relative humidity
PS0 = 87999 # atmospheric pressure in Pa

####### AUTONOMOUS RECORDER UNIT PARAMETERS 
# IN CASE OF THE SM4 :
S = -35  ## Sensibility microphonein dBV
G = 26+16 # total amplification gain : preamplifier gain = 26dB and gain = 16dB
VADC = 2 # Analog to digital conversion (Vpp [peak to peak])
P_REF = 20e-6  # Reference pressure (in the air : 20µPa)

######## Simulation of the detection distance
L0.simu = 80 # sound level in dB SPL at R0
R0.simu = 1  # distance at which the initial sound level L0 was measured (in m)
F0.simu = 0       # min frequency of the call/song/signal
F1.simu = 20      # max frequency of the call/song/signal
DELTA_FBIN.simu = 0.5  # frequency resolution in kHz

# for display
dBMAX_SPL=80
dBMAX_ATT=80
D_MAX=150
CORRECTION_RECORDER = TRUE

######## TIME LINE
marker_duration = 0.5 # marker duration at the beginning of each white noise broadcasting (in s)
ambient_sound_duration = 22 # ambient sound duration before broadcasting white noise (in s)
wn_duration = 22 # white noise duration (in s)
DISTANCES <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100) # distances between the ARU and the loudspeaker (in meters)
total_dur = (marker_duration + ambient_sound_duration + wn_duration)*length(DISTANCES)

###### FILE PARAMETERS
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
timeline <- c(start_10m, start_20m, start_30m, start_40m, start_50m, start_60m, start_70m, start_80m, start_90m, start_100m)

##### PARAMETERS OF THE AUDIO RECORDINGS
# read a portion of the wave file in order to know the sampling frequency (fs_rec)
wavSM4 = readWave(filename, from = 0, to = 0.001, units = 'seconds')
fs_rec = wavSM4@samp.rate
bit = wavSM4@bit

#================= Extract background noise signal (wave)
print('Extracting Background noise...')
bkg.wav  = extractNoise(filename, 
                        am_timeline=timeline, 
                        total_dur=total_dur, 
                        t0=start_10m+marker_duration+1, 
                        dur=ambient_sound_duration-1)
# convert sounds (wave) into SPL and subtract DC offset
bkg.wav = bkg.wav - mean(bkg.wav)
# Power Density Spectrum : PSD
res = meanPSDfromNoise(bkg.wav, fs_rec, NFFT, timeline)
bkg.freq = res$freq
bkg.PSD  = res$PSD

#================= Extract white noise signal (wave)
print('Extracting White noise...')
sig.wav= extractNoise(filename, 
                      am_timeline=timeline, 
                      total_dur=total_dur, 
                      t0=start_10m+marker_duration+ambient_sound_duration+1, 
                      dur=wn_duration-1)
# convert sounds (wave) into SPL and subtract DC offset
sig.wav  = sig.wav  - mean(sig.wav)
# Power Density Spectrum : PSD
res = meanPSDfromNoise(sig.wav, fs_rec, NFFT, timeline)
sig.freq = res$freq
sig.PSD  = res$PSD



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



sig.PSD.Leq = rep(NA,length(DISTANCES))
bkg.PSD.Leq = rep(NA,length(DISTANCES))

for (ii in 1:length(DISTANCES))
{
  sig.PSD.Leq[ii] = psd2leq(sig.PSD[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  bkg.PSD.Leq[ii] = psd2leq(bkg.PSD[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
}

fig_LEQ <-
  plot_ly(
    x = DISTANCES,
    y = sig.PSD.Leq ,
    type = 'scatter',
    mode = 'lines',
    name = 'White noise',
    showlegend = TRUE,
    width = 450,
    height = 250,
  )
fig_LEQ <-
  add_trace(
    fig_LEQ,
    y = bkg.PSD.Leq ,
    type = 'scatter',
    mode = 'lines',
    name = 'Ambient sound',
    showlegend = TRUE
  )

fig_LEQ <- layout(
  fig_LEQ,
  autosize = F,
  paper_bgcolor = 'rgb(255,255,255)',
  plot_bgcolor = 'rgb(225,225,225)',
  legend = list(x = 0.6, y = 0.05),
  xaxis = XAXIS_DIST,
  yaxis = YAXIS_SPECTRUM
)

fig_LEQ


# vector r with all distances
r = DISTANCES

# convert spectrum into specbins
f     = specbin(bkg.PSD, bkg.freq, DELTA_FBIN)$f
P_bkg = specbin(bkg.PSD, bkg.freq, DELTA_FBIN)$s
P     = specbin(sig.PSD, sig.freq, DELTA_FBIN)$s

# Convert psd (P) into dB sound pressure level (L)
# get the original sound level (with noise) Lexp 
Lexp = psd2dBSPL(P, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 

# Plot
fig_spectrum <- plot_ly(width = 450, height = 350,
                        x = f,  y = Lexp[,1], type = 'scatter', mode = 'lines+markers', name = paste(r[1],'m'),showlegend = TRUE, marker=list(opacity=0.75, symbol =1, color= 1))

for (ii in 2:length(r))
{
  fig_spectrum <- add_trace(fig_spectrum, x =f , y =Lexp[,ii], type = 'scatter', mode = 'lines+markers', name = paste(r[ii],'m'),showlegend = TRUE, marker=list(opacity=0.75, symbol =I(ii), color= ii))
}

fig_spectrum =layout(fig_spectrum, 
                     paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
                     legend = list(x = 1, y = 0.9,
                                   automargin = FALSE,
                                   font = list(size = 10)),
                     xaxis = XAXIS_FREQ ,
                     yaxis =  YAXIS_SPECTRUM)
fig_spectrum

# vector r with all distances
r = DISTANCES

# convert spectrum into specbins
f     = specbin(bkg.PSD, bkg.freq, DELTA_FBIN)$f
P_bkg = specbin(bkg.PSD, bkg.freq, DELTA_FBIN)$s
P     = specbin(sig.PSD, sig.freq, DELTA_FBIN)$s

# Convert psd (P) into dB sound pressure level (L)
# get the original sound level (with noise) Lexp and noise level Ln
Lexp = psd2dBSPL(P, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 
Ln = psd2dBSPL(P_bkg, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF)

# subtract the noise level to the original sound level
index = P < P_bkg        
P[index] = NA  
P_corr = P - P_bkg         
Lexp_corr = psd2dBSPL(P_corr, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF)  

# map of dB SPL vs distance and frequency
fig_Lexp      = plot_ly(x = r, y = f, z =Lexp,     type = "heatmap", zmin=0, zmax=60, width = 450, height = 350) %>% hide_colorbar()
fig_Ln        = plot_ly(x = r, y = f, z =Ln,       type = "heatmap", zmin=0, zmax=60) %>% hide_colorbar()
fig_Lexp_corr = plot_ly(x = r, y = f, z =Lexp_corr,type = "heatmap", zmin=0, zmax=60, colorbar=list(title='dB SPL')) 
fig_LdBSPL =subplot(fig_Lexp, fig_Ln, fig_Lexp_corr, shareY = TRUE)
fig_LdBSPL =layout(fig_LdBSPL, 
                   paper_bgcolor='rgb(255,255,255)', 
                   plot_bgcolor='rgb(229,229,229)',
                   xaxis2 = list(title = "Distances [m]"),
                   yaxis = XAXIS_FREQ)  

fig_LdBSPL



