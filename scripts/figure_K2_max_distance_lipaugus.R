# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARY
#
#=======================================================================================================#
if (!require("plotly")) install.packages("plotly")
if (!require("tuneR")) install.packages("signal")
library(tuneR, warn.conflicts = F, quietly = T) 
library(plotly)
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

#=======================================================================================================#
#
#                                           SET GLOBAL FUNCTIONS
#
#=======================================================================================================#

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
                # colorway = mypalette, 
                barmode = 'stack',
                xaxis = list(title ="Attenuation [dB SPL]",
                             range = c(dbmin_att,dbmax_att)),
                yaxis = list(title ="Frequency [kHz]",
                             range = c(fmin, fmax)))
  
  fig2 = plot_ly(x=d, y=f, type = 'bar', orientation = 'h', name="Detection distance")
  fig2 = layout(fig2, 
                # colorway = mypalette, 
                xaxis = list(title='Distance [m]', 
                             range = c(dmin, dmax)),
                yaxis = list(title ="Frequency [kHz]"))
  
  fig_recap = subplot(fig0, fig1, fig2, shareY = TRUE)
  fig_recap = layout(fig_recap,
                     # colorway = mypalette, 
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
#                                           SIMULATION
#
#=======================================================================================================#
#########  white noise @ 80dB SPL between 0-20kHz
# filename root of the data
FILENAME_ROOT = "guiana_svantek_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/" 
FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
CORRECTION_RECORDER = TRUE # SM4 gain correction
TEMP = 24   # Temperature (JURA : 17 / GUIANA : 24 )
RH = 87     # relative humidity (JURA : 67 / GUIANA : 87 )
PS0 = 1.01340e5 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
A0 = 0.019  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019)
L0 = 111.5     # sound level in dB SPL at R0
F0 = 1.5      # min frequency of the call/song/signal
F1 = 5     # max frequency of the call/song/signal
DELTA_FBIN = 0.1  # frequency resolution in kHz
R0 = 1       # distance at which the initial sound level L0 was measured (in m)
# for display
dBMAX_SPL=110
dBMAX_ATT=100
D_MAX=800

# load data : p.bkg, p.exp, TEMP, RH, PS0, A0
filename = "../data/guiana/LIPAUGUS/S4A09154_20191120_134500.wav"
bkg  = readWave(filename, from=1.8, to=3.9,units = 'seconds')

# get the sampling frequency and the signals
fs      = bkg@samp.rate
bkg     = bkg@left  - mean(bkg@left)

# compute the mean spectrogram of the background signal
N = 512
f = meanspec(bkg, f=fs, wl=N, wn="hamming", correction="amplitude", FUN = mean, norm=FALSE, plot=FALSE)[,1] # frequency vector
S = meanspec(bkg, f=fs, wl=N, wn="hamming", correction="amplitude", FUN = mean, norm=FALSE, plot=FALSE)[,2] 
P_bkg = S**2/2 

# convert spectrum into specbins
fbin  = specbin(P_bkg, f, DELTA_FBIN)$f
P_bkg = specbin(P_bkg, f, DELTA_FBIN)$s

# create a vector with the index of the selected frequencies
FREQUENCY_SELECT = fbin>= F0 & fbin<=F1

# select the frequencies and the distances
fbin  = fbin  [FREQUENCY_SELECT]
P_bkg = P_bkg [FREQUENCY_SELECT]

# Experimental value of the background noise P_bkg (mean values (average of all distances) for each frequency bins)
L_bkg  =  psd2dBSPL(P_bkg, gain=30, sensitivity=-35)
# sound level pressure of the source corresponding to the frequency band
L0_per_bin = propa.dBSPL_per_bin(L=L0, f=fbin)$db

# get the maximum listening distance 
dmax = propa.detection_distance(L_bkg=L_bkg, L0=L0_per_bin, f=fbin, r0= R0, delta_r=1, t=TEMP, rh=RH, pa=PS0, a0=A0)

# distance min, max, avg
dmax.min = min(dmax[,2])
dmax.max = max(dmax[,2])
dmax.avg = mean(dmax[,2])

p = propa.plot.active_distance_compiled(L_bkg=L_bkg, L0=L0_per_bin, f=dmax[,1], d=dmax[,2], r0= R0, t=TEMP, rh=RH, pa=PS0, a0=A0, 
                                    showlegend=TRUE, fmin=F0-DELTA_FBIN, fmax=F1+DELTA_FBIN, 
                                    dbmin_spl=0, dbmax_spl=dBMAX_SPL, dbmin_att=0, dbmax_att=dBMAX_ATT, dmin=0, dmax=D_MAX)

layout(p,
       bargap = 0,
       autosize = F, 
       width = 790, 
       height = 230)