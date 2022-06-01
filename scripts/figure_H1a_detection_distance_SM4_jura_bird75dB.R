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
FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/" 
FILENAME = paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep="") # filename
CORRECTION_RECORDER = TRUE # SM4 gain correction
TEMP = 17   # Temperature (JURA : 17 / GUIANA : 24 )
RH = 67     # relative humidity (JURA : 67 / GUIANA : 87 )
PS0 = 87999 # atmospheric pressure in Pa (JURA :87999 / GUIANA : 1.01340e5)
A0 = 0.024  # coef attenuation of the habitat (SM4 JURA : 0.024 / SVANTEK JURA : 0.020 / SM4 GUIANA : 0.011 / SVANTEK GUIANA : 0.019)
L0 = 75     # sound level in dB SPL at R0
F0 = 2      # min frequency of the call/song/signal
F1 = 8     # max frequency of the call/song/signal
DELTA_FBIN = 0.1  # frequency resolution in kHz
R0= 1       # distance at which the initial sound level L0 was measured (in m)
# for display
dBMAX_SPL=100
dBMAX_ATT=90
D_MAX=400

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

p = propa.plot.active_distance_compiled(L_bkg=L_bkg, L0=L0_per_bin, f=dmax[,1], d=dmax[,2], r0= R0, t=TEMP, rh=RH, pa=PS0, a0=A0, 
                                    showlegend=TRUE, fmin=F0-DELTA_FBIN, fmax=F1+DELTA_FBIN, 
                                    dbmin_spl=0, dbmax_spl=dBMAX_SPL, dbmin_att=0, dbmax_att=dBMAX_ATT, dmin=0, dmax=D_MAX)

layout(p,
       bargap = 0,
       autosize = F, 
       width = 790, 
       height = 230)