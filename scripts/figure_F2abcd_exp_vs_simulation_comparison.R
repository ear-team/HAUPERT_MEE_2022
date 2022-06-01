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
#                                           SET GLOBAL PARAMETERS
#
#=======================================================================================================#

# filename root of the data
FILENAME_ROOT = "jura_svantek_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/" 
# SM4 gain correction
CORRECTION_RECORDER = TRUE
# List of frequency bands used for the simulation
F_SIMU = cbind( 
               # c(0,20),
               # c(0,10),
               # c(10,20),
               # c(1,12)
               c(0,5),
               c(5,10),
               c(10,15),
               c(15,20)
               )
# BIN size when transforming PSD into histogram
DELTA_FBIN = 0.1 # in kHz
# Select the minimum distance that will be used for all the calculation : 1m
DISTANCE_MIN = 1

# Set the right environmental values depending on the habitat
if (grepl("jura", FILENAME_ROOT) == TRUE)
{
  # Temperature
  TEMP = 17 # 
  # relative humidity
  RH = 67# 
  # atmospheric pressure in Pa
  PS0 = 87999
  # Initial sound pressure level
  L0 = 78
  # frequency bandwith (kHz)
  FBW = 20
  if (grepl("sm4", FILENAME_ROOT) == TRUE) 
  {
    A0 = 0.024 #0.0027
  } else # svantek
  {
    A0 = 0.020 #0.0023
  }
} else # french guiana
{
  # Temperature
  TEMP = 24 # 
  # relative humidity
  RH = 87# 
  # atmospheric pressure in Pa
  PS0 = 1.01340e5
  # Initial sound pressure level
  L0 = 83  
  # frequency bandwith (kHz)
  FBW = 20
  if (grepl("sm4", FILENAME_ROOT) == TRUE) 
  {
    A0 = 0.011 # 0.0013
  } else # svantek
  {
    A0 = 0.019 # 0.0021
  }
}

#====================================================================================================#
#                                       LOAD Data
#               variables : 
#                   - filename, 
#                   - DISTANCES, FREQUENCY,
#                   - CHANNEL, Nrep,
#                   - NFFT, fs_rec,
#                   - S, G, bit, VADC, P_REF, 
#                   - bkg.PSD.mean, bkg.PSD.std, sig.PSD.mean, sig.PSD.std,
#
#====================================================================================================#
load(paste(FILE_DIR, FILENAME_ROOT, '_average', '.Rdata', sep=""))

# create a vector with the index of the selected distances
DISTANCE_SELECT = DISTANCES>=DISTANCE_MIN 
##### EXPERIMENTAL DATA
# select distances use for the calculation : VECTOR of size M                                   
r.exp = DISTANCES[DISTANCE_SELECT]           
# select the initial distance
r0.exp = r.exp[1]
##### SIMULATION DATA
# select distances use for the calculation : VECTOR of size M   
r.simu = seq(r.exp[1], r.exp[length(r.exp)],5)
# select the initial distance
r0.simu = r.simu[1]

fig = list()
for (ff in 1:dim(F_SIMU)[2])
{
  # convert spectrum into specbins
  f         = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$f
  P_bkg     = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$s
  P         = specbin(sig.PSD.mean, FREQUENCY, DELTA_FBIN)$s
  
  # keep only the frequencies corresponding to the frequency range (F_SIMU)
  index_select = f>= F_SIMU[1,ff] & f<=F_SIMU[2,ff]
  f         = f        [index_select]
  P_bkg     = P_bkg    [index_select, DISTANCE_SELECT]
  P         = P        [index_select, DISTANCE_SELECT]  
  
  #************************* VALUES FROM EXPERIMENTAL DATA ********************************************
  # Experimental value of P with background noise
  L.exp      = psd2dBSPL(apply(P,2,sum),     gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  # Experimental value of the background noise P_bkg
  L_bkg.exp =  psd2dBSPL(apply(P_bkg,2,sum), gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  # Get the mean value along all distance of propagation
  L_bkg.exp.mean = power2dB(mean(dB2power(L_bkg.exp)))
  
  # correct the experimental SPL value to remove as much as possible the SM4 frequency response (which is bad for high frequencies)
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
    
    # compute the mean gain over frequencies
    SM4.gain.mean = power2dB(mean(dB2power(SM4.gain)))
    
    # correct the gain
    L.exp = L.exp - SM4.gain.mean
    L_bkg.exp.mean = L_bkg.exp.mean - SM4.gain.mean
  }  
  
  #************************* VALUES FROM SIMULATION  **************************************
  # SPREADING LOSS 
  Ageo.dB = propa.Ageo(r.simu, r0.simu)$db
  # Mean ATMOSPHERIC ATTENUATION overthe frequency bandwidth   
  Aatm.dB = power2dB(apply(dB2power(propa.Aatm(f, r.simu, r0.simu, TEMP, RH, PS0)$db),2,mean))
  # Mean HABITAT ATTENUATION over the frequency bandwidth  
  Ahab.dB = power2dB(apply(dB2power(propa.Ahab(f, r.simu, r0.simu, A0)$db),2,mean))
  
  # Estimate L0 over the frequency bandwidth
  L0.simu = L0 - 10*log10(FBW/(F_SIMU[2,ff]-F_SIMU[1,ff]))
  # Estimate L.simu after attenuation of L0.simu
  L.simu =  L0.simu -Ageo.dB -Aatm.dB -Ahab.dB

  # add mean background noise to L.simu
  L.simu_w_noise = power2dB(dB2power(L.simu) + dB2power(L_bkg.exp.mean))
  
  # and finaly set to NA negative L.simu 
  L.simu[L.simu<=0] = NA  
  
  #====================================================================================================#
  #               PLOT Sound Pressure Level (with bin spectrum)
  #====================================================================================================#
  
  DEFAULT_PLOTLY_COLORS=list('rgb(31, 119, 180)', 'rgb(255, 127, 14)',
                             'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
                             'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
                             'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
                             'rgb(188, 189, 34)', 'rgb(23, 190, 207)')
  
  XAXIS_DIST = list(title = "Distances [m]",
                    gridcolor = 'rgb(255,255,255)',
                    showgrid = TRUE,
                    showline = FALSE,
                    showticklabels = TRUE,
                    tickcolor = 'rgb(127,127,127)',
                    ticks = 'outside',
                    zeroline = FALSE)
  
  YAXIS_SPL = list(title = "Sound Level L [dB SPL]",
                   range = c(-5, 90),
                   gridcolor = 'rgb(255,255,255)',
                   showgrid = TRUE,
                   showline = FALSE,
                   showticklabels = TRUE,
                   tickcolor = 'rgb(127,127,127)',
                   ticks = 'outside',
                   zeroline = FALSE)
  # Plot
  p9 <- plot_ly(width = 450, height = 300, x = r.simu, y =L.simu, type = 'scatter', mode = 'lines', name = paste('Propagated signal (model)',sep=""),
                showlegend = TRUE)
  p9 <- add_trace(p9, x = r.exp, y = L.exp, type = 'scatter', mode = 'markers', name= paste('Propagated signal (exp)',sep=""),
                  marker=list(opacity=1,  size=10),
                  showlegend = TRUE)
  p9 <- add_trace(p9, x=r.simu, y = L0.simu - Ageo.dB, type = 'scatter', mode = 'lines', name ='Spreading Losses (model)' ,
                  # marker=list(opacity=1, symbol=I(5), size=5),
                  line = list(dash="dash"),
                  showlegend = TRUE)
  p9 <- add_trace(p9, x=r.simu, y = L0.simu - Ageo.dB - Aatm.dB, type = 'scatter', mode = 'lines', name ='Spreading Losses + Atmospheric attenuation (model)' ,
                  # marker=list(opacity=1, symbol=I(2), size=5),
                  line = list(dash="dot"),
                  showlegend = TRUE)
  p9 <- add_lines(p9, x = r.simu, y = L_bkg.exp.mean, mode = 'lines', name= paste('Ambient sound level',sep=""),
                  line = list(dash="dot"),
                  marker=list(opacity=0.75, symbol=I(2), size=9),
                  showlegend = TRUE)
  p9 <- add_lines(p9, x = r.simu, y = L.simu_w_noise, mode = 'lines', name= paste('Propagated signal (model) + ambient sound (exp)',sep=""),
                  marker=list(opacity=0.75, symbol=I(1), size=6, color='black'),
                  line = list(dash="dot", color='black'),
                  showlegend = TRUE)
  p9 <-layout(p9,
              title = paste(F_SIMU[1,ff],'-',F_SIMU[2,ff],'kHz',sep=""),
              paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
              legend = list(x = 0.1, y = 0.95),
              showlegend = FALSE,
              xaxis = XAXIS_DIST,
              yaxis = YAXIS_SPL)
  
  print(ff)
  fig[[ff]] = list(p9)
}

fig