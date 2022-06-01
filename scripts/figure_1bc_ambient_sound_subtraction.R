# clear workspace
rm(list = ls())

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
FILENAME_ROOT = "guiana_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/" 
# Save results
SAVE = TRUE
# display results ?
PLOT = TRUE
# SM4 gain correction
CORRECTION_RECORDER = FALSE
# BIN size when transforming PSD into histogram
DELTA_FBIN = 0.5 # in kHz
# Select the minimum distance that will be used for all the calculation : 10m
DISTANCE_MIN= 10

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
load(paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep=""))

#====================================================================================================#
#               convert spectrum into specbins
#====================================================================================================#

f     = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$f
P_bkg = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$s
P_bkg.std = specbin(bkg.PSD.std, FREQUENCY, DELTA_FBIN)$s
P = specbin(sig.PSD.mean, FREQUENCY, DELTA_FBIN)$s
P.std = specbin(sig.PSD.std, FREQUENCY, DELTA_FBIN)$s

#====================================================================================================#
#               Convert psd (P) into dB sound pressure level (L)
#====================================================================================================#
# create a vector with the index of the selected distances
DISTANCE_SELECT = DISTANCES>=DISTANCE_MIN
# select distances use for the calculation : VECTOR of size M                                   
r = DISTANCES[DISTANCE_SELECT]                                                            

#### ENERGY (no unit)
# P :  energy of the propagated signal : MATRIX [N*M]
P  = P[,DISTANCE_SELECT]           
# P.std : standard deviation of the energy of the propagated signal
P.std = P.std[,DISTANCE_SELECT]
# P_bkg : energy of the ambient sound 
P_bkg = P_bkg[,DISTANCE_SELECT]  

#### SOUND PRESSURE LEVEL (dB SPL)
#### ORIGINAL SIGNAL
# Transform the energy into dB SPL : MATRIX [N*M]
# get the original sound level (with noise) L 
L = psd2dBSPL(P, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 
# get the noise level
L_bkg = psd2dBSPL(P_bkg, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 

#### CORRECTED SIGNAL AFTER AMBIENT SOUND SUBRACTION
# Set to NA the data below the ambient sound in order to avoid using them for the calculation
index = P< P_bkg         
P[index] = NA      
# subtract the noise level to the original sound level
P_corr= P-P_bkg 
L_corr = psd2dBSPL(P_corr, gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef=P_REF) 


#====================================================================================================#
#               Correction of the SM4 frequency response in dB SPL
#====================================================================================================#

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
  
  # create a matrix in order to be able to perform direct subtraction on the matrix of L
  SM4.gain.mat = replicate(dim(L)[2],SM4.gain)
  
  # correct the gain
  L = L - SM4.gain.mat
  L_corr = L_corr - SM4.gain.mat
  L_bkg = L_bkg - SM4.gain.mat
} 

#====================================================================================================#
#               PLOT Sound Pressure Level (with bin spectrum)
#====================================================================================================#

if (PLOT==TRUE)
{
  if (!require("plotly")) install.packages("plotly")
  library(plotly)  
  
  DEFAULT_PLOTLY_COLORS=list('rgb(31, 119, 180)', 'rgb(255, 127, 14)',
                             'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
                             'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
                             'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
                             'rgb(188, 189, 34)', 'rgb(23, 190, 207)')
  
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
  
  #====================================================================================================#
  #               Plot sraw pectrum for each distance
  #====================================================================================================#

  fig_spectrum <- plot_ly(autosize = F, width = 450, height = 350,
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
  
  fig_LdBSPL =layout(fig_LdBSPL, autosize = F, width = 450, height = 350,
                     paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
                     xaxis2 = list(title = "Distances [m]"),
                     yaxis = XAXIS_FREQ)
  show(fig_LdBSPL)
}


