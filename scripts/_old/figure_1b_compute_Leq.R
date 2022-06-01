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
FILENAME_ROOT = "jura_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/"
# display results ?
PLOT = TRUE
# correction recorder
CORRECTION_RECORDER = FALSE

#====================================================================================================#
#                                       LOAD Data
#
#               variables : 
#                   - filename, 
#                   - DISTANCES, FREQUENCY,
#                   - CHANNEL, Nrep,
#                   - NFFT, fs_rec,
#                   - S, G, bit, VADC, P_REF, 
#                   - bkg.PSD.mean, bkg.PSD.std, sig.PSD.mean, sig.PSD.std,
#====================================================================================================#
load(paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep=""))

#===========================================
#           Compute Leq
#===========================================

sig.PSD.Leq.mean = rep(NA,length(DISTANCES))
sig.PSD.Leq.std  = rep(NA,length(DISTANCES))
bkg.PSD.Leq.mean= rep(NA,length(DISTANCES))
bkg.PSD.Leq.std = rep(NA,length(DISTANCES))
for (ii in 1:length(DISTANCES))
{
  sig.PSD.Leq.mean[ii] = psd2leq(sig.PSD.mean[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  sig.PSD.Leq.std[ii] = psd2leq(sig.PSD.std[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  
  bkg.PSD.Leq.mean[ii] = psd2leq(bkg.PSD.mean[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
  bkg.PSD.Leq.std[ii] = psd2leq(bkg.PSD.std[,ii], gain=G, sensitivity=S, bit=bit, Vadc=VADC, pRef = P_REF)
}

#====================================================================================================#
#
#         APPLY THE SM4 RECORDER CORRECTION TO FLATTEN (as possible) THE FREQUENCY RESPONSE
#
#====================================================================================================#
if ((CORRECTION_RECORDER == TRUE) && (grepl("sm4", FILENAME_ROOT) == TRUE))
{
  # load the GENERIC frequency response of the SM4. 
  # The gain to correct the frequency response of the SM4 was measured on a single SM4 recorder, assuming that all SM4 have the exact same frequency response
  # It would be better to adapt the correction for each SM4 if possible.
  load('./toolbox/SM4_gain.Rdata')
  FREQ_GAIN_CORR = SM4.G[,1]
  GAIN_CORR = SM4.G[,2] 
  # linear interpolation of SM4.G in order to match the frequency bin
  # SM4.gain contains the gain to add to the result obtained with a SM4 in order to correct as much as possible the frequency response of the SM4.
  SM4.gain = -power2dB(approx(FREQ_GAIN_CORR, dB2power(GAIN_CORR), FREQUENCY, rule=2)$y)
  # mean gain
  SM4.gain.mean = power2dB(mean(dB2power(SM4.gain)))
  
  # correct the gain
  sig.PSD.Leq.mean = sig.PSD.Leq.mean - SM4.gain.mean
  sig.PSD.Leq.std = sig.PSD.Leq.std - SM4.gain.mean
  
  bkg.PSD.Leq.mean = bkg.PSD.Leq.mean - SM4.gain.mean
  bkg.PSD.Leq.std = bkg.PSD.Leq.std - SM4.gain.mean
} 


# print the Leq @1m
sprintf('White noise Leq = %2.1f at 1m',sig.PSD.Leq.mean[DISTANCES==1])  
# print the average background Leq 
sprintf('Mean Background Leq = %2.1f for all distances',power2dB(mean(dB2power(bkg.PSD.Leq.mean))))

#====================================================================================================#
#                Plot Data
#====================================================================================================#
if (PLOT == TRUE)
{
  if (!require("plotly")) install.packages("plotly")
  library(plotly)
 
  BRIGHT_COLORS=list('rgba(35, 130, 200)', #blue     1
                     'rgba(255, 135, 16)', #orange   2
                     'rgba(50, 180, 50)',  # green   3
                     'rgba(235, 45, 45)',  # red     4
                     'rgba(165, 120, 210)', # purple 5
                     'rgba(200, 200, 40)', # yellow  6
                     'rgba(25, 170, 190)', # blue green 7
                     'rgba(240, 135, 220)', # pink   8
                     'rgba(155, 100, 90)',  # brown  9
                     'rgba(150, 150, 150)', # grey    10
                     'rgba(233, 115, 94)', # brick    11
                     'rgba(150, 220, 50)' # greenish    12
                     )
  
  DARK_COLORS=list('rgba(18, 65, 100)', #blue
                    'rgba(128, 70, 8)', #orange
                    'rgba(25, 90, 25)',  # green
                    'rgba(120, 25, 25)',  # red
                    'rgba(85, 60, 105)', # purple
                    'rgba(100, 100, 20)', # yellow
                    'rgba(13, 90, 95)', # blue green
                    'rgba(120, 70, 110)', # pink
                    'rgba(80, 50, 45)',  # brown
                    'rgba(75, 75, 75)', # grey
                   'rgba(115, 57, 45)', # brick    11
                   'rgba(75, 110, 25)' # greenish    12
                    )
  
  XAXIS_DIST = list(title = "Distances [m]",
                    gridcolor = 'rgb(255,255,255)',
                    showgrid = TRUE,
                    showline = FALSE,
                    showticklabels = TRUE,
                    tickcolor = 'rgb(127,127,127)',
                    ticks = 'outside',
                    zeroline = FALSE)
  
  YAXIS_LEQ = list(title = "SPL [dB]",
                   range = c(0, 100),
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
  sig_mean = sig.PSD.Leq.mean
  sig_std = sig.PSD.Leq.std
  bkg_mean = bkg.PSD.Leq.mean
  bkg_std = bkg.PSD.Leq.std
  
  sig_up = power2dB(dB2power(sig_mean) + dB2power(sig_std))
  #sig_lo = power2dB(dB2power(sig_mean) - dB2power(sig_std))
  sig_lo = sig_mean - (sig_up-sig_mean)

  bkg_up = power2dB(dB2power(bkg_mean) + dB2power(bkg_std))
  #bkg_lo = power2dB(dB2power(bkg_mean) - dB2power(bkg_std))
  bkg_lo = bkg_mean - (bkg_up-bkg_mean)
  
  # # spherical attenuation (select the reference at 1m, not 0.5m)
  # sph_att = 20*log10(DISTANCES/DISTANCES[10]) 
  # theo_att = sig_mean[10] - sph_att
  
  p <- plot_ly(x = DISTANCES, y =sig_mean , type = 'scatter', mode = 'lines', name ='White noise',connectgaps = TRUE, 
               line = list(color=gsub(")", " ,10)", BRIGHT_COLORS[10])),showlegend = TRUE)
  p <- add_trace(p, y = sig_up, type = 'scatter', mode = 'lines', connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.10)", BRIGHT_COLORS[10]), dash="dot"), 
                 showlegend = FALSE)
  p <- add_trace(p, y = sig_lo, type = 'scatter', mode = 'lines',
                 fill = 'tonexty', fillcolor= list(color=gsub(")", " ,10)", BRIGHT_COLORS[10])), connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.10)", BRIGHT_COLORS[10]), dash="dot"),
                 showlegend = FALSE)
  
  p <- add_trace(p, y =bkg_mean , type = 'scatter', mode = 'lines', name = 'Ambient sound', connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,10)", DARK_COLORS[10])),showlegend = TRUE)
  p <- add_trace(p, y = bkg_up, type = 'scatter', mode = 'lines', connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.10)", DARK_COLORS[10]), dash="dot"), 
                 showlegend = FALSE)
  p <- add_trace(p, y = bkg_lo, type = 'scatter', mode = 'lines',
                 fill = 'tonexty', fillcolor= list(color=gsub(")", " ,10)", DARK_COLORS[10])), connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.10)", DARK_COLORS[10]), dash="dot"),
                 showlegend = FALSE)
  
  # p <- add_trace(p, y = theo_att, type = 'scatter', mode = 'lines', name = "Spreading loss [theory]",connectgaps = TRUE, 
  #                line = list(color=gsub(")", " ,1)", '#000000')), 
  #                showlegend = TRUE)
  
  p3 <- layout(p, 
               autosize = F, width = 450, height = 250,
               # title = 'Equivalent Continuous Sound level (Leq)',
               paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
               legend = list(x = 0.6, y = 0.05),
               xaxis = XAXIS_DIST,
               yaxis = YAXIS_LEQ)
  
  p3

}

