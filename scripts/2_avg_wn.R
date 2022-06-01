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
source("./toolbox/toolbox_propagation_wav_extraction.R", chdir=T)

#=======================================================================================================#
#
#                                           SET GLOBAL PARAMETERS
#
#=======================================================================================================#

# filename root of the data
FILENAME_ROOT = "guiana_sm4_wn" # guiana_svantek_wn guiana_sm4_wn jura_svantek_wn jura_sm4_wn
FILE_DIR = "../data/psd/"
# Display ?
PLOT = TRUE
# save ?
SAVE = TRUE
# number of repetition per distance
Nrep = 3

#====================================================================================================#
#                                 LOAD Data
#
#               variables : 
#                   - filename, 
#                   - DISTANCES, FREQUENCY,
#                   - CHANNEL, Nrep,
#                   - NFFT, fs_rec,
#                   - S, G, bit, VADC, P_REF, 
#                   - bkg.PSD, sig.PSD,
#====================================================================================================#

for (ii in 1:Nrep)
{
  # Load data saved with 1_extraction_wn_XXX_YYY.R
  load(paste(FILE_DIR,FILENAME_ROOT,'_timeline',ii,'.Rdata', sep=""))
  
  if (ii ==1)
  {
    bkg.PSD.mat= bkg.PSD
    sig.PSD.mat = sig.PSD
  }
  else
  {
    bkg.PSD.mat= cbind(bkg.PSD.mat,bkg.PSD)
    sig.PSD.mat = cbind(sig.PSD.mat,sig.PSD)
  } 
} 

#====================================================================================================#
#                AVERAGE Data
#====================================================================================================#

bkg.PSD.mean = matrix(NA,nrow=NFFT/2,ncol=length(DISTANCES))
bkg.PSD.std = matrix(NA,nrow=NFFT/2,ncol=length(DISTANCES))
sig.PSD.mean = matrix(NA,nrow=NFFT/2,ncol=length(DISTANCES))
sig.PSD.std = matrix(NA,nrow=NFFT/2,ncol=length(DISTANCES))

for (ii in 1:length(DISTANCES))
{
  idx = seq(ii,Nrep*length(DISTANCES),length(DISTANCES))
  
  # PSD
  bkg.PSD.mean[,ii] = rowMeans(bkg.PSD.mat[,idx])
  bkg.PSD.std[,ii] = apply(bkg.PSD.mat[,idx], 1, sd)
  
  sig.PSD.mean[,ii] = rowMeans(sig.PSD.mat[,idx])
  sig.PSD.std[,ii] = apply(sig.PSD.mat[,idx], 1, sd)
}

########## Save data
if (SAVE==TRUE)
{
  save(filename, 
       DISTANCES, FREQUENCY,
       CHANNEL, Nrep,
       NFFT, fs_rec,
       S, G, bit, VADC, P_REF, 
       bkg.PSD.mean, bkg.PSD.std,sig.PSD.mean, sig.PSD.std,
       file=paste(FILENAME_ROOT, '_average','.Rdata', sep=""))
}

#====================================================================================================#
#                Plot Data
#====================================================================================================#

# Load data previously saved
load(paste(FILENAME_ROOT, '_average', '.Rdata', sep=""))

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
  
  YAXIS_SPL = list(title = "Power Density Spectrum [dB]",
                   range = c(-80, 5),
                   gridcolor = 'rgb(255,255,255)',
                   showgrid = TRUE,
                   showline = FALSE,
                   showticklabels = TRUE,
                   tickcolor = 'rgb(127,127,127)',
                   ticks = 'outside',
                   zeroline = FALSE)

  #====================================================================================================#
  #               PLOT Sound Pressure Level (with the entire spectrum)
  #====================================================================================================#
  # convert PSD into decibel (dB)
  sig_mean = power2dB(sig.PSD.mean/max(abs(sig.PSD.mean))) 
  sig_std  = power2dB(sig.PSD.std/max(abs(sig.PSD.std)))  
  # Create the upper bound of the signal which is +1 std
  sig_up = power2dB(dB2power(sig_mean) + dB2power(sig_std))
  # for display, pseudo std in dB. 
  sig_lo = sig_mean - (sig_up-sig_mean)
  # real values
  #sig_lo = power2dB(dB2power(sig_mean) - dB2power(sig_std))
  
  p <- plot_ly(x = FREQUENCY, y =sig_mean[,1] , type = 'scatter', mode = 'lines', name = paste(DISTANCES[1],'m'),connectgaps = TRUE, 
               line = list(color=gsub(")", " ,1)", BRIGHT_COLORS[1])),showlegend = TRUE)
  p <- add_trace(p, y = sig_up[,1], connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[1]), dash="dot"),
                 showlegend = FALSE)
  p <- add_trace(p,  y = sig_lo[,1], fill = 'tonexty', fillcolor= list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[1])), connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[1]), dash="dot"),
                 showlegend = FALSE)
  
  for (ii in 2:length(DISTANCES))
  {
    p <- add_trace(p,y =sig_mean[,ii] , type = 'scatter', mode = 'lines', name = paste(DISTANCES[ii],'m'),connectgaps = TRUE, 
                 line = list(color=gsub(")", " ,1)", BRIGHT_COLORS[ii])),showlegend = TRUE)
    p <- add_trace(p, y = sig_up[,ii],connectgaps = TRUE,
                   line = list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[ii]), dash="dot"),
                   showlegend = FALSE)
    p <- add_trace(p, y = sig_lo[,ii],
                   fill = 'tonexty', fillcolor= list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[ii])), connectgaps = TRUE,
                   line = list(color=gsub(")", " ,0.1)", BRIGHT_COLORS[ii]), dash="dot"),
                   showlegend = FALSE)
  }
  p <-layout(p, 
         title = 'Power Density Spectrum of white noise after different distances of propagation',
         paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(225,225,225)',
         legend = list(orientation = 'v', 
                       automargin = TRUE,
                       font = list(size = 12)),
         xaxis = XAXIS_FREQ,
         yaxis = YAXIS_SPL)
  p
}

