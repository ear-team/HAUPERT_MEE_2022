# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARIES
#
#=======================================================================================================#

if (!require("plotly")) install.packages("plotly")
if (!require("robustbase")) install.packages("robustbase")
library(plotly)
library(robustbase)

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
# display results ?
PLOT = TRUE
# correction of the background noise
CORRECTION_BKG = TRUE
# BIN size when transforming PSD into histogram
DELTA_FBIN = 0.5 # in kHz

# Select the frequency range
F_RANGE = c(0.5,15)

# Set the right environmental values depending on the habitat
if (grepl("jura", FILENAME_ROOT) == TRUE)
{
  # Temperature
  TEMP = 17.0 # 
  # relative humidity
  RH = 67# 
  # atmospheric pressure in Pa
  PS0 = 87999
  if (grepl("sm4", FILENAME_ROOT) == TRUE) 
  {
    # List of r0 to compute a0 (select only r0 which gives an a0 with R²>0.5)
    LIST_R0 = c(10,20,30,40)
  } else # svantek
  {
    # List of r0 to compute a0 (select only r0 which gives an a0 with R²>0.5)
    LIST_R0 = c(10,20,30,40,50)
  }
  
} else # french guiana
{
  # Temperature
  TEMP = 23.8 # 
  # relative humidity
  RH = 87# 
  # atmospheric pressure in Pa
  PS0 = 1.01340e5
  if (grepl("sm4", FILENAME_ROOT) == TRUE) 
  {
    # List of r0 to compute a0 (select only r0 which gives an a0 with R²>0.5)
    LIST_R0 = c(10,20,30)
  } else # svantek
  {
    # List of r0 to compute a0 (select only r0 which gives an a0 with R²>0.5)
    LIST_R0 = c(10,20,30,40)
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
load(paste(FILE_DIR,FILENAME_ROOT, '_average', '.Rdata', sep=""))

#====================================================================================================#
#               convert spectrum into specbins
#====================================================================================================#

f         = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$f
P_bkg     = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$s
P         = specbin(sig.PSD.mean, FREQUENCY, DELTA_FBIN)$s

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
    index = P.select < P_bkg.select         
    P.select[index] = NA      
    # subtract the noise level to the original sound level
    P.select = P.select-P_bkg.select 
  }

  # spherical or spreading loss attenuation : MATRIX [N*M]
  Ageo.dB = propa.Ageo(r=r, r0=r0)$db
  Ageo.dB = matrix(rep(Ageo.dB,each=length(f)),nrow=length(f))
  # atmospheric attenuation 
  Aatm.dB = propa.Aatm(f=f, r=r, r0=r0, t=TEMP, rh=RH, pa=PS0)$db
  
  # Convert psd into dB scale
  L = power2dB(P.select)
  # create a matrix with L0 (for distance r0)
  L0 = t(matrix(rep(L[,1],each=length(r)),nrow=length(r)))
  
  # excess attenuation : MATRIX [N*M]
  EA =  L0 -L - Ageo.dB  - Aatm.dB
  
  # create matrix of r and f in order to be able to divide the matrix EA
  r_mat = matrix(rep(r,each=length(f)),nrow=length(f))
  f_mat = t(matrix(rep(f,each=length(r)),nrow=length(r)))*1000
  
  if (PLOT==TRUE)
  {
    # plot the average (EA) for distances vs frequency
    if (first == 1)
    {
      EA_per_d = power2dB(rowMeans(dB2power(EA/(r_mat-r0)), na.rm=TRUE))
      pl1 <- plot_ly(x=f,y=EA_per_d, type='scatter', mode = 'markers', name = paste('r0=',r0,'m',sep=''),showlegend = TRUE, marker=list(opacity=0.75, symbol =I(10-dd), color= dd))
      EA_per_d.mean = matrix(EA_per_d)
      first = 0
    }
    else
    {
      
      EA_per_d = power2dB(rowMeans(dB2power(EA/(r_mat-r0)), na.rm=TRUE))
      pl1 = add_trace(pl1,x=f,y=EA_per_d, type='scatter', mode = 'markers', name = paste('r0=',r0,'m', sep=''),showlegend = TRUE, marker=list(opacity=0.75, symbol =I(10-dd), color= dd))
      pl1 = layout(pl1,paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)',
                   xaxis = list(title = "Frequency [kHz]"), #,type = "log"),
                   yaxis = list(title = "EA attenuation [dB/m] "))
      EA_per_d.mean = cbind(EA_per_d.mean, EA_per_d)
    }
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
x= as.vector(x)
data = data.frame(f=x, EA=y)
model <- lmrob(EA~f+0, data, setting = "KS2014", na.rm=TRUE)
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
  sprintf('a0 = %2.3f',a0)
  r.squared = summary(model)$r.squared
  sprintf('r² = %1.2f',r.squared)
  pvalue = "  p<0.001"
  sprintf('p = %1.4f',summary(model)$coefficients[4])
}

predicted = data.frame(f =f)
predicted$EA = predict(model,predicted)
predicted$ci_binf <- predict(model,predicted, interval="confidence",level = 0.95)[,2]
predicted$ci_bsup <- predict(model,predicted, interval="confidence",level = 0.95)[,3]
predicted$pi_binf <- predict(model,predicted, interval="prediction",level = 0.95)[,2]
predicted$pi_bsup <- predict(model,predicted, interval="prediction",level = 0.95)[,3]

if (PLOT==TRUE)
{
  pl1 <- add_trace(pl1, x=f,y=power2dB(rowMeans(dB2power(EA_per_d.mean),na.rm=TRUE)), type='scatter', mode = 'markers', marker=list(opacity=0.75, symbol =I(1), color= 'black'), name='avg')
  
  # pl1 = add_trace(pl1, x=f,y=predicted$EA, type='scatter', mode='lines',name='fit', line=list(dash='solid',color='black',size=3),marker=list(opacity=0), opacity=1)
  # pl1 = add_trace(pl1, x=f,y=predicted$ci_binf, type='scatter', mode='lines',name='ci_binf', line=list(dash='dash',color='red'),marker=list(opacity=0), opacity=1)
  # pl1 = add_trace(pl1, x=f,y=predicted$ci_bsup, type='scatter', mode='lines',name='ci_bsup', line=list(dash='dot',color='red'),marker=list(opacity=0), opacity=1)
  # pl1 = add_trace(pl1, x=f,y=predicted$pi_binf, type='scatter', mode='lines',name='pi_binf', line=list(dash='dash',color='navy'),marker=list(opacity=0), opacity=1)
  # pl1 = add_trace(pl1, x=f,y=predicted$pi_bsup, type='scatter', mode='lines',name='pi_bsup', line=list(dash='dot',color='navy'),marker=list(opacity=0), opacity=1)
  # pl1 = add_annotations(pl1, x=5, y=0.65, text =paste("a0=",round(a0*1000)/1000," dB/kHz/m", sep=""), showarrow = FALSE,marker=list(opacity=0), opacity=1)
  # pl1 = add_annotations(pl1, x=4.8, y=0.57, text =paste("r²=",round(r.squared*100)/100, pvalue, sep=""), showarrow = FALSE,marker=list(opacity=0), opacity=1)
  pl1 = layout(pl1,
               autosize = F, width = 450, height = 300,
               paper_bgcolor='rgb(255,255,255)', 
               plot_bgcolor='rgb(229,229,229)',
               xaxis = list(title = "Frequency [kHz]"), #,type = "log"),
               yaxis = list(title = "EA/(r-r0) [dB/m] ", range=c(-0.5, 0.7)
               ),
               showlegend = TRUE
  )
  pl1
}





