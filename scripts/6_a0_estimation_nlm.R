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
# Save results
SAVE = TRUE
# display results ?
PLOT = TRUE
# correction of the background noise
CORRECTION_BKG = TRUE
# BIN size when transforming PSD into histogram
DELTA_FBIN = 0.5 # in kHz
# Select the minimum distance that will be used for all the calculation : 10m
DISTANCE_MIN=10
# Select the frequency range
F_RANGE = c(0.5,15)

# Set the right environmental values depending on the habitat
if (grepl("jura", FILENAME_ROOT) == TRUE)
{
  # Temperature
  TEMP = 17 # 
  # relative humidity
  RH = 67# 
  # atmospheric pressure in Pa
  PS0 = 87999
} else # french guiana
{
  # Temperature
  TEMP = 23.8 # 
  # relative humidity
  RH = 87# 
  # atmospheric pressure in Pa
  PS0 = 1.01340e5
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

#====================================================================================================#
#               convert spectrum into specbins
#====================================================================================================#

f         = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$f
P_bkg     = specbin(bkg.PSD.mean, FREQUENCY, DELTA_FBIN)$s
P_bkg.std = specbin(bkg.PSD.std, FREQUENCY, DELTA_FBIN)$s
P         = specbin(sig.PSD.mean, FREQUENCY, DELTA_FBIN)$s
P.std     = specbin(sig.PSD.std, FREQUENCY, DELTA_FBIN)$s

# keep only the frequencies corresponding to the frequency range (F_RANGE)
index_select = f>=F_RANGE[1] & f<=F_RANGE[2]
f         = f[index_select]
P_bkg     = P_bkg[index_select,]
P_bkg.std = P_bkg.std[index_select,]
P         = P[index_select,]
P.std     = P.std[index_select,]

#====================================================================================================#
#               Convert psd (P) into dB sound pressure level (L)
#====================================================================================================#
# create a vector with the index of the selected distances
DISTANCE_SELECT = DISTANCES>=DISTANCE_MIN
# select distances use for the calculation : VECTOR of size M                                   
r = DISTANCES[DISTANCE_SELECT]           
# select the initial distance
r0 = r[1]

#### ENERGY (no unit)
# P :  energy of the propagated signal : MATRIX [N*M]
P  = P[,DISTANCE_SELECT]           
# P.std : standard deviation of the energy of the propagated signal
P.std = P.std[,DISTANCE_SELECT]
# P_bkg : energy of the ambient sound 
P_bkg = P_bkg[,DISTANCE_SELECT]  

if (CORRECTION_BKG == TRUE)
{
  #### CORRECTED SIGNAL AFTER AMBIENT SOUND SUBRACTION
  # Set to NA the data below the ambient sound in order to avoid using them for the calculation
  index = P< P_bkg         
  P[index] = NA      
  # subtract the noise level to the original sound level
  P= P-P_bkg 
}

#### signal (amplitude)
p = sqrt(P)
p0 = p[,1]

# spherical or spreading loss attenuation : MATRIX [N*M]
Ageo.dB = propa.Ageo(r=r, r0=r0)$db
Ageo.dB = matrix(rep(Ageo.dB,each=length(f)),nrow=length(f))
# atmospheric attenuation 
Aatm.dB = propa.Aatm(f=f, r=r, r0=r0, t=TEMP, rh=RH, pa=PS0)$db

# create a dataframe with attenuation parameters to use for the model
df.model <- data.frame(r=as.vector(t(replicate(length(f), r))),
                       f = rep(f, length(r)),
                       Ageo.dB = as.vector(Ageo.dB),
                       Aatm.dB = as.vector(Aatm.dB),
                       r0= rep(r0,length(f)*length(r)),
                       p = as.vector(p),
                       p0= rep(p0, length(r)),
                       pr =as.vector(p * matrix(rep(r,each=length(f)),nrow=length(f))),
                       p0r0 = rep(p0 *r0,length(r)), 
                       t = rep(TEMP,length(f)*length(r)),
                       rh= rep(RH,length(f)*length(r)),
                       pa= rep(PS0,length(f)*length(r)))

#********************************************************************************#
# Attenuation law : Exponential law with frequency dependency
#********************************************************************************#
#********  atmospheric attenuation
coeff.Aatm <- function(f, t=20, rh=60, pa=101325)
{   
  " 
  Partially from http://www.sengpielaudio.com/AirdampingFormula.htm
  "
  
  pr = 101.325e3 # reference ambient atmospheric pressure: 101.325 kPa
  To1 = 273.16 # triple-point isotherm temp: 273.16 K
  To = 293.15 #  reference temperature in K: 293.15 K (20°C)
  t = t+273.15 # celcius to farenheit
  f = f*1000 #  kHz to Hz
  
  psat = pr*10**(-6.8346 * (To1/t)**1.261 + 4.6151) #saturation vapor pressure equals
  h = rh * (psat / pa) # molar concentration of water vapor, as a percentage
  frO = (pa / pr) * (24 + 4.04e4 * h * ((0.02 + h) / (0.391 + h))) # oxygen relaxation frequency
  frN = (pa / pr) * sqrt(t / To) * (9 + 280 * h * exp(-4.170*((t/To)**(-1/3) -1))) # nitrogen relaxation frequency
  
  z = 0.1068 * exp (-3352/t) / (frN+f**2 /frN)
  y = (t/To)**(-5/2) * (0.01275 * exp(-2239.1/t) * 1/(frO+f**2/frO) + z)
  
  coef.Aatm = f**2 * ((1.84e-11 * 1/(pa/pr) * sqrt(t/To)) + y)

  return (coef.Aatm)
}

# Full attenuation model (exponential law) ready to be linearized with log
func_exp <- function(f, r, r0, p0, t, rh, pa, a0)
{
  # atmospheric attenuation coefficient
  b0 = coeff.Aatm(f, t, rh, pa)
  # full attenuation model 
  pr = p0*r0 * exp(-(b0 + a0*f) * (r-r0))
  return (pr)
}

# nonlinear least square regression
# take the log of the exponential function in order to linearize the expression.
model.nonlinear = nls(log(pr) ~ log(func_exp(f, r, r0, p0, t, rh, pa, a0)), df.model, start = list(a0=0.001), control=list(maxiter = 500), na.action=na.exclude)

# test if the p-value < 0.001
if (summary(model.nonlinear)$coefficients[4] < 1e-3)
{
  a0 = round(coef(model.nonlinear)['a0'] * 10000)/10000
  a0_dB = 20*log10(exp(1)) * a0
  pvalue = summary(model.nonlinear)$coefficients[4] 
} else
{
  a0 = NA
  a0_dB = NA
  pvalue = summary(model.nonlinear)$coefficients[4] 
}

sprintf('a0 = %2.4f / a0_dB = %2.4f / pvalue = %2.5f',a0, a0_dB, pvalue)
