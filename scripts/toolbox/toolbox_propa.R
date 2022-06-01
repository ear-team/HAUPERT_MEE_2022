###**********************************************************************###
# Functions to convert wav file into Sound Pressure Level dB SPL
# Functions to calculate outdoor propagation in natural habitat
# Functions to estimate the detection or active distance
#
# Author : Sylvain HAUPERT, MNHN, sylvain.haupert@mnhn.fr
# Date : 14/12/2021
# Licence : MIT
####**********************************************************************###

if (!require("signal")) install.packages("signal")
library(signal)

# Constant : minimum value possible
MIN = 2.2250738585072014e-308

#####################
dB2amplitude <- function (db)
{
  return (10**(db/20))
}

dB2power <- function (db)
{
  return (10**(db/10))
}

amplitude2dB <- function (s)
{
  return (20*log10(s))
}

power2dB <- function (s)
{
  return (10*log10(s))
}

####################

##### wav2volt
wav2volt <- function (wave, bit=16, Vadc=2) 
{ 
  # convert in Volt. Vadc is 2Vpp [peak to peak] for the SM4
  volt = (wave/(2^(bit))) * Vadc
  return (volt)
}

#####  volt2pressure
volt2pressure <- function (volt, gain, sensitivity=-35) 
{  
  #volt amplitude to instantaneous sound pressure (Pa)
  coef = 1/10^(sensitivity/20) 
  p = volt * coef / 10^(gain/20)
  
  return(p)
}

#####  wav2pressure
wav2pressure <- function (wave, gain, sensitivity=-35, bit=16, Vadc=2) 
{
  # Vadc => 2Vpp [peak to peak] for the SM4
  p = volt2pressure(wav2volt(wave, bit, Vadc), gain, sensitivity) 
  
  return (p)
}

###################################"

#####  pressure2dBSPL
pressure2dBSPL <- function(p, pRef=20e-6)
{
  # if p <=0 set to MIN
  p[p==0] = MIN
  # Take the log of the ratio pressure/pRef
  L = 20*log10(p/pRef) 
  return (L)
}

#####  dBSPL2pressure
dBSPL2pressure <- function(L, pRef=20e-6)
{
  # dB SPL to pressure
  p = 10**(L/20)*pRef
  return (p)
}                

#####  wav2dBSPL
wav2dBSPL <- function (wave, gain, sensitivity=-35, bit=16, Vadc=2, pRef=20e-6) 
{
  p = wav2pressure(wave, gain, sensitivity, bit, Vadc) 
  L = pressure2dBSPL(p, pRef)
  return (L)
}

#####  psd2dBSPL
psd2dBSPL <- function (P, gain, sensitivity=-35, bit=16, Vadc=2, pRef=20e-6) 
{
  # convert power (energy) to amplitude
  s = sqrt(P)
  # convert amplitude to dB sPL
  L = wav2dBSPL(s, gain, sensitivity, bit, Vadc, pRef)
  return (L)
}

########################

#####  wav2Leq
wav2leq <- function (wave, f, gain, dt=1, sensitivity=-35, bit=16, Vadc=2, pRef= 20e-6) 
{  
  # convert into pressure
  p = wav2pressure(wave, gain, sensitivity, bit, Vadc) 
  
  # p to Leq (Equivalent Continuous Sound level)
  Leq = pressure2leq(p, f, dt, pRef) 
  
  return(Leq)
}

#####  pressure2Leq
pressure2leq <- function (p, f, dt=1, pRef = 20e-6) 
{  
  # wav to RMS
  dN = dt*f # integration period in number of points
  N_RMS = floor(length(p)/dN)
  
  p_RMS = matrix(NA,nrow=1,ncol=N_RMS)
  
  for (ii in 1:N_RMS)
  {
    p_mean = mean(p[(1+(ii-1)*dN):(ii*dN)]^2)
    p_RMS[ii] = sqrt(p_mean)
  }
  
  # if p_RMS ==0 set to MIN
  p_RMS[p_RMS==0] = MIN
  
  # RMS to Leq (Equivalent Continuous Sound level)
  Leq = 20*log10(p_RMS/pRef)
  
  return(Leq)
}

#####  psdLeq
psd2leq <- function (psd, gain, sensitivity=-35, bit=16, Vadc=2, pRef=20e-6) 
{  
  # convert P (amplitude²) to pressure²
  P = wav2pressure (sqrt(psd), gain, sensitivity, bit, Vadc)**2
  
  # if P ==0 set to MIN
  P[P==0] = MIN
  
  # Energy (pressure^2) to Leq 
  # => Leq (Equivalent Continuous Sound level) if the sum is performed on the whole PSD
  leq = 10*log10(sum(P)/pRef**2)
  return (leq)
}

########################

#******** Convert a full PSD into bins PSD with a resolution = res
specbin <- function (spec, f, res = 0.5)
{
  " 
  Convert a full PSD into bins PSD with a resolution = res
  
  INPUTS :
  spec : power spectrum density (amplitude²) [VECTOR or MATRIX]
  f : frequency vector in kHz [VECTOR]
  res : frequency resolution of each frequency bin in kHz
  
  OUTPUT
  return a list with
  - freq_bin : the frequency bin 
  - spec_bin : the bins PSD
  "
  # create the frequency bins vector
  vbin = round(seq(from=f[1]*1000, to=f[length(f)]*1000,by=res*1000))/1000
  freq_bin = round(vbin[2:length(vbin)]*1000 - (res/2)*1000)/1000
  for (jj in 1:(length(vbin)-1))
  {
    # Selected frequencies
    fselect = (f>=vbin[jj] & f<vbin[jj+1])
    
    # ratio number of f / number of bins
    ratio = length(f)/length(vbin)
    
    # sum of each point of the PSD (10*log[sum (PSD) / pref²]) corresponding to the frequency range fselect
    if (jj==1)
    {
      if (is.null(dim(spec)))
      {
        spec_bin= mean(spec[fselect]) * ratio
      }
      else
      {
        spec_bin= apply(spec[fselect,], 2, mean) * ratio
      }
    }
    else
    {
      if (is.null(dim(spec)))
      {
        spec_bin= rbind(spec_bin,mean(spec[fselect])* ratio) 
      }
      else
      {
        spec_bin= rbind(spec_bin,apply(spec[fselect,], 2, mean)* ratio) 
      }      
    }
  }

  return (list(f = freq_bin, s = spec_bin))
}

########################

#******** geometrical (or spherical) attenuation
propa.Ageo <- function(r,r0) 
{
  " 
  Get the attenuation due to spreading loss (also known as spherical or geometrical attenuation)
  
  INPUTS :
  r : propagation distances in m [SCALAR or VECTOR]
  r0 : reference distance in m [SCALAR]
  
  OUTPUT
  return a list with 
  - factor : the geometrical (or spherical) attenuation factor of an acoustic pressure 
  => Multiply this value with the effective reference acoustic pressure p0 measured at r0 to estimate the pressure after attenuation
  - db : the geometrical (or spherical) attenuation of an acoustic pressure in dB 
  => subtract this value to the reference sound pressure level (dB SPL) 
  "
  Ageo = r0/r
  Ageo.dB  = -20*log10(Ageo)
  return (list(factor= Ageo, db= Ageo.dB))
}

#********  atmospheric attenuation
propa.Aatm <- function(f, r=NA, r0=NA, t=20, rh=60, pa=101325)
{   
  " 
  Get the atmospheric attenuation factor and in dB and coefficients in Neper/m and dB/m  
  
  INPUTS :
  f: frequency in kHz [SCALAR]
  r : propagation distances in m [SCALAR or VECTOR]
  r0 : reference distance in m [SCALAR]
  t: temperature in °C [SCALAR]
  rh: relative humidity in % [SCALAR]
  pa: atmospheric pressure in Pa [SCALAR]
  
  OUTPUT
  return a list with 
  - factor : atmospheric attenuation factor depending on the frequency, the temperature, the atmospheric pressure, the relative humidity and the distance [MATRIX] 
  => Multiply this value with the effective reference acoustic pressure p0 measured at r0 to estimate the pressure after attenuation
  - db : atmospheric attenuation in dB depending on the frequency, the temperature, the atmospheric pressure, the relative humidity and the distance [MATRIX] 
  => subtract this value to the reference sound pressure level (dB SPL) for each frequency and distance
  - coef : atmospheric attenuation coefficient in Neper/m
  - coeff.db : atmospheric attenuation coefficient in dB/m
  
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
  # atmospheric attenuation coefficient in Neper/m 
  coef.Aatm = f**2 * ((1.84e-11 * 1/(pa/pr) * sqrt(t/To)) + y) 
  # atmospheric attenuation coefficient in dB/m  
  coef.Aatm.dB = 20*log10(exp(1)) * coef.Aatm

  if (sum(is.na(r) == 0) >0 & sum(is.na(r0) == 0) >0)
  {
    # atmospheric attenuation in dB
    Aatm.dB = coef.Aatm.dB %*% t(r-r0)
    # atmospheric attenuation factor
    Aatm = exp(-coef.Aatm%*% t(r-r0))
  }
  else
  {
    Aatm.dB = NA
    Aatm = NA
  }

  return (list(factor = Aatm, db = Aatm.dB, coef = coef.Aatm, coef.db = coef.Aatm.dB))
}

########################

#********  Habitat attenuation
propa.Ahab<- function(f, r=NA, r0=NA, a0=0.02)
{
  " 
  Get the habitat attenuation factor and in dB and coefficients in Neper/m and dB/m  
  
  INPUTS :
  f: frequency in kHz [SCALAR]
  r : propagation distances in m [SCALAR or VECTOR]
  r0 : reference distance in m [SCALAR]
  a0 : attenuation coefficient of the habitat in dB/kHz/m [SCALAR]
  
  OUTPUT
  return a list with 
  - factor : habitat attenuation factor depending on the frequency, the habitat attenuation coefficient a0 and the distance [MATRIX] 
  => Multiply this value with the effective reference acoustic pressure p0 measured at r0 to estimate the pressure after attenuation
  - db : habitat attenuation in dB depending on the frequency, the habitat attenuation coefficient a0 and the distance [MATRIX] 
  => subtract this value to the reference sound pressure level (dB SPL) for each frequency and distance to estimate the sound pressure level after attenuation
  - coef : habitat attenuation coefficient in Neper/m
  - coeff.db : habitat attenuation coefficient in dB/m
  "
  
  coef.Ahab.dB = a0*f 
  coef.Ahab = coef.Ahab.dB / 20*log10(exp(1))

  if (sum(is.na(r) == 0) >0 & sum(is.na(r0) == 0) >0)
  {
    # habitat attenuation in dB
    Ahab.dB = coef.Ahab.dB %*% t(r-r0)
    # habitat attenuation factor
    Ahab = exp(-coef.Ahab%*% t(r-r0))
  }
  else
  {
    Ahab.dB = NA
    Ahab = NA
  }

  return (list(factor = Ahab, db = Ahab.dB, coef = coef.Ahab, coef.db = coef.Ahab.dB))
}

########################

#************ full attenuation
propa.Atotal  <- function(f, r, r0, t=20, rh=60, pa=101325, a0=0.02)
{
  " 
  get full attenuation factor and dB taking into account the geometric, atmospheric and habitat attenuation contributions
  
  INPUTS:
  f: frequency in kHz [SCALAR]
  r : propagation distances in m [SCALAR or VECTOR]
  r0 : reference distance in m [SCALAR]
  t: temperature in °C [SCALAR]
  rh: relative humidity in % [SCALAR]
  pa: atmospheric pressure in Pa [SCALAR]
  a0 : attenuation coefficient of the habitat in dB/kHz/m [SCALAR]
  
  OUTPUT:
  factor : return the total attenuation factor of an acoustic pressure [Pa]
  => Multiply this value with the effective reference acoustic pressure p0 measured at r0 to estimate the pressure after attenuation
  db : return the total attenuation in dB of an acoustic pressure [Pa]
  => subtract this value to the reference sound pressure level (dB SPL) for each frequency and distance to estimate the sound pressure level after attenuation
  "
  
  # Replicate the geometric attenuation vector (length of r) the number of time the length of f
  Ageo.matrix = matrix(rep(propa.Ageo(r,r0)$factor, each=length(f)), nrow=length(f))
  Atotal = Ageo.matrix * propa.Aatm(f,r,r0,t,rh,pa)$factor * propa.Ahab(f,r,r0,a0)$factor

  # Replicate the geometric attenuation vector (length of r) the number of time the length of f
  Ageo.matrix.dB = matrix(rep(propa.Ageo(r,r0)$db, each=length(f)), nrow=length(f))  
  Atotal.dB = Ageo.matrix.dB + propa.Aatm(f,r,r0,t,rh,pa)$db + propa.Ahab(f,r,r0,a0)$db  
  
  return (list(factor = Atotal, db = Atotal.dB))
}

########################

#************** Repartition of energy along the frequencies
propa.dBSPL_per_bin <- function(L, f)
{
  " 
  Function to spread the sound pressure level (Energy in dB) along a frequency vector (bins)   
  
  INPUTS:
  L : Sound Pressure Level in dB
  f: frequency in kHz [SCALAR or VECTOR]
  OUTPUT :
  Two vectors : 1st column is the frequency vector and the 2nd column is the sound pressure level corresponding to the frequency number of bins 
  "
  # init
  L_per_bin = rep(0, length(f))
  
  #  if f is a single value [SCALAR]
  if (length(f)==1)
  {
    nb_bin= 1
    # dB SPL for the frequency bandwidth
    L_per_bin = L - 10*log10(nb_bin) 
  }
  else # uf f is a vector [VECTOR]
  {  
    nb_bin= length(f)
    # dB SPL for the frequency bandwidth
    L_per_bin = L - 10*log10(nb_bin) 
  }
  
  vector = cbind(f,L_per_bin)
  
  return (list(f=vector[,1], db=vector[,2]))
}

########################

#************** detection distance
propa.detection_distance  <- function(L_bkg, L0, f, r0= 1, delta_r=1, t=20, rh=60, pa=101325, a0=0.02)
{
  " 
  get full attenuation factor taking into account the geometric, atmospheric and habitat attenuation contributions
  
  INPUTS:
  L_bkg : sound pressure level of the background in dB SPL [SCALAR, VECTOR]
  L0 : sound pressure level of the sound that is propagated
  f : frequency vector in kHz [VECTOR]
  r0 : distance at which L0 was measured (generally @1m)
  delta_r : distance resolution in m [SCALAR]
  t: temperature in °C [SCALAR]
  rh: relative humidity in % [SCALAR]
  pa: atmospheric pressure in Pa [SCALAR]
  a0 : attenuation coefficient of the habitat in dB/kHz/m [SCALAR]
  
  OUTPUT:
  distance_max : maximum distance of propagation before the sound pressure level is below the background [SCALAR or VECTOR]
  "
  # set f0 and f1 depending if freq is a single value or a vector
  if (length(f)>1) {f0 = f[1]; f1=f[length(f)]}
  else {f0 = f; f1 = NULL}
  
  # set the distance vector
  r = seq(1,10000,delta_r) 
  
  # set the distance max vector to store the result
  distance_max = rep(NA,length(f))

  # test for each frequency when the sound level L at distance r is below the background level L_bkg
  for (ii in (1:length(f)))
  {
    # Get the sound level after subtracting the full attenuation
    L = L0[ii] - propa.Atotal(f[ii], r, r0, t, rh, pa, a0)$db
    
    # distance max
    if (sum((L - L_bkg[ii])>0)>1) {
      distance_max[ii] = r[which.min((L - L_bkg[ii])[(L - L_bkg[ii])>0])] 
    } 
    else {
      distance_max[ii] = 0
    }
  }
  
  # return the frequency vector associated with the distance max
  return (cbind(f,distance_max))
}

########################

# #************* apply attenuation
propa.apply.att <- function(p0, fs, r, r0= 1, t=20, rh=60, pa=101325, a0=0.02)
{
  " 
  Apply attenuation of a temporal signal p0 after propagation between the reference distance r0 and the final distance r 
  taken into account the geometric, atmospheric and habitat attenuation contributions
  
  INPUTS:
  p0 : temporal signal (time domain) [VECTOR]
  fs: sampling frequency Hz [SCALAR]
  r : propagation distances in m [SCALAR or VECTOR]
  r0 : reference distance in m [SCALAR]
  t: temperature in °C [SCALAR]
  rh: relative humidity in % [SCALAR]
  pa: atmospheric pressure in Pa [SCALAR]
  a0 : attenuation coefficient of the habitat in dB/kHz/m [SCALAR]
  
  OUTPUT:
  p : temporal signal (time domain) after attenuation [VECTOR]
  "
  
  # Fourier domain
  P0 = fft(p0)/length(p0)
  f = seq(0,length(P0)-1) / length(P0) * fs /2
  # apply attenuation
  P = P0  * propa.Atotal(f/1000, r, r0, t, rh, pa, a0)$factor
  # Go back to the time domain
  p = fft(P, inverse=TRUE)
  # keep the real part
  p= Re(p)
  
  return (p)
}

########################

#************************** Get the pressure at r0 from simulation
propa.p_at_r0 <- function(f, r, p, r0=1, t=20, rh=60, pa=101325, a0=0.02)
{
  " 
    get the pressure at distance r0 from experimental values of pressure p measured at distance r.
    This function takes into account the geometric, atmospheric and habitat attenuations

  INPUTS:
    f : frequency vector in kHz [VECTOR]
    r : distance vector in m [VECTOR]
    p : pressure vector in Pa [VECTOR]
    r0 : distance where the pressure will be evaluated [SCALER]
    t: temperature in °C [SCALAR]
    rh: relative humidity in % [SCALAR]
    pa: atmospheric pressure in Pa [SCALAR]
    a0 : attenuation coefficient of the habitat in dB/kHz/m [SCALAR]

  OUTPUT:
    p0 : estimated pressure at distance r0 [SCALAR]
  "
  # Replicate the geometric attenuation vector (length of r) the number of time the length of f
  Ageo.matrix = matrix(rep(propa.Ageo(r,r0)$factor,each=length(f)),nrow=length(f))
  
  p0 = p * Ageo.matrix**(-1)
  # p0 = p0 * exp(func.coeff.Aa(f, t, rh, pa) %*% t(r-r0))
  p0 = p0 * (propa.Aatm(f, t, rh, pa)$factor)**(-1)
  p0 = p0 * (propa.Ahab(f, r, r0, a0)$factor)
  return (p0)  
} 