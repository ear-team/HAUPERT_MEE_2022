#=======================================================================================================#
#
#                                           LOAD LIBRARIES
#
#=======================================================================================================#

if (!require("signal")) install.packages("signal")
if (!require("fftw")) install.packages("fftw")
if (!require("rgl")) install.packages("rgl")
if (!require("rpanel")) install.packages("rpanel")
if (!require("seewave")) install.packages("seewave", repos="http://cran.at.r-project.org/")
if (!require("tuneR")) install.packages("tuneR")

library(seewave)
library(tuneR)
library(signal)

#=======================================================================================================#
#
#                                           My functions 
#
#=======================================================================================================#

##### Spectrum into bins
specBins <- function (spec, nbins = 256) 
{
  # Transform amplitude to energy in order to be able to change the number of points of the spectrum without changing the total energy
  PSD = spec^2
  N = floor(length(PSD)/nbins)
  res <- numeric(nbins)
  for (i in 1:nbins) {
    res[i] <- sum(PSD[((i-1)*N+1):(i*N)])
  }
  # SQRT : back to amplitude spectrum
  res = sqrt(res)
  return(res)
}

##### Tuckey window
win.tukey <- function(n, a = 0.5)
{
  M <- length(n) - 1
  n <- 0:(length(n) - 1)
  w <- rep(0, M)
  c <- n < a * (M) / 2
  w[c] <- 0.5 * (1 + cos(pi * (2 * n[c] / (a * (M)) - 1)))
  c <- a * (M) / 2 <= n & n <= (M) * (1 - a / 2)
  w[c] <- 1
  c <- (M) * (1 - a / 2) <= n
  w[c] <- 0.5 * (1 + cos(pi * (2 * n[c] / (a * (M)) - 2 / a + 1)))
  return(w)
}

#=================================================
##### Get noise signal (background noise or white noise)
extractNoise <- function(filename, channel ='left', am_timeline, total_dur, t0, dur)
{
  # Initialization of the matrix. Need to know the length which depends on the duration
  print(1)
  wavSM4 = readWave(filename, 
                    from = am_timeline[1],
                    to   = am_timeline[1]+total_dur,
                    units = 'seconds')
  if (channel=='left') {wav_temp= c(extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@left)}
  if (channel == 'right' )
  { 
    if (wavSM4@stereo == TRUE) {wav_temp= c(extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@right)}
    else {stop("This is not a stereo wavefile, so there is no right channel")}
  }
  wav_vector = matrix(NA, nrow=length(wav_temp), ncol=length(am_timeline))
  wav_vector[,1] = wav_temp

  for (ii in 2:length(am_timeline))
  {
    print(ii)
    wavSM4 = readWave(filename, 
                      from = am_timeline[ii],
                      to   = am_timeline[ii]+total_dur,
                      units = 'seconds')

    if (channel=='left') {wav_vector[,ii] = c(extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@left)}
    if (channel == 'right' )
    { 
      if (wavSM4@stereo == TRUE) {wav_vector[,ii] = c(extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@right)}
      else {stop("This is not a stereo wavefile, so there is no right channel")}
    }
  }
  return(wav_vector)
}

#=================================================
##### Get noise signal (background noise or white noise)
extractTemplate <- function(filename, channel ='left', am_timeline, total_dur = Inf, t0, dur, template, delay=0, Nrep=1, display=FALSE)
{
  # time beefore and after the sweep that we are looking for
  preT = 0.2
  postT = 0.2
  
  lag = c(NA)
  corrval = c(NA)
  dist = c(NA)
  
  # extract the portion of the template
  template = extractWave(template, from = 0, to = (preT+postT), xunit = "time")
  
  # Initialization of the matrix. Need to know the length which depends on the duration
  wavSM4 = readWave(filename, 
                    from = am_timeline[1],
                    to   = am_timeline[1]+total_dur,
                    units = 'seconds')
  fs_rec = wavSM4@samp.rate
  
  if (channel=='left') {wav_temp= extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@left}
  if (channel == 'right' )
  { 
    if (wavSM4@stereo == TRUE) {wav_temp= extractWave(wavSM4, from=t0, to=(t0+dur), xunit = "time")@right}
    else {stop("This is not a stereo wavefile, so there is no right channel")}
  }

  wav_vector = matrix(NA, nrow=length(wav_temp), ncol=length(am_timeline)*Nrep)  

  # loop along the timeline
  for (ii in 1:length(am_timeline))
  {
    # load the portion of the wavefile corresponding to the current timeline
    print(ii)
    wavSM4 = readWave(filename, 
                      from = am_timeline[ii],
                      to   = am_timeline[ii]+total_dur,
                      units = 'seconds')
    
    # loop depending on the repetition of the template
    for (kk in 1:Nrep)
    {
      t1 = t0 - preT + (kk-1)*(dur+delay)
      t2 = t0 + postT +(kk-1)*(dur+delay)
      
      # Extract the portion of the signal to analyse
      wav_tmp = extractWave(wavSM4, from=t1, to=t2, xunit = "time")

      if (display==TRUE) {plot(wav_tmp)}
      
      # Cross-correlation : template matching.
      if (channel=='left') {corr = ccf(wav_tmp@left, template@left, lag.max=length(wav_tmp@left), type="correlation", plot=display)}
      if (channel == 'right' )
      { 
        if (wavSM4@stereo == TRUE) {corr = ccf(wav_tmp@right, template@left, lag.max=length(wav_tmp@right), type="correlation", plot=display)}
        else {stop("This is not a stereo wavefile, so there is no right channel")}
      }

      # index
      idx = kk+(ii-1)*Nrep
      # lag
      lag[idx] = corr$lag[which.max(abs(corr$acf))] / fs_rec
      # correlation value. If corrval < 0.01, there is no correlation.
      corrval[idx] = max(corr$acf)

      # extraction of the signal with the corresponding lag to align all the signals together
      if (channel=='left') {wav_vector[,idx] = c(extractWave(wavSM4, from = (t1+lag[idx]), to = (t1+lag[idx]+dur), xunit = "time")@left)}
      if (channel == 'right' )
      { 
        if (wavSM4@stereo == TRUE) {wav_vector[,idx] = c(extractWave(wavSM4, from = (t1+lag[idx]), to = (t1+lag[idx]+dur), xunit = "time")@right)}
        else {stop("This is not a stereo wavefile, so there is no right channel")}
      }
      if (display==TRUE)
      {
        plot(extractWave(wavSM4, from = (t1+lag[idx]), to = (t1+lag[idx]+dur), xunit = "time"))
        spec(wav_vector[,idx], f=fs_rec, wn="rectangle", correction="amplitude", norm=FALSE, plot=TRUE)        
      }
    }
  }
  
  output <- list("wav" = wav_vector,"corrval"= corrval, "lag"=lag)
  return (output)
}


#=================================================
meanPSDfromNoise <- function(wav_vector, fs, NFFT, am_timeline)
{
  # Get the frequency vector
  freq = meanspec(wav_vector[,1], f=fs, wl=NFFT, wn="hamming", correction="amplitude", norm=FALSE, plot=FALSE)[,1]
  # Get the number of repetitions to average
  Nrep = dim(wav_vector)[2] / length(am_timeline) 
  # Get the length of the wave
  ll = dim(wav_vector)[1]
  
  # spectrum
  spec= matrix(NA,nrow=NFFT/2,ncol=length(am_timeline)*Nrep)
  for (ii in 1:dim(wav_vector)[2])
  {
    spec[,ii]= meanspec(wav_vector[,ii], f=fs, wl=NFFT, wn="hamming", correction="amplitude", norm=FALSE, plot=FALSE)[,2] 
  }
  
  # mean spectrum if number of repetition is >1
  if (Nrep>1)
  {
    spec.mean= matrix(NA,nrow=NFFT/2,ncol=length(am_timeline))
    for (ii in 1: length(am_timeline))
    {
      spec.mean[,ii] = rowMeans(spec[,((ii-1)*Nrep+1):(ii*Nrep)] )
    }
  }
  else
  {
    spec.mean = spec
  }
  
  #### mean spectrogram to Power Spectra Density (PSD)
  # Mean spectrogram is already corrected by the window used for the SDFT and corrected for taking into acount 1/2 of the spectrum
  PSD.mean = spec.mean^2/2

  output <- list("freq" = freq, "PSD"= PSD.mean)
  return (output)
}


#=================================================
meanPSDfromTemplate <- function(wav_vector, fs, NFFT, am_timeline)
{
  # Get the frequency vector
  freq = meanspec(wav_vector[,1], f=fs, wl=NFFT, wn="rectangle", correction="amplitude", norm=FALSE, plot=FALSE)[,1]
  # Get the number of repetitions to average
  Nrep = dim(wav_vector)[2] / length(am_timeline) 
  # Get the length of the wave
  ll = dim(wav_vector)[1]
  
  # spectrum
  spec= matrix(NA,nrow=NFFT/2,ncol=length(am_timeline)*Nrep)
  for (ii in 1:dim(wav_vector)[2])
  {
    tmp_spec = spec(wav_vector[,ii], f=fs, wn="rectangle", correction="amplitude", norm=FALSE, plot=FALSE)[,2]/ll
    # reduce the number of points in frequency in order to match the results obtain with meanspec 
    spec[,ii] = specBins(tmp_spec, nbins=(NFFT/2))
  }
  
  # mean spectrum if number of repetition is >1
  if (Nrep>1)
  {
    spec.mean= matrix(NA,nrow=NFFT/2,ncol=length(am_timeline))
    for (ii in 1: length(am_timeline))
    {
      spec.mean[,ii] = rowMeans(spec[,((ii-1)*Nrep+1):(ii*Nrep)] )
    }
  }
  else
  {
    spec.mean = spec
  }

  #### mean spectrogram to Power Spectra Density (PSD)
  # Mean spectrogram is already corrected by the window used for the SDFT and corrected for taking into acount 1/2 of the spectrum
  PSD.mean = spec.mean^2/2
  
  output <- list("freq" = freq, "PSD"= PSD.mean)
  return (output)
}

