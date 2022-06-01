library(seewave)
library(tuneR)


##### Change working directory to the current script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

wl = 512
 
##################################################
## CALIBRATION with soundlevel meter Svantek
##################################################
sin <- readWave("../../data/jbl_calibration/white_noise.wav")
sout <- readWave("../../data/jbl_calibration/jbl_frequency_response.wav")

## remove possible offset
sin <- rmoffset(sin, output="Wave")
sout <- rmoffset(sout, output="Wave")

## spectra
sin.spec <- meanspec(sin, wl=wl, plot=FALSE)
sout.spec <- meanspec(sout, wl=wl, plot=FALSE)

freq = sin.spec[,1]

## correction
sin.cor <- fir(sin, custom=sout.spec, bandpass=FALSE, output="Wave")
sin.cor.spec <- meanspec(sin.cor, wl=wl, plot=FALSE)

## save
savewav(sin.cor, file="../../data/jbl_calibration/noise_corrected.wav")
save(sin.cor.spec, file="../../data/jbl_calibration/JBL_frequency_response.Rdata")

# find the amplitude of the signal @1kHz in order to transform into dB
sin.spec_dB_ref = sin.spec[findInterval(1, sin.spec[,1])+1,2]
sout.spec_dB_ref = sout.spec[findInterval(1, sout.spec[,1])+1,2]
sin.cor.spec_dB_ref = sin.cor.spec[findInterval(1, sin.cor.spec[,1])+1,2]

# signal into dB
sin.spec_dB = 20*log10(sin.spec[,2]/sin.spec_dB_ref)
sout.spec_dB = 20*log10(sout.spec[,2]/sout.spec_dB_ref)
sin.cor.spec_dB = 20*log10(sin.cor.spec[,2]/sin.cor.spec_dB_ref)

# recreate a matrix, 1st col = Freq, 2nd col = dB
sin.spec_dB = cbind(freq, sin.spec_dB)
sout.spec_dB = cbind(freq, sout.spec_dB)
sin.cor.spec_dB = cbind(freq, sin.cor.spec_dB)

# mean spectrum (windowing : flattop in order to avoid ripples and average)
plot(x=NULL,xlab="Frequency (kHz)", xaxs="i", xlim=c(0,20),ylim=c(-40,6),
     ylab="Amplitude dB")
grid(col = "lightgray", lty = "dotted")
lines(sin.spec_dB, col=1, type="l", lty=1,lwd=2)
lines(sout.spec_dB, col=2, type="l", lty=1,lwd=2)
lines(sin.cor.spec_dB, col=4, type="l", lty=1,lwd=2)

rect(0,-6,max(freq),6, border = NA,col = rgb(0.5,0.5,0.5,1/4))
text(18, -1.2, '+/- 6dB')

legend("bottomleft", legend=c("input : white noise", "output : JBL frequency response before correction", "output : JBL frequency response after correction"), col=c(1,2,4), lty=1, bty="n")
