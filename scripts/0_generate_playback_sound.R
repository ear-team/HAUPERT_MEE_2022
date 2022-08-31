library(seewave)
library(tuneR)


##### Change working directory to the current script directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

wl=512

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

##### Sampling Frequency (Recording, Generation)
fs_rec = 44100
NN = 8
fs_gene = fs_rec*NN

###### load frequency response
load('/home/haupert/DATA/mes_articles/2022/MEE_PROPA/HAUPERT_2022.git/data/jbl_calibration/JBL_frequency_response.Rdata')
freq_response = sin.cor.spec
tf = freq_response[,2]
#plot(tf)
tf = c(tf, rep(1,(NN-1)*length(tf)))
#plot(tf, type='l')

## frequency resolution
(freq <- seq(500, 22000, by=500))
freq <- c(freq, 1000)
length(freq)

############# amorce
# with JBL correction
am <- synth(f=fs_gene, d=6, cf=500, harmonics=rep(1,43), output="Wave")
am_comp <- fir(am, custom=tf,wl =length(tf)*2, bandpass = TRUE, output="Wave")
am_comp <- normalize(am_comp)
am_comp <- win.tukey(am_comp, a=0.01)*am_comp

########### bruit blanc
# with JBL correction
wn <- noisew(f=fs_gene, d=6, output="Wave")
wn_comp <- fir(wn, custom=tf,wl =length(tf)*2, bandpass = TRUE, output="Wave")
wn_comp <- normalize(wn_comp) # trick to convert into 32bits
wn_comp <- win.tukey(wn_comp, a=0.01)*wn_comp

# plot each spectrum
## visualise spectra
## Full spectrum (no windowing and no average)
plot(spec(wn_comp, wn="rectangle", plot=FALSE),
     xlab="Frequency (kHz)", xaxs="i", xlim=c(0,22),
     ylab="Amplitude", ylim=c(0,1.2), yaxt="n",
     type="l")
legend("topright", legend=c("white noise"), col=1:3, lty=1, bty="n")


####################################
#
#       EQUALIZER
#
######################################
dB_am <- c(67.6, 67.4, 67.4, 67.4, 67.4)
am_wav <- convSPL(dB_am, d=0.5)$p
am_wav <- mean(am_wav)

dB_wn <- c(65.4, 65.3, 65.4, 65.5, 65.4)
wn_wav <- convSPL(dB_wn, d=0.5)$p
wn_wav <- mean(wn_wav)

### calibration
wav <- c(am_wav,
         wn_wav)
calibration <- s1kHz_wav / wav
calibration
#### calibration=> (7.740514 2.306320 1.000000 4.934258 2.306320)

## equalizer Amorce : should be set at the max (max(calibration), instead of calibration[1])
am_comp_calib <- calibration[1]*am_comp
## equalizer white noise
wn_comp_calib = calibration[2]*wn_comp

## silence in between
silbetween <- silence(duration=0.5, samp.rate=fs_gene, xunit="time", bit=32, pcm=TRUE) 

## assemblage des sons
assemblage <-bind(am_comp_calib, silbetween,
                  wn_comp_calib, silbetween,
)

plot(assemblage/max(assemblage@left), col=4, xlab="Time (s)",  xaxs="i", ylab="Amplitude")
grid(col = "lightgray", lty = "dotted")

#savewav(assemblage, file="assemblage_JBL_comp_cal.wav")
