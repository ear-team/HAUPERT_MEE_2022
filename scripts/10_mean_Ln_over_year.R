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

CORRECTION_RECORDER = TRUE # Adjust the theoritical L0 to the SM4 in order to simulate the value as it is "seen" by the SM4
SITE = 'GUIANA' # JURA or GUIANA

# simulation with a white noise L0 = 80 dB
L0 = 80 # Initial Sound pressure level
L0_per_bin = propa.dBSPL_per_bin(L=L0, f=seq(1,20,1))$db[1] # repartition of the initial sound level over the whole frequency band (0-20kHz)
L0 = L0_per_bin

# all data
# select_bandwidth <- cbind(c('LEQf_0to1kHz','LEQf_1to2kHz', 'LEQf_2to3kHz','LEQf_3to4kHz','LEQf_4to5kHz','LEQf_5to6kHz','LEQf_6to7kHz','LEQf_7to8kHz','LEQf_8to9kHz','LEQf_9to10kHz'),
#             c('L0_0to1kHz','L0_1to2kHz', 'L0_2to3kHz','L0_3to4kHz','L0_4to5kHz','L0_5to6kHz','L0_6to7kHz','L0_7to8kHz','L0_8to9kHz','L0_9to10kHz'),
#             c('1kHz','2kHz','3kHz','4kHz','5kHz','6kHz','7kHz','8kHz','9kHz','10kHz'))
# remove some frequency
select_bandwidth <- cbind(c('LEQf_0to1kHz','LEQf_1to2kHz', 'LEQf_2to3kHz','LEQf_3to4kHz','LEQf_4to5kHz','LEQf_5to6kHz','LEQf_6to7kHz','LEQf_7to8kHz'),
                          c('L0_0to1kHz','L0_1to2kHz', 'L0_2to3kHz','L0_3to4kHz','L0_4to5kHz','L0_5to6kHz','L0_6to7kHz','L0_7to8kHz'),
                          c('1kHz','2kHz','3kHz','4kHz','5kHz','6kHz','7kHz','8kHz'))

fn = seq(1,dim(select_bandwidth)[1]) # Frequency vector in kHz

#=======================================================================================================#
#
#                                           LOAD DATA
#
#=======================================================================================================#

## ===================================================================================================================
#              loop to construct the listening distance of the 2 devices in 2 env when a sound at 80dB is produced
# ===================================================================================================================

# Month
month = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Guyane
if (SITE =='GUIANA')
{
  root = 'guiana_SPL_ground_perHour_'
  file_dir = '../data/guiana/YEAR_AMBIENT_SOUND/'
  day = rep(6,12)
  night = rep(19,12)
  # get the constants
  A0 = 0.011
  ps0 = 101325
}

# Jura
if (SITE =='JURA')
{
  root = 'jura_SPL_perHour_'
  file_dir = '../data/jura/YEAR_AMBIENT_SOUND/'
  day = c(9,8,7,6,6,6,6,7,7,8,8,9)
  night = c(18,18,19,20,21,21,21,20,19,19,18,18)
  # get the constants
  A0 = 0.024
  ps0 = 88000
}


mean_Ln = data.frame()

for (mm in seq(1:12))
{
  filename=paste(file_dir,root,month[mm],'.csv',sep='')
  df= read.csv(filename)
  df_Ln = df[select_bandwidth[,1]]
  if (length(mean_Ln) == 0) 
    {mean_Ln = 10^(df_Ln/20)}
  else 
    {mean_Ln = mean_Ln + 10^(df_Ln/20)}
}
mean_Ln = 20*log10(colMeans(mean_Ln)/12)
  
  