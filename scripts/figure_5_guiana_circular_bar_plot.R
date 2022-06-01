# clear workspace
rm(list = ls())

#=======================================================================================================#
#
#                                           LOAD LIBRARIES
#
#=======================================================================================================#

# https://www.r-graph-gallery.com/299-circular-stacked-barplot.html
# https://www.r-graph-gallery.com/295-basic-circular-barplot.html

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("grid")) install.packages("grid")
library(plotly)  
library(ggplot2)
library(gridExtra)
library(grid)

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


plot_list = list()

for (mm in seq(1:12))
{
  filename=paste(file_dir,root,month[mm],'.csv',sep='')
  
  # load data : filename, habitat, recorder, A0, p.bkg, p.exp, L0.exp, TEMP, RH, PS0
  df= read.csv(filename)
  
  # light
  hour = seq(0, 23)
  light= rep('NA',24)
  hour.day = c(seq(day[mm], night[mm]-1))
  hour.night = c(seq(0,day[mm]-1), seq(night[mm], 23))
  light[hour %in% hour.day] = 'day'
  light[hour %in% hour.night] = 'night'

  df$ps0 = ps0
  df$a0 = A0
  
  # select frequency
  df = df[c('hour', 'temperature', 'HR', 'ps0', 'a0', select_bandwidth[,1])]

  # set the L0 per frequency and hour of day
  df[select_bandwidth[,2]] = rep(L0,24)

  # adjust L0 depending on the SM4 frequency response
  if (CORRECTION_RECORDER == TRUE)
  {
    #====================================================================================================#
    #              To adapt the theoritical sound pressure value to the SM4 frequency response (which is not flat)
    #====================================================================================================#
    # load the GENERIC frequency response of the SM4. 
    # The gain to correct the frequency response of the SM4 was measured on a single SM4 recorder, 
    # assuming that all SM4 have the exact same frequency response
    # It would be better to adapt the correction for each SM4 if possible.
    # load('SM4_corr_frequency_response_dB_N2048.Rdata')
    load('./toolbox/SM4_gain.Rdata')
    FREQ_GAIN_CORR = SM4.G[,1]
    GAIN_CORR = SM4.G[,2] 
    # linear interpolation of SM4.G in order to match the frequency bin
    SM4.GAIN = power2dB(approx(FREQ_GAIN_CORR, dB2power(GAIN_CORR), seq(0.5,length(df[select_bandwidth[,2]]),1), rule=2)$y)
    # Set L0 as it is "seen" by the SM4 (higher sensitivity around (1-8kHz), low sensitivity for >10kHz)
    df[select_bandwidth[,2]]  = t(t(df[select_bandwidth[,2]]) - SM4.GAIN)
  } 
  
  # get the maximum listening distance 
  for (row in (1:nrow(df)))
  {
    dmax = propa.detection_distance(L_bkg=unlist(df[row,select_bandwidth[,1]]), L0=unlist(df[row,select_bandwidth[,2]]), f=fn, r0= 1, delta_r=1, t=unlist(df[row,'temperature']), rh=unlist(df[row,'HR']), pa=unlist(df[row,'ps0']), a0=unlist(df[row,'a0']))

    df$dmax[row] = list(as.vector(dmax[,2]))
  }
  
  
  # set the dmax per frequency and hour of day
  df[select_bandwidth[,3]] = NA
  for (row in (1:nrow(df)))
  {
    for (ff in (1:nrow(select_bandwidth)))
    {
      df[row,select_bandwidth[ff,3]] = df$dmax[[row]][ff]
    }
  }
  
  df = df[c('hour', select_bandwidth[,3])]
  
  # Create dataset
  data =data.frame(
    group = as.factor(df$hour),
    name = paste(df$hour, "h ", sep=""),
    light = light,
    df[select_bandwidth[,3]] )# keep only dmax per freq and hour
  data =data %>% gather(key="observation", value="value", -c(1,2,3))

  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=1
  to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
  colnames(to_add) = colnames(data)
  to_add$group=rep(levels(data$group), each=empty_bar)
  data=rbind(data, to_add)
  data=data %>% arrange(group)
  data$id=seq(1, nrow(data))
  
  # Get the name and the y position of each label
  label_data=data
  number_of_bar=nrow(label_data)
  angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust<-ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data=data %>%
    group_by(group) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  # drop 2 rows every 3 rows
  keeps = seq(1,23,3)

  # prepare a data frame for base lines
  df1kHz = subset(data, data$observation == 'X1kHz', select = c('group','id','value'))
  df1kHz['start'] = df1kHz$id
  df1kHz['end'] = c(df1kHz$start[-1], df1kHz$start[1])
  df1kHz['value_end'] = c(df1kHz$value[-1], df1kHz$value[1])
  df1kHz = head(df1kHz,-1) # remove last row
  
  df2kHz = subset(data, data$observation == 'X2kHz', select = c('group','id','value'))
  df2kHz['start'] = df2kHz$id
  df2kHz['end'] = c(df2kHz$start[-1], df2kHz$start[1])
  df2kHz['value_end'] = c(df2kHz$value[-1], df2kHz$value[1])
  df2kHz = head(df2kHz,-1) # remove last row
  
  df5kHz = subset(data, data$observation == 'X5kHz', select = c('group','id','value'))
  df5kHz['start'] = df5kHz$id
  df5kHz['end'] = c(df5kHz$start[-1], df5kHz$start[1])
  df5kHz['value_end'] = c(df5kHz$value[-1], df5kHz$value[1])
  df5kHz = head(df5kHz,-1) # remove last row
  
  df8kHz = subset(data, data$observation == 'X8kHz', select = c('group','id','value'))
  df8kHz['start'] = df8kHz$id
  df8kHz['end'] = c(df8kHz$start[-1], df8kHz$start[1])
  df8kHz['value_end'] = c(df8kHz$value[-1], df8kHz$value[1])
  df8kHz = head(df8kHz,-1) # remove last row
  
  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$start = grid_data$start -10
  grid_data$end = grid_data$start +11
  grid_data=grid_data[-1,]

  # Make the plot
  p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
    # For day and night annular
    geom_bar(aes(x=as.factor(id), y=rep(20,length(id)), fill=light), stat="identity", alpha=0.75, width = 2) +
    
    # Add a val=200/150/100/50 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=0.75, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey", alpha=0.75, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey", alpha=0.75, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey", alpha=0.75, size=0.3 , inherit.aes = FALSE ) +
  
    # Add text showing the value of each 200/150/100/50 lines
    ylim(-100,275) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +

    # Add base line information
    geom_text(data=base_data[keeps,], aes(x = title, y = -30, label=paste(group,'h',sep='')), hjust=c(rep(0,length(keeps)/2)+0.4, rep(1,length(keeps)/2)-0.4), colour = "black", alpha=0.75, size=2, fontface="bold", inherit.aes = FALSE) +
  
    # Add line as spidergraph to connect same frequencies
    geom_segment(data=df1kHz, aes(x =start, y = value, xend = end, yend = value_end), colour = "darkblue", alpha=0.75, size=0.75, inherit.aes = FALSE ) +
    geom_segment(data=df2kHz, aes(x =start, y = value, xend = end, yend = value_end), colour = "darkgreen", alpha=0.75, size=0.75, inherit.aes = FALSE ) +
    geom_segment(data=df5kHz, aes(x =start, y = value, xend = end, yend = value_end), colour = "darkorange", alpha=0.75, size=0.75, inherit.aes = FALSE ) +
    geom_segment(data=df8kHz, aes(x =start, y = value, xend = end, yend = value_end), colour = "darkred", alpha=0.75, size=0.75, inherit.aes = FALSE ) +
    
    # Add vertical scale legend
    annotate("text", x = rep(max(data$id),4), y = c(50, 100, 150, 200), label = c("50m", "100m", "150m", "200m") , color="black", size=2.5 , angle=0, fontface="bold", hjust=1) +
    # Add month legend
    annotate("text", x = c(0), y = c(200), label = c(month[mm]) , color="black", size=6 , angle=0, fontface="bold", hjust=-1) +
    
    scale_fill_grey(start=1, end=0)
  plot_list[[mm]] = p
}


library(gridExtra)
library(grid)

legd <- legendGrob(c("0-1kHz", "1-2kHz", "4-5kHz","7-8kHz"), nrow=1, do.lines = TRUE,
                   gp=gpar(col = c("darkblue","darkgreen","darkorange","darkred"), cex=0.9))
# draw
grid.arrange(grobs=plot_list,ncol=3, bottom=legd)

# doesn't draw but create an object that can be saved
mp = arrangeGrob(grobs=plot_list,ncol=3, bottom=legd)
ggsave(file=paste(SITE,'_circular_plot.pdf',sep=''), mp, width = 16, height = 24, units = "cm")
