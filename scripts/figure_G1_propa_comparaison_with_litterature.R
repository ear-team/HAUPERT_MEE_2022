# clear workspace
rm(list = ls())

library(plotly)

##### CALL functions from MY TOOLBOX
source("/home/haupert/DATA/mes_projets/_TOOLBOX/R/propagation/toolbox_propa.R")


a0_jura = 0.022   # average value between SVANTEK and SM4
a0_guiana = 0.015 # average value between SVANTEK and SM4

f = seq(1000,10000,1000)
f2 = seq(200,4000,1000)
f3 = seq(2000,10000,1000)
f4 = seq(200,20000,1000)

# references
# [1] M.A.Price, K.Attenborough and N.W.Heap, Sound attenuation through trees: measurements and models, J. Acoust. Soc. Am., 84(5): 1836–1844 (1988)
# [2] Huisman W.H.T.Huisman, Sound propagation over vegetation-covered ground, PhD Thesis, University of Nijmegen, The Netherlands (1990)
# [3] ISO 9613-2
# [4] Ellinger 2003 @2.5m
# [5] Gronberg 2015
# [6] Watanabe JASJapan 1996

p <- plot_ly(x =f, y =0.7*log10(f)-2.03, type = 'scatter', mode = 'lines+markers', name='Mixed conifers [1]', line=list(color='green',dash = 'dot'), marker=list(color='green', symbol='diamond'), opacity=1)
p <- add_trace(p, x = f, y = 0.4*log10(f)-1.2, type = 'scatter', mode = 'lines+markers', name='Mixed oak and spruce [1]', line=list(color='red',dash = 'dot'), marker=list(color='red',symbol='diamond'))

# p <- plot_ly(x = f, y = 0.4*log10(f)-1.2, type = 'scatter', mode = 'lines+markers', name='mixed oak and spruce [1]', line=list(color='red',dash = 'dot'), marker=list(color='red',symbol='diamond'))
p <- add_trace(p, x = f, y = 0.26*log10(f)-0.75, type = 'scatter', mode = 'lines+markers',name='Spruce monoculture [1]', line=list(color='blue',dash = 'dot'), marker=list(color='blue',symbol='diamond-open'))
p <- add_trace(p, x = c(2000,3000,4000,5000,6000), y = 0.39*log10(c(2000,3000,4000,5000,6000))-1.23, type = 'scatter', mode = 'lines+markers',name='Pine forest [2]', line=list(color='orange',dash = 'dot'), marker=list(color='orange',symbol='square')) 

# Watanabe
# p <- add_trace(p, x = f2, y = -10*log10(1-1/8*0.0015*5.45*sqrt(f2)), type = 'scatter', mode = 'lines+markers', name='Laboratory max [6]', line=list(color='grey',dash = 'dash'), marker=list(color='grey',symbol='triangle-up'))
# p <- add_trace(p, x = f2, y = -10*log10(1-1/8*0.002*6.53*sqrt(f2)), type = 'scatter', mode = 'lines+markers', name='Forest in laboratory [6]', line=list(color='grey',dash = 'dash'), marker=list(color='grey',symbol='triangle-down'))

# /(92m-4.5m), frequency in kHz
k = 2.832 # 2.5m : 2.832 5m : 2.655
d = 2.758 # 2.5m: 2.758 5m : 0.534
p <- add_trace(p, x = f3, y = (k*f3/1000+d)/(92-4.5), type = 'scatter', mode = 'lines+markers', name='Neotropical rainforest [4]', line=list(color='green',dash = 'dash'), marker=list(color='green',symbol='diamond')) 
# suppression de l'atténuation atmosphérique
# p <- add_trace(p, x = f3, y = (k*f3/1000+d)/(92-4.5) - propa.get.Aatm.dB(f3/1000,r0=4.5,r=92, t=27, rh=80)/(92-4.5), type = 'scatter', mode = 'lines+markers', name='Neotropical Forest [5]', line=list(color='green',dash = 'dash'), marker=list(color='green',symbol='diamond')) 
# equation linéaire passant par l'origine (d=0)
# p <- add_trace(p, x = f3, y = (k*f3/1000)/(92-4.5) - propa.get.Aatm.dB(f3/1000,r0=4.5,r=92, t=27, rh=80)/(92-4.5), type = 'scatter', mode = 'lines+markers', name='Neotropical Forest [5]', line=list(color='green',dash = 'dash'), marker=list(color='green',symbol='diamond')) 

# Gronberg
p <- add_trace(p, x = f4, y = 0.038+8.199*1e-9*f4**1.808, type = 'scatter', mode = 'lines+markers', name='Mixed eucalyptus and acacias [5]', line=list(color='red',dash = 'dash'), marker=list(color='red',symbol='diamond-open')) 


# ISO
p <- add_trace(p, x = c(250,500, 1000,2000,4000,8000), y = c(0.04,0.05, 0.06,0.08,0.09,0.12), type = 'scatter', mode = 'lines+markers', name='ISO 9613-2 [3]', line=list(color='purple',dash = 'solid'), marker=list(color='purple',symbol='circle'))


# my study
p <- add_trace(p, x = f4, y = f4*a0_guiana/1000, type = 'scatter',  mode = 'lines+markers', name='this study (Neotropical rainforest)', line=list(color='black',dash = 'solid'),marker=list(color='black',symbol='cross'),opacity=1) 
p <- add_trace(p, x = f4, y = f4*a0_jura/1000, type = 'scatter',  mode = 'lines+markers', name='this study (Temperate cold forest)', line=list(color='black',dash = 'dot'),marker=list(color='black',symbol='open circle'),opacity=1) 

# Found with our data when avering data in order to test Price law EA = m * log10(f) - c
# p <- add_trace(p, x = f, y = 0.130*log10(f)-0.36, type = 'scatter', mode = 'lines+markers',name='this study (EA=mlog10(f)-c)', line=list(color='black',dash = 'solid'), marker=list(color='black',symbol='diamond-open'))

p <- layout(p,
            legend = list(orientation="h",
                          yanchor="bottom",
                          y=1.02,
                          # xanchor="left",
                          x=0, 
                          font = list(size=10)),
            xaxis = list(title = 'Frequency [Hz]',
                         gridcolor = 'rgb(255,255,255)',
                         showgrid = TRUE,
                         zeroline = FALSE,
                         type = "log"),
            yaxis = list(title = 'Excess attenuation [dB/m]',
                         gridcolor = 'rgb(255,255,255)',
                         zeroline = FALSE,
                         showgrid = TRUE),
            paper_bgcolor='rgb(255,255,255)', 
            plot_bgcolor='rgb(229,229,229)')
p
