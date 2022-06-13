# Haupert S. & al. Physics-based model to predict the acoustic detection distance of terrestrial autonomous recorder units over the diel cycle and across seasons: insights from an Alpine and a Neotropical forest. Journal to be defined (XXXX)

This repository contains all the functions in 'R' that were used to process the data and create figures of the publication 
[Haupert et. al. XXXX, journal to be defined](https://www.to.be.defined.fr)

The code is available for reproducibility.

This repository is archived on Zenodo: [![DOI](https://zenodo.org/badge/DOI/XXXX/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.3530203)

If the code, even partially, is used for other purpose please cite the paper `Haupert S., Sèbe F., Sueur J. Physics-based model to predict the acoustic detection distance of terrestrial autonomous recorder units over the diel cycle and across seasons: insights from an Alpine and a Neotropical forest. Journal to be defined (XXXX)`

## Setup and usage

Download the `.zip` from Github (click on `code` then `Download Zip`) and extract all folders without changing the name of the folders neither rearrange the folder and sub-folders.

Then, use your favorite R environment (e.g. RStudio). The scripts are ready to be used. They will install the required libraries if they are not already installed in your environment. Here is the list of library that are requested :
* seewave and tuneR
* plotly
* robutslrm
* tidyverse, ggplot2, dplyr

The R scripts to reproduce figures/analyses presented in the paper are: 
* [`figure_1bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_1bc_ambient_sound_subtraction.R) and [`figure_2bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_2bc_ambient_sound_subtraction.R) : Ambient sound and white noise variations according to distance in the Neotropical rainforest (French Guiana) or the Alpine coniferous forest (Jura). (a) Time-frequency representation of the sounds recorded by the ARU of the ambient sound and the white noise between 10 m and 100 m. (b) Mean spectra of the total sound pressure level (dB SPL) (i.e. ambient sound + white noise) measured at each distance. (c) Map of sound levels: overall transmitted signal at each propagation distance (left), ambient sound (centre) and white noise (right). The unit is dB SPL (re20 µPa). The audio signals were collected by the ARU on 18 February 2019 between 2 pm and 5 pm.
* [`figure_3abcd_exp_vs_simulation_comparison.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_3abcd_exp_vs_simulation_comparison.R): Comparison between experimental data acquired with the ARU in a Neotropical rainforest (French Guiana) and theoretical propagation curves obtained from the full model as described by Eq. 4. The comparison was performed for the following bandwidths: (a) 0–5 kHz; (b) 5–10 kHz; (c) 10–15 kHz; and (d) 15–20 kHz. The habitat coefficient attenuation parameter was a0 = 0.011 dB/kHz/m. L (model) = predicted sound level according to the full propagation model; Lexp (exp) = experimental sound level; L - Ageo (model) = predicted sound level according to a partial propagation model taking into account only the geometric attenuation; L - (Ageo + Aatm) (model) = predicted sound level according to a partial propagation model taking into account both, the geometric and atmospheric attenuation; Ln (exp) = experimental ambient sound level: L (model) + Ln (exp) = experimental ambient sound level added to the predicted sound level according to the full propagation model. Similar results are provided for the Alpine coniferous forest and for the sound meter level in Appendix F.
* [`figure_4a_detection_distance_SM4_guiana_80dB.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_4a_detection_distance_SM4_guiana_80dB.R), [`figure_4b_detection_distance_SM4_jura_80dB.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_4b_detection_distance_SM4_jura_80dB.R): Detection distance of the ARU obtained with the model in: (a) a Neotropical rainforest (French Guiana); (b) an Alpine coniferous forest (Jura, France). The source was a wideband (0–20 kHz) white noise broadcast at 80 dB SPL re20 µPa at 1 m. The contribution of geometric (i.e. spreading loss) attenuation (Ageo, green), atmospheric absorption (Aatm, red) and habitat attenuation (Ahab, purple) are depicted with different colours.
* [`figure_5_guiana_circular_bar_plot.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_5_guiana_circular_bar_plot.R) and [`figure_6_jura_circular_bar_plot.R`](https://github.com/shaupert/haupert_2022/blob/master/scripts/figure_6_jura_circular_bar_plot.R) : Detection distance variation of the ARU according to night and day cycle and month in a Neotropical rainforest (French Guiana) or an Alpine coniferous forest (Jura) taken into account the variation of the environmental factors (T °C and RH %) as well as the  ambient sound level Ln. The sound source was a wideband (0–20 kHz) white noise broadcast at 80 dB SPL re20 µPa at 1 m. Each circular plot represents one month, with each circle corresponding to a propagation distance from 50 to 200 m by a step of 50 m over a complete day (white) and night (black) cycle. Line colours refer to four 1 kHz frequency bands (0–1 kHz, 1–2 kHz, 4–5 kHz, 7–8 kHz).

R scripts to obtain the figures in the supplementary are in this [`subdirectory`](https://github.com/shaupert/HAUPERT_2022/tree/master/scripts/supplementary) while the R scripts and the correction factors to correct the SM4's frequency response are in the subfolder [`toolbox`](https://github.com/shaupert/HAUPERT_2022/tree/master/scripts/toolbox) 

## Google Colab
<p align="center">
  <img src="https://s2.qwant.com/thumbr/474x190/f/9/aae347431a927c9b5deb63431ea29c0dd6fceb9210443fdd6bb9b3dba23146/th.jpg?u=https%3A%2F%2Ftse2.mm.bing.net%2Fth%3Fid%3DOIP.IVRAF7_KdEVWUFq1wmDvmQHaC-%26pid%3DApi&q=0&b=1&p=0&a=0g" />
</p>

Google Colab is a very convenient way to test a code written in Python or R without having to install anything on your computer. Everything is run iin the cloud, directly on computers in Google. The required libraries are installed on the fly on the distant computers.

For an interactive way to test our code, we provide most of the R scripts from this repository as Notebook ready to be used on Google Colab. Some adjustment had 
to be done to be able to run the scripts as some R libraries such as Plotly for R don't work properly on Google Colab. In such a case, we mixed R and Python. R 
is used for data processing while Python is used to build and display the figures.

All the notebooks are available [here](https://drive.google.com/drive/folders/1p_xJDaCP2ynVswfaWJLCYIJLMLw9NRic?usp=sharing). 

## step by step tutorial

1. The loudspeaker used for the field experiment should be calibrated in order to correct its frequency response and obtain a frequency curve as flat as possible. The full process is explained in Appendix A of the paper and the R script is [`here`](https://github.com/shaupert/HAUPERT_2022/blob/master/scripts/supplementary/figure_A1_frequency_response_jbl.R) 


## Integration of the propagation functions into libraries

For convenience, we added the propagation functions into the library [seewave](https://rug.mnhn.fr/seewave/)  (R language) and into the package [scikit-maad](https://scikit-maad.github.io/) (Python language)

<p align="center">
  <img src="https://rug.mnhn.fr/seewave/PICT/seewave_logo.png" />
</p>
<p align="center">
  <img src="https://scikit-maad.github.io/_images/maad_key_visual_black.png" />
</p>

