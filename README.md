# Haupert S. & al. Ambient noise level drives the acoustic range of detection of the biophony over seasons and diel cycles Journal to be defined (XXXX)

This repository contains all the functions in 'R' that were used to process the data and create figures of the publication 
[Haupert et. al. XXXX, journal to be defined](https://www.to.be.defined.fr)

The code is available for reproducability.

This repository is archived on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/XXXX/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.3530203)

If the code, even partially, is used for other purpose please cite the paper `Haupert S., Sèbe F., Sueur J. Ambient noise level drives the acoustic range of detection of the biophony over seasons and diel cycles journal to be defined, XXXX`

## Setup and usage

Download the `.zip` from Github (click on `code` then `Download Zip`) and extract all folders without changing the name of the folders neither rearrange the folder
and sub-folders.

Then, use your favorite R environment (e.g. RStudio). The scripts are ready to be used. They will install the required libraries if they are not already installed
in your environment. Here is the list of library that are requested :
* seewave and tuneR
* plotly
* robutslrm
* tidyverse, ggplot2, dplyr

The R scripts to reproduce figures/analyses from our paper are:
 
* [`figure_1a_compute_Leq.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_1a_compute_Leq.R) and [`figure_1b_compute_Leq.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_1b_compute_Leq.R): Equivalent continuous sound level (Leq) of ambient sound and white noise, measured with the ARU in (a) neotropical rain forest (French Guiana, France) and (b) Alpine coniferous forest (Jura, France)). The initial source level of the white noise was respectively 83 dB and 78 dB at 1m
* [`figure_2bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_2bc_ambient_sound_subtraction.R), [`figure_3bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_3bc_ambient_sound_subtraction.R), [`figure_D1bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_D1bc_ambient_sound_subtraction.R) and [`figure_D2bc_ambient_sound_subtraction.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_D2bc_ambient_sound_subtraction.R) : Ambient sound and white noise variation according to distance in the neotropical rain forest (French Guiana). (a) Time-frequency representation of the sounds recorded by the ARU of the ambient sound and the white noise between 10 m to 100 m. (b) Mean spectra of the total sound pressure level (dB SPL) (i.e. ambient sound + white noise) measured at each distance. (c) Maps sound levels: overall transmitted signal at each propagation distance (left),  ambient sound (center) and white noise (right). The unit is dB SPL (re 20 µPa).
* [`figure_4a_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_4a_EA_per_m_vs_frequency.R), [`figure_4b_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_4b_EA_per_m_vs_frequency.R), [`figure_4c_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_4c_EA_per_m_vs_frequency.R), [`figure_4d_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_4d_EA_per_m_vs_frequency.R), [`figure_F2a_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_F2a_EA_per_m_vs_frequency.R), [`figure_F2b_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_F2b_EA_per_m_vs_frequency.R), [`figure_F2c_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_F2c_EA_per_m_vs_frequency.R) and [`figure_F2d_EA_per_m_vs_frequency.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_F2d_EA_per_m_vs_frequency.R) : Frequency dependence of the excess attenuation EA normalized by the propagation distance (EA/(r-r0)) (see Eq.12) for different reference distances r0 in the Alpine coniferous forest (Jura) for the ARU (a) with and (b) without subtraction of ambient sound and for the sound meter level (c) with and (d) without subtraction of ambient sound. Similar figures for the neotropical rain forest (French Guiana) are available in appendix F.
* [`figure_5abcd_exp_vs_simulation_comparison.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_5abcd_exp_vs_simulation_comparison.R), [`figure_babcd_exp_vs_simulation_comparison.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_babcd_exp_vs_simulation_comparison.R), [`figure_E1abcd_exp_vs_simulation_comparison.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_E1abcd_exp_vs_simulation_comparison.R), and [`figure_E2abcd_exp_vs_simulation_comparison.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts) : Comparison between experimental data acquired with the ARU or sound meter level, in a neotropical rainforest (French Guiana) and simulated data obtained from the complete model as described by Eq.12. The comparison was performed for the following bandwidths: (a) 0-5 kHz, (b) 5-10 kHz, (c) 10-15 kHz, and (d) 15-20 kHz. The habitat coefficient attenuation parameter was a0 = 0.011 dB/kHz/m. Similar results are provided for the Alpine coniferous forest in appendix E.
* [`figure_7a_detection_distance_SM4_guiana_80dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_7a_detection_distance_SM4_guiana_80dB.R), [`figure_7b_detection_distance_SM4_jura_80dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_7b_detection_distance_SM4_jura_80dB.R),[`figure_H1a_detection_distance_Svantek_guiana_80dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_H1a_detection_distance_Svantek_guiana_80dB.R) and [`figure_G1b_detection_distance_SM4_jura_bird100dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_G1b_detection_distance_SM4_jura_bird100dB.R) :  Detection distance of the ARU or the sound meter level obtained with the model in (a) a neotropical rainforest (French Guiana) and in (b) an Alpine coniferous forest (Jura, France). The source was a wideband (0-20kHz) white noise broadcasted at 80 dB SPL re20µPa @1m. The contribution of geometric (i.e. spreading loss) attenuation (Ageo, green), atmospheric absorption (Aatm, red) and habitat attenuation (Ahab, purple) are depicted with different colors.
* [`figure_G1a_detection_distance_SM4_jura_bird75dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_G1a_detection_distance_SM4_jura_bird75dB.R) and [`figure_G1b_detection_distance_SM4_jura_bird100dB.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_G1b_detection_distance_SM4_jura_bird100dB.R) : Detection distance of an SM4 ARU obtained by the model in an Alpine coniferous forest (Jura, France) for the quietest bird species (75 dB SPL re20µPa @1m) and the loudest one (100 dB SPL re20µPa @1m) with frequency range between 2 kHz to 8 kHz. The contribution of geometric (i.e. spreading loss) attenuation (Ageo, green), atmospheric absorption (Aatm, red) and habitat attenuation (Ahab, purple) are depicted with different colors. T°C, RH, Ps and a0 values are set to the experimental values found during the Jura transect
* [`figure_8_guiana_circular_bar_plot.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_8_guiana_circular_bar_plot.R) and [`figure_9_jura_circular_bar_plot.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_9_jura_circular_bar_plot.R) : Detection distance variation according to day/night cycle and month in the neotropical rainforest (French Guiana, France) and an alpine coniferous forest (Jura, France). The sound source was a wideband (0-20kHz) white noise broadcasted at 80 dB SPL re20µPa @1m. Each circular plot represents one month with each circle corresponding to a propagation distance from 50 to 250 m by a step of 50 m over a complete day (white) and night cycle (black). Line colors refer to four 1 kHz frequency bands (0-1 kHz, 1-2 kHz, 4-5 khz, 7-8 kHz). For a propagation independent of ambient sound, frequency bands curves should appear as concentric circles.
* [`figure_I1abcd_propagation_lipaugus.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_I1abcd_propagation_lipaugus.R) : Spectrograms of the Lipaugus vociferans song in the neotropical rainforest (French Guiana) at increasing distances of propagation (200 m, 300 m, 400 m, 600 m and 800 m). Original raw audio recording (experimental data) was used to compute the spectrogram at 200 m while the other spectrograms were computed after applying attenuation law corresponding to the distance of propagation (simulated data)
* [`figure_I2_max_distance_lipaugus.R`](https://github.com/shaupert/haupert_MEE_2022/blob/master/scripts/figure_I2_max_distance_lipaugus.R) : Active space obtained by the model in a neotropical rainforest (French Guiana) for the screaming piha (Lipaugus vociferans) which calls at 111.5dB (SPL re20µPa @1m). The active space may depend on the ambient sound level as well as geometric or spreading loss attenuation (Ageo), atmospheric absorption (Aatm) and habitat or environmental attenuation (Ahab). ambient sound level Ln, T°C, RH, PS and a0 values are set to the experimental values found during the French Guiana transect

## Google Colab
<p align="center">
  <img src="https://s2.qwant.com/thumbr/474x190/f/9/aae347431a927c9b5deb63431ea29c0dd6fceb9210443fdd6bb9b3dba23146/th.jpg?u=https%3A%2F%2Ftse2.mm.bing.net%2Fth%3Fid%3DOIP.IVRAF7_KdEVWUFq1wmDvmQHaC-%26pid%3DApi&q=0&b=1&p=0&a=0g" />
</p>

Google Colab is a very convenient way to test a code written in Python or R without having to install anything on your computer. Everything is run iin the cloud, directly on computers in Google. The required libraries are installed on the fly on the distant computers.

For an interactive way to test our code, we provide most of the R scripts from this repository as Notebook ready to be used on Google Colab. Some adjustment had 
to be done to be able to run the scripts as some R libraries such as Plotly for R don't work properly on Google Colab. In such a case, we mixed R and Python. R 
is used for data processing while Python is used to build and display the figures.

All the notebooks are available [here](https://drive.google.com/drive/folders/1p_xJDaCP2ynVswfaWJLCYIJLMLw9NRic?usp=sharing). 

## Integration of the propagation functions into libraries

For convenience, we added the propagation functions into the library [seewave](https://rug.mnhn.fr/seewave/)  (R language) and into the package [scikit-maad](https://scikit-maad.github.io/) (Python language)

<p align="center">
  <img src="https://rug.mnhn.fr/seewave/PICT/seewave_logo.png" />
</p>
<p align="center">
  <img src="https://scikit-maad.github.io/_images/maad_key_visual_black.png" />
</p>

