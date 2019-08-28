#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:08:17 2019

@author: leon hauser (HAUZERO)
"""
### IMPORT LIBRARIES ###
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import prosail
import pyprosail
from sewar.full_ref import sam as sam

import seaborn as sns
import statsmodels.formula.api as smf
import statsmodels.stats.multicomp as multi

### LOAD DATA ########
### SPECTRA + TRAITS FROM SENTINEL-2 IMAGE #########
data = pd.read_csv('./Export_Output_2.txt')

### MODE OF TRAIT VALUES USED IN SNAP BIOPHYSICAL PROCESSOR ###
modes = pd.read_csv('./atbdmodes.csv')

### WAVELENGTH SENSITIVITY ACROSS BANDS FOR SENTINEL-2A - ONLY 9 BANDS CONSIDERED (10-20m res) ###
spectralsensitivityfile = './S2-SRF_COPE-GSEG-EOPG-TN-15-0007_3.0.csv'
s2sens = pd.read_csv(spectralsensitivityfile)
bandslist = ['B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12']
data['SAM'] = np.nan

### SET DEFAULT VALUES FOR PROSAIL MODELLING #####
#LAI    = modes[modes['Variable'] == 'LAI']['Mode']
#chloro    = modes[modes['Variable'] == 'Cab']['Mode']
N      = modes[modes['Variable'] == 'N']['Mode']
brown    = modes[modes['Variable'] == 'Cpb']['Mode']
brown    = 0
#EWT     = modes[modes['Variable'] == 'LAI']['Mode']
#LMA     = modes[modes['Variable'] == 'Cdm']['Mode']
LMA     = 0.005
irsoil  = modes[modes['Variable'] == 'Soil-BS']['Mode']
hot_spot  = modes[modes['Variable'] == 'Hot']['Mode']
caroten = 8
psoil = 0.5
LIDF = pyprosail.Spherical

### HOW THE PROSAIL COMMANDS WORK ####
#rho_canopy = prosail.run_prosail(n, cab, car, cbrown, cw, cm, lai, lidfa, hspot, tts, tto, psi, \
#                    ant=0.0, alpha=40.0, prospect_version='5', typelidf=2, lidfb=0.0, \
#                    factor='SDR', rsoil0=None, rsoil=None, psoil=None, \
#                    soil_spectrum1=None, soil_spectrum2=None)
#pyprosail = N, chloro, caroten, brown, EWT, LMA, psoil, LAI, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF

### ITERATE PROSAIL MODELLING FOR EACH PIXEL IN THE DATA FILE. IT READS THE LAI, CHLORO AND EWT TRAIT COMBINATIONS FROM THE FILE ###
### GENERATES SPECTRA USING PROSAIL. COMPARES WITH THE ACTUAL SPECTRA FROM THE SENTINEL-2 SCENCE USING SPECTRAL ANGLE MAPPER ###
### OPTIONS TO PLOT, PAUSE AND PRINT RESULTS ARE CURRENTLY COMMENTED ### 
for irow in range(np.shape(data)[0]):
        rowvars =  data.iloc[irow]
        if rowvars['flagsbyfiveall_1'] == 1: 
            if rowvars['B8A_1'] > 0.1:
                LAI             = rowvars['LAI_1']
                chloro          = rowvars['CAB_1']*rowvars['LAI_1']
                EWT             = rowvars['CWC_1']*rowvars['LAI_1']
                solar_zenith    = rowvars['sun_zenith_1']
                solar_azimuth   = rowvars['sun_azimuth_1']
                view_zenith     = rowvars['view_zenith_mean_1'] 
                view_azimuth    = rowvars['view_azimuth_mean_1']
                specout         = pyprosail.run(N, chloro, caroten, brown, EWT, LMA, psoil, LAI, hot_spot, solar_zenith, solar_azimuth, view_zenith, view_azimuth, LIDF)
#                print LAI, chloro, EWT, solar_zenith, solar_azimuth, view_zenith, view_azimuth
                s2out = pd.DataFrame(columns={'B2','B3','B4', 'B5', 'B6', 'B7', 'B8A', 'B11', 'B12'})
                s2outt = np.array([[520,560,665,705,740,783,865,1610,2190],np.zeros(9)])
                for ib, band in enumerate(bandslist):
                    weigthedvalue = 0
                    bandname = 'S2A_SR_AV_' + band
                    sumvalue        = np.sum(s2sens[bandname])
                    for ix in np.where(s2sens[bandname] != 0)[0]:
                        iwl = np.where((specout[:,0]*1000) == int(s2sens['SR_WL'].iloc[ix]))[0][0]
                        weigthedvalue   = weigthedvalue + (float(s2sens[bandname].iloc[ix]) * specout[iwl,1])
                    s2out['band'] = (weigthedvalue/sumvalue)     
                    s2outt[1,ib] =  (weigthedvalue/sumvalue)
#                    print band, (weigthedvalue/sumvalue)
                specin          = np.array([[520,560,665,705,740,783,865,1610,2190],[rowvars['B2_1'], rowvars['B3_1'], rowvars['B4_1'], rowvars['B5_1'], rowvars['B6_1'], rowvars['B7_1'], rowvars['B8A_1'], rowvars['B11_1'], rowvars['B12_1']]])
#                plt.plot(s2outt[0], s2outt[1]); plt.plot((specout[:,0]*1000),specout[:,1]); plt.plot(specin[0], specin[1])
#                plt.legend()
#                plt.show()
                samerror = sam(specin, s2outt)
                data['SAM'].iloc[irow] = samerror
#                bla = raw_input("Press any key")
                if irow in range(0,22251,222):
                    print float(irow/np.shape(data)[0])
         
### MEANS PER LAND USE ####            
np.nanmean(data['SAM'][data['lcselectedint'] == 1])
np.nanmean(data['SAM'][data['lcselectedint'] == 2])
np.nanmean(data['SAM'][data['lcselectedint'] == 6])
np.nanmean(data['SAM'][data['lcselectedint'] == 7])
np.nanmean(data['SAM'][data['lcselectedint'] == 8])

### SAVE RESULTS WITH OUT NAN-VALUES
data = data[data['SAM'] >= 0]
outname = '/media/leon/FREECOM HDD/Data/Scaling-Borneo/Forwardmodelling/samdata.csv'
data.to_csv(outname, index=False)

### LOAD RESULTS AGAIN ###
data = pd.read_csv('/media/leon/FREECOM HDD/Data/Scaling-Borneo/Forwardmodelling/samdata.csv')

### PLOTS SAM DISTRIBUTIONS PER LAND USE TYPE ###
unique_classes  =   [1, 2, 6, 7, 8]
colors          =   ['green','skyblue','pink','brown', 'red', 'black']
classnames      =   ["IntactForest", "LoggedForest", "OilPalmPlan", "TimberPlan", "MosaicCropland"]

fig = plt.figure(figsize=(8,8))
plt.xlabel('SAM', fontsize=25)
plt.ylabel('KDE', fontsize=25)
for i,lc in enumerate(unique_classes):    
        traitlc = data['SAM'][data['lcselectedint'] == lc]
        ax1 = sns.kdeplot(traitlc, linewidth=2.5, color=colors[i], label=classnames[i])
plt.tight_layout()
plt.legend()
plt.show()

### RUN AVONA AND TUKEY HSD TEST OF SAM RESULTS AGAINST LAND USE ###
anova = smf.ols(formula='SAM ~ C(lcselectedint)', data=data).fit()
print anova.summary()
tukey = multi.pairwise_tukeyhsd(data['SAM'],data['lcselectedint'])
print tukey.summary()



