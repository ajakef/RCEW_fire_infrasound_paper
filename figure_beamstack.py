#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 14:47:52 2023

@author: jake
"""

import cleanbf
import obspy, pickle
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_bf_result(b, style = 'k.'):
    if type(b) is str:
        b = pd.read_csv(b)
    t = (b.iloc[:,0] - np.floor(b.iloc[0,0]))*24-6
    titles = ['backazimuth', 'horiz. slowness']
    for i in range(1, 3):
        plt.subplot(2,1,i)
        plt.plot(t, b.iloc[:,i+2], style)
        for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
            plt.axvline(x)
        plt.title(titles[i-1])

#%% TOP

station = 'TOP'
plt.figure(figsize = (6.5, 8))
plt.subplot(3,1,1)

bf_fn = f'beamform_results/10-06T12:00_10-06T23:59_{station}_fl4_fh8_winlen60_winfrac1.csv'
bf = pd.read_csv(bf_fn)

yticks = [0, 90, 180, 270, 360]
plt.yticks(yticks)
for y in yticks:
    plt.axhline(y, color='gray', linestyle='--')
plt.plot((bf.t-bf.t[0])*24+6, bf.backazimuth % 360, color = 'darkred', marker = '.', linestyle = 'none')
plt.xlim([6, 18])
#plt.xlabel('Time (MDT)')
plt.ylabel('Backazimuth (degrees)')
plt.title(f'4-8 Hz infrasound backazimuth, {station} array (40 sensors, 210 m x 130 m)')

#%%
filename = f'beamform_results/beamstack_10-06T12:00_10-06T23:59_{station}_fl4_fh8_welch5.pkl'
with open(filename, 'rb') as f:
    sg = pickle.load(f)

plt.subplot(3,1,2)
cleanbf.image(np.log10(sg['power']), ((bf['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False)
cbar = plt.colorbar()
cbar.set_label('log power (arbitrary units)')
#ticks = np.array([0.03, 0.1, 0.3, 1])
#cbar.set_ticks(np.log10(ticks))
#cbar.set_ticklabels(ticks)
plt.ylim(-1.0, 1.65)
plt.xlim([6, 18])
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Beam-Stack Power Spectrogram')
plt.ylabel('Frequency (Hz)')
cbar.set_label('log power (arbitrary units)')
#for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
#    plt.axvline(x)

#for y in [0.3, 1, 3, 10, 30]:
for y in [1, 10]:
    plt.axhline(np.log10(y), color = 'gray', linestyle = '--')

plt.subplot(3,1,3)
cleanbf.image(np.log10(np.abs(sg['semblance'])), ((bf['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False, zmin = -1.5, zmax = 0)
cbar = plt.colorbar()
ticks = np.array([0.03, 0.1, 0.3, 1])
cbar.set_ticks(np.log10(ticks))
cbar.set_ticklabels(ticks)
plt.ylim(-1.0, 1.65)
plt.xlim([6, 18])
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Beam-Stack Semblance Spectrogram')
plt.xlabel('Local Time (hours)')
plt.ylabel('Frequency (Hz)')
cbar.set_label('semblance (unitless)')
#for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
#    plt.axvline(x)

#for y in [0.3, 1, 3, 10, 30]:
for y in [1, 10]:
    plt.axhline(np.log10(y), color = 'gray', linestyle = '--')

plt.tight_layout()
#%%
plt.savefig(f'figures/beamstack_{station}.png')
#####################################################
#%% QST
station = 'QST'
plt.figure(figsize = (6.5, 8))
plt.subplot(3,1,1)

bf_fn = f'beamform_results/10-06T12:00_10-06T23:59_{station}_fl4_fh8_winlen60_winfrac1.csv'
bf = pd.read_csv(bf_fn)

yticks = [0, 90, 180, 270, 360]
plt.yticks(yticks)
for y in yticks:
    plt.axhline(y, color='gray', linestyle='--')
plt.plot((bf.t-bf.t[0])*24+6, bf.backazimuth % 360, color = 'darkred', marker = '.', linestyle = 'none')
plt.xlim([6, 18])
#plt.xlabel('Time (MDT)')
plt.ylabel('Backazimuth (degrees)')
#lt.title(f'4-8 Hz infrasound backazimuth, {station} array (40 sensors, 210 m x 130 m)')

#%%
filename = f'beamform_results/beamstack_10-06T12:00_10-06T23:59_{station}_fl4_fh8_welch5.pkl'
with open(filename, 'rb') as f:
    sg = pickle.load(f)

plt.subplot(3,1,2)
cleanbf.image(np.log10(sg['power']), ((bf['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False)
plt.ylim(-1.0, 1.65)
plt.xlim([6, 18])
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Beam-Stack Power Spectrogram')
plt.ylabel('Frequency (Hz)')
cbar = plt.colorbar()
cbar.set_label('log power (arbitrary units)')
#for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
#    plt.axvline(x)

#for y in [0.3, 1, 3, 10, 30]:
for y in [1, 10]:
    plt.axhline(np.log10(y), color = 'gray', linestyle = '--')

plt.subplot(3,1,3)
cleanbf.image(np.log10(np.abs(sg['semblance'])), ((bf['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False, zmin = -1.5, zmax = 0)
plt.ylim(-1.0, 1.65)
plt.xlim([6, 18])
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Beam-Stack Semblance Spectrogram')
plt.xlabel('Local Time (hours)')
plt.ylabel('Frequency (Hz)')
cbar = plt.colorbar()
ticks = np.array([0.03, 0.1, 0.3, 1])
cbar.set_ticks(np.log10(ticks))
cbar.set_ticklabels(ticks)

#for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
#    plt.axvline(x)

#for y in [0.3, 1, 3, 10, 30]:
for y in [1, 10]:
    plt.axhline(np.log10(y), color = 'gray', linestyle = '--')

plt.tight_layout()
#%%
plt.savefig(f'figures/beamstack_{station}.png')