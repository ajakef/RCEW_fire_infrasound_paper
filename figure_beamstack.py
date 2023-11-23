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
import sys
base_dir = '/home/jake/Dropbox/fire_infrasound'
sys.path.append(f'{base_dir}/code')
from fire_day_utils import get_t1_t2, get_semblance_ticks, get_power_ticks

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

base_dir = '/home/jake/Dropbox/fire_infrasound/'
#%% make a back-azimuth/beamstack figure for a station

station = 'VPT'



#semblance_ticks = {'TOP':[0.03, 0.1, 0.3, 1], 'QST': [0.1, 0.3, 1], 'VPT':[0.1, 0.3, 1]}.get(station, [0.25, 0.5, 1])
#t1_t2 = {'CHKB':['15:00', '23:59'], 'VPT':['18:00','22:50']}.get(station, ['12:00', '23:59'])
semblance_ticks = get_semblance_ticks(station)
t1_t2 = get_t1_t2(station)
power_ticks = get_power_ticks(station)

bf_fn = f'{base_dir}/beamform_results/10-06T{t1_t2[0]}_10-06T{t1_t2[1]}_{station}_fl4_fh8_winlen60_winfrac1.csv'

bf = pd.read_csv(bf_fn)
    
fig, axes = plt.subplots(
    3,2, sharex='col', figsize=(6.5, 8),
    gridspec_kw = {'width_ratios' : (20,1)},
)

yticks = [0, 90, 180, 270, 360]
axes[0,0].set_yticks(yticks, labels = ['N', 'E', 'S', 'W', 'N'])
for y in yticks:
    axes[0,0].axhline(y, color='gray', linestyle='--')
w = np.where(bf.slowness < 3.5)[0]
x = ((bf['t'] % 1) * 24)-6
xlim = [x.min(), x.max()]
axes[0,0].plot(x[w], bf.backazimuth[w] % 360, color = 'darkred', marker = '.', linestyle = 'none')

axes[0,0].set_xlim(xlim)
#axes[0].xlabel('Time (MDT)')
axes[0,0].set_ylabel('Direction of Arrival')
axes[0,0].set_title(f'4-8 Hz Direction of Arrival')

#second_axis = axes[0,0].twinx()
#second_axis.plot((bf.t-bf.t[0])*24+6, bf.backazimuth % 360, alpha=0)
#second_axis.set_yticks(yticks)
#second_axis.set_ylabel('Direction of Arrival (degrees)')

axes[0,1].remove() # remove colorbar slot for top panel


filename = f'{base_dir}/beamform_results/beamstack_10-06T{t1_t2[0]}_10-06T{t1_t2[1]}_{station}_fl4_fh8_welch5.pkl'
with open(filename, 'rb') as f:
    sg = pickle.load(f)

im = cleanbf.image(np.log10(sg['power']), x, np.log10(sg['freqs'][0,:]), 
              crosshairs = False, ax = axes[1,0], zmax = power_ticks[-1], zmin = power_ticks[0])
cbar = plt.colorbar(im, cax = axes[1,1])
cbar.set_label('log power (arbitrary units)')
cbar.set_ticks(power_ticks)




axes[1,0].set_ylim(-1.0, 1.65)
axes[1,0].set_xlim(xlim)

axes[1,0].set_yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
axes[1,0].set_title('Beam-Stack Power Spectrogram')
axes[1,0].set_ylabel('Frequency (Hz)')
cbar.set_label('log power (arbitrary units)')

for y in [1, 10]:
    axes[1,0].axhline(np.log10(y), color = 'gray', linestyle = '--')

im = cleanbf.image(np.log10(np.abs(sg['semblance'])), x, np.log10(sg['freqs'][0,:]), 
              crosshairs = False, zmin = np.log10(semblance_ticks[0]), zmax = 0, ax = axes[2,0])
cbar = plt.colorbar(im, cax = axes[2,1])
#ticks = np.array([0.03, 0.1, 0.3, 1])
cbar.set_ticks(np.log10(semblance_ticks))
cbar.set_ticklabels(semblance_ticks)
axes[2,0].set_ylim(-1.0, 1.65)
axes[2,0].set_xlim(xlim)
axes[2,0].set_yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
axes[2,0].set_title('Beam-Stack Semblance Spectrogram')
axes[2,0].set_xlabel('Local Time (hours)')
axes[2,0].set_ylabel('Frequency (Hz)')
cbar.set_label('semblance (unitless)')
#for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
#    axes[2].axvline(x)

#for y in [0.3, 1, 3, 10, 30]:
for y in [1, 10]:
    axes[2,1].axhline(np.log10(y), color = 'gray', linestyle = '--')

fig.tight_layout()

fig.savefig(f'{base_dir}/figures/beamstack_{station}.png')
#####################################################
#%%



















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