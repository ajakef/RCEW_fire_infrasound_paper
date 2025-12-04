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
base_dir = '..'
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

#%% loop through stations and make a back-azimuth/beamstack figure for each
for station in ['JNA', 'JNB', 'JSA', 'JSB', 'QST', 'TOP']: #, 'CHKB', 'VPT']
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
    axes[0,0].set_yticks(yticks)
    for y in yticks:
        axes[0,0].axhline(y, color='gray', linestyle='--')
    w = np.where(bf.slowness < 3.5)[0]
    x = ((bf['t'] % 1) * 24)-6
    xlim = [x.min(), x.max()]
    axes[0,0].plot(x[w], bf.backazimuth[w] % 360, color = 'darkred', marker = '.', linestyle = 'none')
    
    axes[0,0].set_xlim(xlim)
    axes[0,0].set_ylabel('Direction of Arrival ($\degree$ from N)')
    axes[0,0].set_title(f'4-8 Hz Direction of Arrival')

    
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

    for y in [1, 10]:
        axes[2,1].axhline(np.log10(y), color = 'gray', linestyle = '--')
    
    fig.tight_layout()
    
    fig.savefig(f'{base_dir}/figures/beamstack_{station}.png')
    print(f'{base_dir}/figures/beamstack_{station}.png')