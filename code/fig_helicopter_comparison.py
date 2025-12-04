import matplotlib.pyplot as plt
import pandas as pd
import obspy
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
import sys
base_dir = '..' # run this script inside the 'code/' folder
sys.path.append(f'{base_dir}/code')
from fire_day_utils import get_t1_t2, get_semblance_ticks

def time_hours(t, t1 = obspy.UTCDateTime('2023-10-06:06:00')):
    return np.array([(obspy.UTCDateTime(tt) - t1)/3600 for tt in t])

def annotate_plot(text = False):
    eps = 10/60
    if text:
        plt.text(13+12/60+eps, 360, 'Main ignition begins', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
        plt.text(14+23/60+eps, 360, 'Helicopter leaves', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
        plt.text(15+12/60+eps, 360, 'Helicopter returns', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
        plt.text(16+41/60+eps, 360, 'Helicopter leaves', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
        plt.text(7.25, 360, 'Distant infrasound\n (quiet morning)', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
        plt.text(11, 360, 'Mid-day \nturbulence', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
    for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
        plt.axvline(x)
#%%
for station in ['JNA', 'JNB', 'JSA', 'JSB', 'TOP']: #, 'QST', 'VPT', 'CHKB']
    t1_t2 = get_t1_t2(station)
    heli_track = pd.read_csv(f'{base_dir}/data/helicopter_data/helicopter_station_azimuths.csv')
    bf_fn = f'{base_dir}/beamform_results/10-06T{t1_t2[0]}_10-06T{t1_t2[1]}_{station}_fl4_fh8_winlen60_winfrac1.csv'
    bf = pd.read_csv(bf_fn)
    w_heli = np.where(heli_track[f'{station}_dist'] < 10000)[0]
    
    plt.figure(figsize = (7.5, 8.5))
    plt.subplot(2,1,1)
    yticks = [0, 90, 180, 270, 360]
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    w_inf = np.where(bf.slowness < 3.5)[0]
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    if station == 'TOP':
        plt.fill_between([6, 18], [-90, -90], [450,450], color = 'lightgray')
        plt.fill_between([6, 18], [33, 33], [187,187], color = 'white')
    plt.ylim([0,360])
    plt.xlim([6, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival ($\degree$ from N)')
    plt.legend([l1, l2], ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'A. 4-8 Hz infrasound backazimuth, {station} array', loc = 'left')# (40 sensors, 210 m x 130 m)
    
    plt.subplot(2,1,2)
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    if station == 'TOP':
        plt.fill_between([6, 18], [-90, -90], [450,450], color = 'lightgray')
        plt.fill_between([6, 18], [33, 33], [187,187], color = 'white')
    plt.ylim([0,360])
    plt.xlim([12.45, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival ($\degree$ from N)')
    plt.legend((l1, l2), ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'B. 4-8 Hz infrasound backazimuth, {station} array (zoom-in)', loc = 'left')# (40 sensors, 210 m x 130 m)
    plt.tight_layout()
    plt.savefig(f'{base_dir}/figures/helicopter_comparison_{station}_fl4_fh8.png')
