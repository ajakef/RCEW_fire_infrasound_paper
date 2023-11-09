import matplotlib.pyplot as plt
import pandas as pd
import xmltodict, obspy, glob
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

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

heli_track = pd.read_csv('helicopter_data/heli_track.csv')
bf_fn = 'beamform_results/10-06T12:00_10-06T23:59_TOP_fl4_fh8_winlen60_winfrac1.csv'
bf = pd.read_csv(bf_fn)
w = np.where(heli_track['TOP_dist'] < 10000)[0]

plt.figure(figsize = (9, 6.5))
plt.subplot(2,1,1)
yticks = [0, 90, 180, 270, 360]
plt.yticks(yticks)
for y in yticks:
    plt.axhline(y, color = 'gray', linestyle = '--')
plt.plot(time_hours(heli_track['t'][w]), heli_track['TOP_az'][w], 'k.')
plt.plot((bf.t-bf.t[0])*24+6, bf.backazimuth % 360, 'r.')
plt.xlim([6, 18])
plt.xlabel('Local Time (hours)')
plt.ylabel('Backazimuth (degrees)')
plt.legend(['Helicopter', 'Infrasound'])
plt.title('4-8 Hz infrasound backazimuth, TOP array (40 sensors, 210 m x 130 m)')

plt.subplot(2,1,2)
plt.yticks(yticks)
for y in yticks:
    plt.axhline(y, color = 'gray', linestyle = '--')
plt.plot(time_hours(heli_track['t'][w]), heli_track['TOP_az'][w], 'k.')
plt.plot((bf.t-bf.t[0])*24+6, bf.backazimuth % 360, 'r.')
plt.xlim([13, 18])
plt.xlabel('Local Time (hours)')
plt.ylabel('Backazimuth (degrees)')
plt.legend(['Helicopter', 'Infrasound'])
plt.title('4-8 Hz infrasound backazimuth, TOP array (40 sensors, 210 m x 130 m)')
#annotate_plot()
plt.tight_layout()
plt.savefig('figures/TOP_helicopter_comparison_fl4_fh8.png')