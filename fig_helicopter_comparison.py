import matplotlib.pyplot as plt
import pandas as pd
import xmltodict, obspy, glob
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

def time_hours(t, t1 = obspy.UTCDateTime('2023-10-06:06:00')):
    return np.array([(obspy.UTCDateTime(tt) - t1)/3600 for tt in t])

heli_track = pd.read_csv('helicopter_data/heli_track.csv')
bf_fn = 'beamform_results/fire_TOP_win60_fl0.5_fh2.csv'
bf = pd.read_csv(bf_fn)

w = np.where(heli_track['TOP_dist'] < 10000)[0]
plt.plot(time_hours(heli_track['t'][w]), heli_track['TOP_az'][w], 'k.')
plt.plot((bf.t-bf.t[0])*24-6, bf.backazimuth % 360, 'r.')
plt.xlim([7, 18])
plt.xlabel('Time (MDT)')
plt.ylabel('Backazimuth (degrees)')
plt.legend(['Helicopter', 'Infrasound'])
plt.title('0.5-2 Hz infrasound backazimuth, TOP array (40 sensors, 210 m x 130 m)')
eps = 10/60
plt.text(13+12/60+eps, 360, 'Main ignition begins', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
plt.text(14+23/60+eps, 360, 'Helicopter leaves', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
plt.text(15+12/60+eps, 360, 'Helicopter returns', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
plt.text(16+41/60+eps, 360, 'Helicopter leaves', rotation = 'vertical', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
plt.text(7.25, 360, 'Distant infrasound\n (quiet morning)', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
plt.text(11, 360, 'Mid-day \nturbulence', horizontalalignment = 'left', verticalalignment = 'top', bbox = dict(facecolor='white', alpha = 0.6, linewidth=0))
for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
    plt.axvline(x)
