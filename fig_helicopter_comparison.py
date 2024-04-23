import matplotlib.pyplot as plt
import pandas as pd
import xmltodict, obspy, glob
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
import sys
base_dir = '/home/jake/Dropbox/fire_infrasound'
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
for station in ['CHKB', 'JDNA', 'JDNB', 'JDSA', 'JDSB', 'QST', 'TOP', 'VPT']:
    t1_t2 = get_t1_t2(station)
    heli_track = pd.read_csv('helicopter_data/heli_track.csv')
    bf_fn = f'beamform_results/10-06T{t1_t2[0]}_10-06T{t1_t2[1]}_{station}_fl4_fh8_winlen60_winfrac1.csv'
    bf = pd.read_csv(bf_fn)
    w_heli = np.where(heli_track[f'{station}_dist'] < 10000)[0]
    
    plt.figure(figsize = (9, 6.5))
    plt.subplot(2,1,1)
    yticks = [0, 90, 180, 270, 360]
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    w_inf = np.where(bf.slowness < 3.5)[0]
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    
    plt.xlim([6, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival ($\degree$ from N)')
    plt.legend([l1, l2], ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'4-8 Hz infrasound backazimuth, {station} array')# (40 sensors, 210 m x 130 m)
    
    plt.subplot(2,1,2)
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    plt.xlim([13, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival ($\degree$ from N)')
    plt.legend((l1, l2), ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'4-8 Hz infrasound backazimuth, {station} array (zoom-in)')# (40 sensors, 210 m x 130 m)
    #annotate_plot()
    plt.tight_layout()
    plt.savefig(f'figures/helicopter_comparison_{station}_fl4_fh8.png')
    #%% compare to 25-32 Hz band
    
for station in ['TOP']:
    t1_t2 = get_t1_t2(station)
    heli_track = pd.read_csv('helicopter_data/heli_track.csv')
    bf_fn = f'beamform_results/10-06T{t1_t2[0]}_10-06T{t1_t2[1]}_{station}_fl25_fh32_winlen60_winfrac1.csv'
    bf = pd.read_csv(bf_fn)
    w_heli = np.where(heli_track[f'{station}_dist'] < 10000)[0]
    
    plt.figure(figsize = (9, 6.5))
    plt.subplot(2,1,1)
    yticks = [0, 90, 180, 270, 360]
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    w_inf = np.where(bf.slowness < 3.5)[0]
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    
    plt.xlim([6, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival (degrees)')
    plt.legend([l1, l2], ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'4-8 Hz infrasound backazimuth, {station} array')# (40 sensors, 210 m x 130 m)
    
    plt.subplot(2,1,2)
    plt.yticks(yticks)
    for y in yticks:
        plt.axhline(y, color = 'gray', linestyle = '--')
    l1, = plt.plot(time_hours(heli_track['t'][w_heli]), heli_track[f'{station}_az'][w_heli], 'k.')
    l2, = plt.plot((bf.t[w_inf] % 1)*24-6, bf.backazimuth[w_inf] % 360, 'r.')
    plt.xlim([13, 18])
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Direction of Arrival (degrees)')
    plt.legend((l1, l2), ['Helicopter', 'Infrasound'], loc = 'upper right', framealpha = 1)
    plt.title(f'4-8 Hz infrasound backazimuth, {station} array (zoom-in)')# (40 sensors, 210 m x 130 m)
    #annotate_plot()
    plt.tight_layout()
    
#%%%

from fire_day_utils import *
#from  import beam_stack_spectrum, beam_stack_spectrogram, apply_function_windows, calc_num_windows
import cleanbf, pickle
def calc_meanfreq(specgram, freqs):
    return np.einsum('ij,j->i', specgram, freqs) / specgram.sum(1)

station = 'TOP'
t1_t2 = get_t1_t2(station)
t1 = obspy.UTCDateTime(f'2023-10-06T{t1_t2[0]}:00')
t2 = obspy.UTCDateTime(f'2023-10-06T{t1_t2[1]}:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
#%%
st = obspy.read(path).merge()
st.trim(t1-20, t2+20) # pad to eliminate potential filter artifacts
st = st.select(station=f'{station}*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333, corners = 6)
st.trim(t1, t2)
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')
cleanbf.add_inv_coords(st, inv)
S = pd.read_csv(f'beamform_results/{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl25_fh32_winlen60_winfrac1.csv')

#%%
for welch_ratio in [5]:#, 10, 15, 20, 30]: 
    filename = f'beamform_results/beamstack_{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl25_fh32_welch{welch_ratio}.pkl'
    if True:
        sg = beam_stack_spectrogram(st, S.backazimuth+180, S.slowness, win_length_sec = 60, welch_ratio=welch_ratio)
        with open(filename, 'wb') as f:
            pickle.dump(sg, f)
    else:
        with open(filename, 'rb') as f:
            sg = pickle.load(f)
    print(sg['freqs'][0,0])
#%%
w = (sg['freqs'][0,:] > 25 ) & (sg['freqs'][0,:] < 32)
specgram = sg['power'][:,w]
freqs = sg['freqs'][0,w]
meanfreq = calc_meanfreq(sg['semblance'][:,w], freqs)
peakfreq = freqs[specgram.argmax(1)]
#%%
import riversound
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'

st = obspy.read(path).merge()
st.trim(t1-20, t2+20) # pad to eliminate potential filter artifacts
st = st.select(station=f'{station}*')[2:3]
#%%
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 20, corners = 6)
st.trim(t1, t2)
sg = riversound.spectrum(st[0], criterion_function = None)
wf = (sg['freqs'] > 25) & (sg['freqs'] < 32)
wt = np.arange(4850,5850)

#%%
specgram = sg['specgram'][wf,wt]
freqs = sg['freqs'][wf]
times = sg['times'][wt]/86400
peakfreq = freqs[specgram.argmax(0)]

#%%
plt.figure()
plt.subplot(2,1,1)
cleanbf.image(np.log10(specgram[1:,:].T), ((times % 1) * 24)-6, np.log10(freqs[1:]), crosshairs = False)
#plt.ylim(-1.0, 1.65)
#plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Spectrogram (beam-stack power)')
plt.ylabel('Frequency (Hz)')
plt.plot(((times % 1) * 24)-6, np.log10(peakfreq))

    #%%
plt.subplot(2,1,2)
cleanbf.image(np.log10(specgram[wf,wt].T), ((times % 1) * 24)-6, np.log10(freqs[wf]), crosshairs = False)
plt.title('Spectrogram (beam-stack power)')
plt.ylabel('Frequency (Hz)')


plt.tight_layout()



#%% calculate helicopter range
# r = k/A
# dr/dt = c * (1-f/f0)
A = np.sqrt(specgram.sum(0))
plt.figure()
plt.subplot(2,1,1)
plt.plot((peakfreq[wt]-27.5)/27.5)
plt.subplot(2,1,2)
plt.semilogy(1/A[wt])
plt.figure()
plt.plot(peakfreq[wt[1:]], np.diff(1/A[wt])*A[wt[1:]], 'k.')
plt.axhline(0)
