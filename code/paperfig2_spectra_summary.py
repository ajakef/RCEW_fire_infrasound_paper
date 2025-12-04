import cleanbf
import obspy, pickle
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import riversound
base_dir = '..'
sys.path.append(f'{base_dir}/code')
from fire_day_utils import get_t1_t2, get_semblance_ticks, get_power_ticks
#%%
arrays = ['JNA', 'JSB', 'QST', 'TOP']#, 'CHKB', 'VPT']
array_reps = ['JNA2', 'JSB1', 'QST01', 'TOP01']#, 'CHKB1', 'VPTA.01']

t1 = obspy.UTCDateTime('2023-10-06_10:00:00')
t2 = obspy.UTCDateTime('2023-10-06_23:59:59')
starttime_morning = obspy.UTCDateTime('2023-10-06_12:00')
starttime_noon = obspy.UTCDateTime('2023-10-06_18:00')

st = obspy.Stream()
for station in array_reps:
    tr = obspy.read(f'{base_dir}/data/infrasound/2023-10-06*{station}*DF.mseed')[0]
    tr.data = tr.data * 3.5012e-3
    st += tr
st_filt = st.copy()
st_filt.filter('highpass', freq = 2)
    
time_interval_length = 1800

colors = ['#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9', '#999999']

def rescale(x):
    x = x - np.mean(x)
    return x/np.max(np.abs(x))

plt.figure(figsize = (7.5, 7))
## A (top): one time series for each array
plt.subplot(2,1,1)
for i, station in enumerate(array_reps):
    tr = st_filt.select(station = station.split('.')[0]).slice(t1, t2)[0]
    t = np.arange(len(tr.data)) * tr.stats.delta / 3600 + 4 # hours MDT
    plt.plot(t, tr.data/5 + len(array_reps)-i, color = colors[i], marker = ',')
    plt.text(4.5, len(array_reps)-i + 0.1, station)

ylim = [0,5]
plt.plot([6,6], ylim, 'k-')
plt.plot([6.5,6.5], ylim, 'k-')
plt.plot([12,12], ylim, 'k:')
plt.plot([12.5,12.5], ylim, 'k:')
plt.yticks([])
plt.plot([10, 10], [0.8, 1.2], 'k-', linewidth = 2)
plt.text(10.2, 1.1, '2 Pa')
plt.ylim(ylim[0], ylim[1])
plt.xlabel('Local Time (hours, MDT)')
plt.title('A. Infrasound Time Series', loc = 'left')

## B (bottom right): Semblance spectra for early morning and noon for each array
plt.subplot(2,2,3)
def F(s): return (s - 1/N) * N/(N-1) # "Adjusted coherence"

for i, station in enumerate(arrays):
    filename = f'{base_dir}/beamform_results/beamstack_10-06T12:00_10-06T23:59_{station}_fl4_fh8_welch5.pkl'
    with open(filename, 'rb') as f:
        sg = pickle.load(f)
    w = np.where((sg['t_mid'] >= starttime_morning) & (sg['t_mid'] <= (starttime_morning + time_interval_length)))[0]
    #N = getN(station, starttime_morning)
    N = sg['N'][0]
    print((station, N))
    plt.semilogx(sg['freqs'][0,:], F(np.median(sg['semblance'][w,:], axis=0)), color = colors[i], linestyle = '-', label = station)
    w = np.where((sg['t_mid'] >= starttime_noon) & (sg['t_mid'] <= (starttime_noon + time_interval_length)))[0]
    plt.semilogx(sg['freqs'][0,:], F(np.median(sg['semblance'][w,:], axis=0)), color = colors[i], linestyle = ':')

plt.plot([0.5, 0.5], [0.5, 0.5], color = 'gray', linestyle = '-', label = '6:00-6:30')
plt.plot([0.5, 0.5], [0.5, 0.5], color = 'gray', linestyle = ':', label = '12:00-12:30')
plt.xlim([0.09, 30])
plt.ylabel('Adj. Semblance (unitless)')
plt.xlabel('Frequency (Hz)')
plt.legend()
plt.title('B. Pre-fire Semblance', loc = 'left')

## C (bottom left): Spectra for early morning and noon for each array
plt.subplot(2,2,4)
for i, station in enumerate(array_reps):
    tr = st.select(station = station).slice(starttime_morning, starttime_morning + time_interval_length)[0]
    s = riversound.spectrum(tr, criterion_function = None, nfft = 4096, window = 'hamming')
    plt.loglog(s['freqs'], s['median'], color = colors[i], linestyle = '-', label = station)
    
    tr = st.select(station = station).slice(starttime_noon, starttime_noon + time_interval_length)[0]
    s = riversound.spectrum(tr, criterion_function = None, nfft = 4096, window = 'hamming')
    plt.loglog(s['freqs'], s['median'], color = colors[i], linestyle = ':')

plt.xlim([0.09, 30])
plt.ylim([1e-7, 1])
plt.ylabel('PSD (Pa$^2$/Hz)')
plt.xlabel('Frequency (Hz)')
plt.title('C. Pre-fire Spectra', loc = 'left')

## below 2 Hz, these spectra basically span the LNM-HNM range. For 2-8 Hz, they are firmly in the middle of the two. Gem self-noise is below actual noise for f<2 Hz (TOP and QST morning), f<20 Hz (others in morning), and always at noon.

plt.tight_layout()
plt.savefig(f'{base_dir}/figures/paperfig_spectra_summary.png', dpi = 300)