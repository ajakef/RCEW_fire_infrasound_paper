import cleanbf
import obspy
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import sys
sys.path.append(f'{base_dir}/code')

from fire_day_utils import get_t1_t2, get_semblance_ticks, beam_stack_spectrum, beam_stack_spectrogram, apply_function_windows, calc_num_windows
base_dir = '..' # run this script from the 'code/' folder
#%% define stream to process
station = 'QST'
for station in ['TOP', 'QST', 'JDNA', 'JDNB', 'JDSA', 'JDSB']:
    t1_t2 = get_t1_t2(station)
    t1 = UTCDateTime(f'2023-10-06T{t1_t2[0]}:00')
    t2 = UTCDateTime(f'2023-10-06T{t1_t2[1]}:59')
    path = f'{base_dir}/data/infrasound/2023-10-06*'
    st = obspy.read(path).merge()
    st.trim(t1-20, t2+20) # pad to eliminate potential filter artifacts
    st = st.select(station=f'{station}*')
    st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
    st.filter('highpass', freq = 0.333, corners = 6)
    st.trim(t1, t2)
    inv = obspy.read_inventory(f'{base_dir}/coordinates/RCEW_inventory.xml')
    cleanbf.add_inv_coords(st, inv)


    #%% drop high-noise traces, balancing large number of sensors vs dropping highest-noise sensors
    sd = np.sort(st.std())
    total_power = sd**2
    num_stations = np.arange(len(st))+1
    assumed_signal_power = (sd**2)[0] * 0.8 # for 0.3, 0.5, and 0.8, the result is the same: drop the 4 loudest sensors
    noise_power = np.cumsum(sd**2 - assumed_signal_power)/num_stations # if all sensors have the same power, this will be the same for all of them
    noise_power_after_stack = noise_power/num_stations
    max_sd_to_keep = sd[noise_power_after_stack.argmin()]
    
    w = np.where(np.array(st.std()) <= max_sd_to_keep)[0]
    st = obspy.Stream([st[i] for i in w])
    
    #%% run the beam stack spectrogram
    #for station in ['TOP', 'QST']
    
    S = pd.read_csv(f'{base_dir}/beamform_results/{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl4_fh8_winlen60_winfrac1.csv')
    
    ## scaling check. Note that cleanbf spectra scaling is tested in tests.py.
    if False:
        sg = beam_stack_spectrogram(st, S.backazimuth+180, S.slowness, win_length_sec = 60, welch_ratio=1)
        plt.loglog(sg['freqs'][0,:],  np.mean(sg['power'], axis=0)) # plot mean beam-stack power spectrum
        cs = cleanbf.calc_cross_spectrum(st, win_length_sec = 60)
        for i in range(cs[0].shape[0]):
            plt.loglog(cs[2], cs[0][i,i,:])

    ## Following test shows using a higher welch ratio (shorter window) loses frequency resolution but 
    ## doesn't visibly improve noise reduction. So, stick with welch ratio 5
    #for welch_ratio in [5]:#, 10, 15, 20, 30]:
    welch_ratio = 5
    filename = f'{base_dir}/beamform_results/beamstack_{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl4_fh8_welch{welch_ratio}.pkl'
    if True:
        sg = beam_stack_spectrogram(st, S.backazimuth+180, S.slowness, win_length_sec = 60, welch_ratio=welch_ratio)
        with open(filename, 'wb') as f:
            pickle.dump(sg, f)
    else:
        with open(filename, 'rb') as f:
            sg = pickle.load(f)
    print(sg['freqs'][0,0])
    #%%
    plt.figure()
    plt.subplot(2,1,1)
    cleanbf.image(np.log10(sg['power']), ((S['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False)
    plt.ylim(-1.0, 1.65)
    plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
    plt.title('Spectrogram (beam-stack power)')
    plt.ylabel('Frequency (Hz)')

    plt.subplot(2,1,2)
    cleanbf.image(np.log10(np.abs(sg['semblance'])), ((S['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False, zmin = -1.5, zmax = 0)
    plt.ylim(-1.0, 1.65)
    plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
    plt.title('Spectrogram (beam-stack semblance)')
    plt.xlabel('Local Time (hours)')
    plt.ylabel('Frequency (Hz)')

    plt.tight_layout()
    for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
        plt.axvline(x)

    for y in [0.3, 1, 3, 10, 30]:
        plt.axhline(np.log10(y))

