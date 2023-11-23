import cleanbf
import obspy
from obspy import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import sys
base_dir = '/home/jake/Dropbox/fire_infrasound'
sys.path.append(f'{base_dir}/code')
from fire_day_utils import get_t1_t2, get_semblance_ticks

def beam_stack_spectrum(st, sx, sy, win_length_sec = 10):
    print(f'({sx:0.4f}, {sy:0.4f})')
    cross_spec, FT, freqs, dfn, dfd = cleanbf.calc_cross_spectrum(st, win_length_sec = win_length_sec)
    weight = cleanbf.calc_weights(cleanbf.make_steering_vectors(st, freqs, [sx], [sy]))[:,0,0,:]
    P = np.abs(np.einsum('ki,ijk,kj->k', weight.conj(), cross_spec, weight)) # i, j: stations; k: freq.
    semblance = P/np.einsum('iik->k', cross_spec)
    return {'freqs':freqs, 'power':P, 'semblance':semblance}

def beam_stack_spectrogram(st, forward_az, slowness, win_length_sec = 40, welch_ratio = 4, overlap = 0):
    num_windows = calc_num_windows(st[0].stats.endtime - st[0].stats.starttime, win_length_sec, overlap)
    sx = np.sin(forward_az * np.pi/180) * slowness
    sy = np.cos(forward_az * np.pi/180) * slowness
    if (len(sx) > 1) and (len(sx) != num_windows):
        raise Exception(f'sx and sy are the wrong length (should be 1 or {num_windows})')
                                
    f_inputs = []
    for sxx, syy in zip(sx, sy):
        f_inputs.append([sxx, syy, win_length_sec / welch_ratio])
    return apply_function_windows(st, beam_stack_spectrum, win_length_sec, overlap = overlap, f_inputs = f_inputs)

def apply_function_windows(st, f, win_length_sec, overlap = 0.5, f_inputs = [], verbose = True):
    """
    Run an analysis (or suite of analyses) on overlapping windows for some dataset
    
    Parameters:
    -----------
    st : obspy.Stream
    Stream including data to divide into windows and analyze

    f : function
    Accepts single variable "st" (obspy.Stream), returns dictionary of results

    win_length_sec : float
    Length of analysis windows [seconds]

    overlap : float
    Proportion of a window that overlaps with the previous window [unitless, 0-1]
    f_inputs: list
    List of inputs to provide to f() in addition to st (1 for each time window)

    Returns:
    --------
    dictionary with following items:
    --t_mid (obspy.UTCDateTime): mean time of each window
    --all elements of the output of "f", joined into numpy arrays

    Note:
    -----
    For each time window, the stream is trimmed to fit, but not tapered, detrended, or otherwise 
    processed. If those steps are necessary, be sure they are carried out in f().
"""
    # f must input a stream and return a dict
    eps = 1e-6
    t1 = st[0].stats.starttime
    t2 = st[0].stats.endtime
    data_length_sec = t2 - t1
    #num_windows = 1 + int(np.ceil((data_length_sec - win_length_sec) / (win_length_sec * (1 - overlap)) - eps))
    num_windows = calc_num_windows(data_length_sec, win_length_sec, overlap)
    if (len(f_inputs) > 0) and (len(f_inputs) != num_windows):
        raise Exception(f'len(f_inputs) is not the same as num_windows {num_windows}')
    print(num_windows)
    for i in range(num_windows):
        if verbose and ((i % (num_windows//100)) == 0):
            print(f'{i} of {num_windows}')
        win_start = t1 + i*(data_length_sec - win_length_sec) / (num_windows-1)
        st_tmp = st.slice(win_start-eps, win_start + win_length_sec - eps, nearest_sample = False)
        if len(f_inputs) > 0:
            win_dict = f(st_tmp, *(f_inputs[i]))
        else:
            win_dict = f(st_tmp)
        if i == 0:
            output_dict = {key:[] for key in win_dict.keys()}
            output_dict['t_mid'] = []
        for key in win_dict.keys():
            output_dict[key].append(win_dict[key])
        output_dict['t_mid'].append(win_start + win_length_sec/2)
    output_dict = {key:np.array(output_dict[key]) for key in output_dict.keys()}
    return output_dict

def calc_num_windows(data_length_sec, win_length_sec, overlap):
    eps = 1e-3
    return 1 + int(np.floor((data_length_sec - win_length_sec) / (win_length_sec * (1 - overlap)) - eps))

#%% define stream to process
station = 'VPT'
t1_t2 = get_t1_t2(station)
t1 = UTCDateTime(f'2023-10-06T{t1_t2[0]}:00')
t2 = UTCDateTime(f'2023-10-06T{t1_t2[1]}:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
st.trim(t1-20, t2+20) # pad to eliminate potential filter artifacts
st = st.select(station=f'{station}*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333, corners = 6)
st.trim(t1, t2)
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')
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

S = pd.read_csv(f'beamform_results/{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl4_fh8_winlen60_winfrac1.csv')

## Following test shows using a higher welch ratio (shorter window) loses frequency resolution but 
## doesn't visibly improve noise reduction. So, stick with welch ratio 5
for welch_ratio in [5]:#, 10, 15, 20, 30]: 
    filename = f'beamform_results/beamstack_{t1.strftime("%m-%dT%H:%M")}_{t2.strftime("%m-%dT%H:%M")}_{station}_fl4_fh8_welch{welch_ratio}.pkl'
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

