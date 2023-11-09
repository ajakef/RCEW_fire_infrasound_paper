from obspy import UTCDateTime
import obspy, glob
from obspy.signal.array_analysis import array_processing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def store(x, y, offset_samples):
    if (offset_samples % 360000) == 0:
        print(offset_samples * 0.01/3600)
        
def ap_to_df(x, fn):
    y = pd.DataFrame(x, columns = ['t', 'semblance', 'power', 'backazimuth', 'slowness'])
    y.to_csv(fn, index = False)
    return y

def plot_bf_result_all(b, style = 'k.', slowness_max = 4):
    if type(b) is str:
        b = pd.read_csv(b)
    b = b.loc[b.iloc[:,4] < slowness_max, :]
    if b.shape[0] == 0:
        return
    t = (b.iloc[:,0] - np.floor(b.iloc[0,0]))*24-6
    titles = ['semblance', 'power', 'backazimuth', 'horiz. slowness']
    for i in range(1, 5):
        plt.subplot(2,2,i)
        if i in [1,2]:
            plt.semilogy(t, b.iloc[:,i], style)
        else:
            plt.plot(t, b.iloc[:,i] % 360, style)
        for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
            plt.axvline(x)
        plt.title(titles[i-1])

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


def add_inv_coords(st, inv):
    """
    For each trace tr in stream st, adds a 'coordinates' AttribDict to tr.stats using location info
    from an inventory. See application in obspy.signal.array_analysis.array_processing.

    Parameters
    ----------
    st: stream that locations should be added to
    inv: obspy.Inventory or pandas.DataFrame containing locations for all traces in st

    Returns: None, changes st in place
    """
    if type(inv) is obspy.Inventory:
        for tr in st:
            loc = obspy.core.AttribDict(inv.get_coordinates(tr.get_id()))
            tr.stats = obspy.core.Stats({**tr.stats, 'coordinates': loc})
    elif type(inv) is pd.DataFrame:
        for tr in st:
            id = tr.get_id()
            w = np.where(inv.SN == tr.stats.station)[0][0]
            loc = obspy.core.AttribDict(
                {'latitude':inv.loc[w,'lat'],
                 'longitude':inv.loc[w,'lon'],
                 'elevation':0,
                 'local_depth':0}
                )
            tr.stats = obspy.core.Stats({**tr.stats, 'coordinates': loc})
#%%
######################################
######################################

#path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/*TOP*'
path = '/home/jake/2023-10-20_OSIRIS-REx_download/mseed_to_share/*TOP*'
st_all = obspy.read(path)
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st_all, inv)



### OSIRIS-REx: 
t1 = UTCDateTime('2023-09-24T15:00:00')
t2 = UTCDateTime('2023-09-24T16:30:00')

st_OR = st_all.slice(t1, t2)
st_OR.select(station='TOP*')

x10l = array_processing(st_OR, win_len = 10, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 0.5, frqhigh = 5, stime = t1, etime = t2, prewhiten = False)

x30l = array_processing(st_OR, win_len = 30, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 0.5, frqhigh = 5, stime = t1, etime = t2, prewhiten = False)

x30h = array_processing(st_OR, win_len = 30, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 5, frqhigh = 20, stime = t1, etime = t2, prewhiten = False)


for i in range(5):
    plt.subplot(2,3,i+1)
    plt.plot(x30h[:,i], 'k,')


### Helicopter example: high SNR, great beamforming. The array works and the coordinates are good.
t1 = UTCDateTime('2023-10-06T21:30:00')
t2 = UTCDateTime('2023-10-06T22:00:00')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
#path = '/home/jake/2023-10-20_OSIRIS-REx_download/mseed_to_share/*'
st_H = obspy.read(path).merge()
st_H = st_H.slice(t1, t2)
st_H = st_H.select(station='TOP*')
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st_H, inv)


H20 = ap_to_df(array_processing(st_H, win_len = 20, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 10, frqhigh = 30, stime = t1, etime = t2, prewhiten = False))

H20.to_csv('helicopter_win20_fl10_fh30.csv', index = False)

for i in range(5):
    plt.subplot(2,3,i+1)
    plt.plot(H20.iloc[:,i], 'k.')

#################################
#################################
### Burn day, low freq (attempt no helicopter)
t1 = UTCDateTime('2023-10-06T00:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
st = st.slice(t1, t2)
st = st.select(station='TOP*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333)
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st, inv)

fn = 'fire_TOP_win60_fl0.5_fh10.csv'
if False:
    F = ap_to_df(array_processing(st, win_len = 60, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 0.5, frqhigh = 10, stime = t1, etime = t2, prewhiten = False), fn)
else:
    F = pd.read_csv(fn)
t = (F.iloc[:,0] - np.floor(F.iloc[0,0]))*24-6
for i in range(1, 5):
    plt.subplot(2,2,i)
    if i in [1,2]:
        plt.semilogy(t, F.iloc[:,i], 'k.')
    else:
        plt.plot(t, F.iloc[:,i], 'k.')
plt.show()
###############################
## Burn day, higher freq (focus on helicopter)
fn = 'fire_TOP_win60_fl10_fh35.csv'
if False:
    F = ap_to_df(array_processing(st, win_len = 60, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 10, frqhigh = 35, stime = t1, etime = t2, prewhiten = False), fn)
else:
    F = pd.read_csv(fn)

plt.show()

### Burn day, low freq narrow 0.5-2 Hz (attempt no helicopter)
t1 = UTCDateTime('2023-10-06T00:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
st = st.slice(t1, t2)
st = st.select(station='TOP*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333)
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st, inv)

fn = 'fire_TOP_win60_fl0.5_fh2.csv'
if True:
    F = ap_to_df(array_processing(st, win_len = 60, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 0.5, frqhigh = 2, stime = t1, etime = t2, prewhiten = False, store = store), fn)
else:
    F = pd.read_csv(fn)


### Burn day, high freq narrow 25-35 Hz
t1 = UTCDateTime('2023-10-06T00:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
st = st.slice(t1, t2)
st = st.select(station='TOP*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333)
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st, inv)

fn = 'fire_TOP_win60_fl25_fh35.csv'
if False:
    F = ap_to_df(array_processing(st, win_len = 60, win_frac = 1, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 25, frqhigh = 35, stime = t1, etime = t2, prewhiten = False, store = store), fn)
else:
    F = pd.read_csv(fn)
t = (F.iloc[:,0] - np.floor(F.iloc[0,0]))*24-6
for i in range(1, 5):
    plt.subplot(2,2,i)
    if i in [1,2]:
        plt.semilogy(t, F.iloc[:,i], 'k.')
    else:
        plt.plot(t, F.iloc[:,i], 'k.')
plt.show()



### Burn day, low freq moderate 1-5 Hz
t1 = UTCDateTime('2023-10-06T00:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
st = st.slice(t1, t2)
st = st.select(station='TOP*')
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.333)
inv = obspy.read_inventory('RCEW_inventory.xml')
add_inv_coords(st, inv)

fn = 'fire_TOP_win60_fl1_fh5.csv'
if True:
    F = ap_to_df(array_processing(st, win_len = 60, win_frac = 0.5, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = 1, frqhigh = 5, stime = t1, etime = t2, prewhiten = False, store = store), fn)
else:
    F = pd.read_csv(fn)
plot_bf_result('fire_TOP_win60_fl1_fh5.csv', 'k.')




#########
plot_bf_result('fire_TOP_win60_fl1_fh5.csv', 'ko')
plot_bf_result('fire_TOP_win60_fl25_fh35.csv', 'r.')


#####################################
## basic spectrograms--no stacking
import obspy, riversound
tr = obspy.read('mseed_to_share/*JDNA2*')[0]
#tr = obspy.read('mseed_to_share/*JDSB2*')[0]
#tr = obspy.read('mseed_to_share/*TOP11*')[0]
t1 = obspy.UTCDateTime('2023-10-06T08:00:00')
t2 = obspy.UTCDateTime('2023-10-07T00:00:00')
tr.filter('highpass', freq=0.1)
tr.trim(t1, t2)
S = riversound.spectrum(tr, criterion_function = None, nfft = 2**12, runmed_radius_t = 1)
plt.figure()
cleanbf.image(np.log(S['specgram'].T), S['times'], np.log10(S['freqs']+1e-2))
#%% loop through multiple stations and frequency bands: first, read in the data
####################################
####################################
## loop through frequency bands
t1 = UTCDateTime('2023-10-06T12:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge().slice(t1, t2)
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.05, corners = 6)
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')
add_inv_coords(st, inv)
#%%
#station_list = ['JDNA', 'JDNB', 'JDSA', 'JDSB', 'QST', 'TOP']
station_list = ['VPT', 'CHKB']
fl_list = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 25]
fh_list = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 32]
#fl_list = [25]
#fh_list = [32]
win_len = 60
win_frac = 1
skip = []
for station in station_list:
    if station[:3] == 'VPT':
        t1 = obspy.UTCDateTime('2023-10-06T18:00:00')
        t2 = obspy.UTCDateTime('2023-10-06T22:50:00')
    elif station[:4] == 'CHKB':
        t1 = obspy.UTCDateTime('2023-10-06T15:00:00')
        t2 = UTCDateTime('2023-10-06T23:59:59')
    else:
        t1 = obspy.UTCDateTime('2023-10-06T12:00:00')
        t2 = UTCDateTime('2023-10-06T23:59:59')

    st_sta = st.select(station=f'{station}*')
    for (fl, fh) in zip(fl_list, fh_list):
        if [fl, fh, station] in skip:
            continue
        print([station, fl, fh])
        filename = f'beamform_results/{t1.isoformat()[5:-3]}_{t2.isoformat()[5:-3]}_{station}_fl{fl}_fh{fh}_winlen{win_len}_winfrac{win_frac}.csv'
        F = ap_to_df(array_processing(st_sta, win_len = win_len, win_frac = win_frac, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = fl, frqhigh = fh, stime = t1, etime = t2, prewhiten = False, store = store), filename)    
#%%
########
## plot results

for station in station_list:
    fn_head = glob.glob(f'beamform_results/10-06*{station}*csv')[0].split(station)[0]
    plt.figure()
    plot_bf_result_all(f'{fn_head}{station}_fl0.25_fh0.5_winlen60_winfrac1.csv', 'k.')
    plot_bf_result_all(f'{fn_head}{station}_fl0.5_fh1_winlen60_winfrac1.csv', 'r.')
    plot_bf_result_all(f'{fn_head}{station}_fl1_fh2_winlen60_winfrac1.csv', 'g.')
    plot_bf_result_all(f'{fn_head}{station}_fl2_fh4_winlen60_winfrac1.csv', 'b.')
    plot_bf_result_all(f'{fn_head}{station}_fl4_fh8_winlen60_winfrac1.csv', 'c.')
    plot_bf_result_all(f'{fn_head}{station}_fl8_fh16_winlen60_winfrac1.csv', 'm.')
    plot_bf_result_all(f'{fn_head}{station}_fl16_fh32_winlen60_winfrac1.csv', 'y.')
    plt.title(station)
#%%
#######################################
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge()
t1 = UTCDateTime('2023-10-06T19:12:00')

t2 = UTCDateTime('2023-10-06T19:16:00')
st = obspy.read(path).merge().slice(t1-30, t2+30)
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.05, corners = 6)
st.trim(t1, t2)
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')
add_inv_coords(st, inv)
station = 'TOP'

st = st.select(station=f'{station}*')
w = np.where(np.array(st.std()) < 60)[0]
st = obspy.Stream([st[i] for i in w])
win_length_sec = 10
sx = 0
sy = 3 # forward slowness

freqs, power, semblance = beam_stack(st, sx, sy, 5)

plt.figure()
plt.subplot(2,1,1)
plt.loglog(freqs, power)
plt.subplot(2,1,2)
plt.loglog(freqs, semblance)

##################
#beam stack spectrogram
t1 = UTCDateTime('2023-10-06T12:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
path = '/home/jake/2023-10-20_FireDataDownload/mseed_to_share/2023-10-06*'
st = obspy.read(path).merge().slice(t1-30, t2+30)
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.1, corners = 6)
st.trim(t1, t2)
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')
cleanbf.add_inv_coords(st, inv)

station = 'TOP'

st = st.select(station=f'{station}*')
w = np.where(np.array(st.std()) < 60)[0]
st = obspy.Stream([st[i] for i in w])

S = pd.read_csv('beamform_results/10-06T12:00_10-06T23:59_TOP_fl4_fh8_winlen60_winfrac1.csv')

sg = beam_stack_spectrogram(st, S.backazimuth+180, S.slowness, win_length_sec = 60, welch_ratio=6)

plt.subplot(3,1,2)
cleanbf.image(np.log10(sg['power']), ((S['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False)
plt.ylim(-1.0, 1.65)
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Spectrogram (beam-stack power)')
plt.ylabel('Frequency (Hz)')

plt.subplot(3,1,3)
cleanbf.image(np.log10(np.abs(sg['semblance'])), ((S['t'] % 1) * 24)-6, np.log10(sg['freqs'][0,:]), crosshairs = False, zmin = -1.5, zmax = 0)
plt.ylim(-1.0, 1.65)
plt.yticks(np.log10([0.3, 1, 3, 10, 30]), [0.3, 1, 3, 10, 30])
plt.title('Spectrogram (beam-stack semblance)')
plt.xlabel('Local Time (hours)')
plt.ylabel('Frequency (Hz)')

plt.tight_layout()
for x in [16+41/60, 15+12/60, 14+23/60, 13+12/60]:
    plt.axvline(x)
