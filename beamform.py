from obspy import UTCDateTime
import obspy, glob
from obspy.signal.array_analysis import array_processing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
base_dir = '..'
sys.path.append(f'{base_dir}/code')

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
            
#%% beamforming for all stations on 10/6
# loop through multiple stations and frequency bands: first, read in the data
####################################
####################################
## loop through frequency bands
t1 = UTCDateTime('2023-10-06T12:00:00')
t2 = UTCDateTime('2023-10-06T23:59:59')
infrasound_files = f'{base_dir}/data/infrasound/2023-10-06*'
st = obspy.read(infrasound_files).merge().slice(t1, t2)
st = obspy.Stream([tr for tr in st if type(tr.data) is np.ndarray])
st.filter('highpass', freq = 0.05, corners = 6)
inv = obspy.read_inventory(f'{base_dir}/data/coordinates/RCEW_inventory.xml')
add_inv_coords(st, inv)
#%%
station_list = ['JNA', 'JNB', 'JSA', 'JSB', 'QST', 'TOP']
#station_list = ['VPT', 'CHKB'] # fire infrasound not observed at VPT & CHKB
fl_list = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 25]
fh_list = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 32]
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
        filename = f'{base_dir}/beamform_results/{t1.isoformat()[5:-3]}_{t2.isoformat()[5:-3]}_{station}_fl{fl}_fh{fh}_winlen{win_len}_winfrac{win_frac}.csv'
        F = ap_to_df(array_processing(st_sta, win_len = win_len, win_frac = win_frac, sll_x = -4, slm_x = 4, sll_y = -4, slm_y = 4, sl_s = 0.1, semb_thres = 0, vel_thres = 0, frqlow = fl, frqhigh = fh, stime = t1, etime = t2, prewhiten = False, store = store), filename)    
#%%
########
## plot results--this is just for quick visualization and does not make figures for the paper

for station in station_list:
    fn_head = glob.glob(f'../beamform_results/10-06*{station}*csv')[0].split(station)[0]
    plt.figure()
    plot_bf_result_all(f'{fn_head}{station}_fl0.125_fh0.25_winlen60_winfrac1.csv', 'y+')
    plot_bf_result_all(f'{fn_head}{station}_fl0.25_fh0.5_winlen60_winfrac1.csv', 'k.')
    plot_bf_result_all(f'{fn_head}{station}_fl0.5_fh1_winlen60_winfrac1.csv', 'r.')
    plot_bf_result_all(f'{fn_head}{station}_fl1_fh2_winlen60_winfrac1.csv', 'g.')
    plot_bf_result_all(f'{fn_head}{station}_fl2_fh4_winlen60_winfrac1.csv', 'b.')
    plot_bf_result_all(f'{fn_head}{station}_fl4_fh8_winlen60_winfrac1.csv', 'c.')
    plot_bf_result_all(f'{fn_head}{station}_fl8_fh16_winlen60_winfrac1.csv', 'm.')
    plot_bf_result_all(f'{fn_head}{station}_fl16_fh32_winlen60_winfrac1.csv', 'y.')
    plt.title(station)
