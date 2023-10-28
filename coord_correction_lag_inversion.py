#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 19:01:13 2023

@author: jake
"""

import obspy, gemlog
import numpy as np
import matplotlib.pyplot as plt
from gemlog.xcorr import xcorr_function, get_coordinates, apply_function_windows
from obspy.signal.cross_correlation import correlate, xcorr_max
from scipy.optimize import minimize

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
    
    
    
def xcorr_function(st, args):
    maxshift_seconds = args.get('maxshift_seconds')
    if maxshift_seconds is None:
        maxshift_seconds = 1

    quiet = args.get('quiet')
    if quiet is None:
        quiet = False
    #if not quiet:
    #    print(st)
    dt = st[0].stats.delta
    st.detrend('linear')
    st.taper(0.05, 'hann') # Hann window, default taper for SAC
    output_dict = {'mean_coef':0}
    consistency = 0
    ## simplest loop: go through sensors in order, O(N). Not robust if a single sensor misbehaves.
    for i in range(0, len(st)):
        j = (i+1) % len(st)
        tr1 = st[i]
        tr2 = st[j]
        ID1 = f'{tr1.stats.network}.{tr1.stats.station}.{tr1.stats.location}'
        ID2 = f'{tr2.stats.network}.{tr2.stats.station}.{tr2.stats.location}'
        output_dict[f'rms_{ID1}'] = tr1.std()
        pair_name = f'{ID1}_{ID2}'
        shift = correlate(tr1.data, tr2.data, int(np.round(maxshift_seconds / dt)))
        output_dict[f'lag_{pair_name}'] = shift * dt
        #output_dict[f'r_{pair_name}'] = value
        consistency += shift/len(st)
        #output_dict['mean_coef'] += value/len(st)
    ## consistency: allow up to one sample error per cross-correlation
    output_dict['consistency'] = np.abs(consistency) <= len(st)
    return output_dict

def scaled_sum_square(dx, dy, t_resid):
    return np.sum((dx/0.003)**2 + (dy/0.003)**2 + (t_resid/0.01)**2)

#%% read stream
st = obspy.read('/home/jake/Dropbox/Aftershocks/all_PARK_daily_mseed_renamed/PARK.XP.*.HDF.2020.163')
t1 = obspy.UTCDateTime('2020-06-11T04:00:00.000000Z')
t2 = t1 + 3600 * 1
#%%
st = obspy.read('/home/jake/Dropbox/Aftershocks/all_PARK_daily_mseed_renamed/PARK.XP.*.HDF.2020.146')
t1 = obspy.UTCDateTime('2020-05-25T09:00:00.000000Z')
t2 = t1 + 3600 * 0.1
#%%
st.filter('bandpass', freqmin = 2, freqmax = 8, corners = 4)
st = obspy.Stream([tr for tr in st if tr.stats.location != '13'])

#%% trim to quiet period

st.trim(t1, t2)
N=len(st)

#%% find lags (old)
d = xcorr_function(st, {})
lags = [d[k] for k in d.keys() if k.find('lag')==0]
print(np.sum(lags))
#%% find lags (window averaging)
d = apply_function_windows(st, xcorr_function, win_length_sec = 5)
#%%
for k in d.keys(): 
    if k.find('lag')==0:
        print(k)
        plt.plot(d[k].sum(0))
#%%
#lags = np.array([xcorr_max(np.median(d[k],0))[0] for k in d.keys() if k.find('lag')==0]) * 0.01
lags = np.array([xcorr_max(np.sum(d[k],0))[0] for k in d.keys() if k.find('lag')==0]) * 0.01

#%%

#%%
inv = obspy.read_inventory('/home/jake/Dropbox/Aftershocks/XP_PARK_inventory.xml')
add_inv_coords(st, inv)
coords = get_coordinates(st)
x = coords.x
y = coords.y
plt.plot(x, y, 'k.')
for i in range(N):
    plt.text(x[i], y[i], coords.location[i])
#%%
#%% define G matrix

def objective(p):
    dx = p[:(len(p)//2)]
    dy = p[(len(p)//2):]
    
    G = np.zeros((N,2))

    for i in range(0, N):
        G[i,0] = (x+dx)[(i+1)%N] - (x+dx)[i]
        G[i,1] = (y+dy)[(i+1)%N] - (y+dy)[i]

    H = np.dot(np.linalg.inv(np.dot(G.T, G)), G.T)

    s = np.dot(H, lags)
    #s = 3 * np.array([np.sin(232*np.pi/180), np.cos(232*np.pi/180)]) #LFF
    t_resid = np.dot(G, s) - lags
    return scaled_sum_square(dx, dy, t_resid)
#%% do a search for the coordinate perturbations that minimize the objective function
print(objective(np.zeros(N*2)))
a = minimize(objective, np.zeros(N*2))
print(objective(a.x))
dx = a.x[:N]
dy = a.x[N:]
plt.plot(x,y,'k.')
plt.plot(x+dx, y+dy, 'go')