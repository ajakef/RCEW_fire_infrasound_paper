#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 21:32:48 2023

@author: jake
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import obspy, datetime
#%% plot surface observations

for station in ['124']:#, 'jdt3b', '124']: #124 is above JD on road, 500 m from middle of TOP
    filename = f'met_data/RCEW_met_data/{station}/q{station}met20231011.dat'
    names = pd.read_csv(filename, skiprows=1, nrows=1).iloc[:1,:].keys()
    df = pd.read_csv(filename, skiprows = 4, names = names)
    df['time'] = [datetime.datetime.fromisoformat(t) for t in df.TIMESTAMP]

#%%
    w = (df.time > datetime.datetime.fromisoformat('2023-10-05')) & \
        (df.time < datetime.datetime.fromisoformat('2023-10-08'))
    df_trim = df.loc[w, :]

    plt.figure()
    plt.subplot(4,1,1)
    plt.plot(df_trim.time, df_trim.wnd3sa)
    plt.subplot(4,1,2)
    plt.plot(df_trim.time, df_trim.wnd3d)
    plt.subplot(4,1,3)
    plt.plot(df_trim.time, df_trim.tmpj)
    plt.subplot(4,1,4)
    plt.plot(df_trim.time, df_trim.hum3_Avg)
#%%    
    plt.figure()
    plt.subplot(3,1,1)
    plt.plot(df_trim.time, df_trim.wnd3sa * np.sin(df_trim.wnd3d * np.pi/180))
    plt.plot(df_trim.time, df_trim.wnd3sa * np.cos(df_trim.wnd3d * np.pi/180))
    plt.legend(['East-West', 'North-South'])
    plt.axhline(0)
    plt.subplot(3,1,2)
    plt.plot(df_trim.time, df_trim.tmpj)
    plt.subplot(3,1,3)
    plt.plot(df_trim.time, df_trim.hum3_Avg)

#%% read BOI radiosonde data
rs = pd.read_csv(f'met_data/radiosonde/BOI_20231006T1200_radiosonde_data.txt', names = ['x', 'pressure', 'z', 'temp', 'dew_pt', 'wind_dir', 'wind_speed'], skiprows = 4, sep = '\s+', na_values=99999).iloc[:,1:]

rs['dew_pt'] = rs['dew_pt'] * 0.1
rs['temp'] = rs['temp'] * 0.1
rs['wind_speed'] = rs.wind_speed * 0.1
rs['wx'] = rs['wind_speed'] * np.sin(rs['wind_dir'] * np.pi/180)
rs['wy'] = rs['wind_speed'] * np.cos(rs['wind_dir'] * np.pi/180)
plt.figure()
plt.subplot(1,3,1)
plt.plot(rs['wx'], rs['z'], '.')
plt.plot(rs['wy'], rs['z'], '.')
plt.plot(rs['wind_speed'], rs['z'], 'k.')
plt.axvline(0)
plt.subplot(1,3,2)
plt.plot(rs['wind_dir'], rs['z'], 'k.')
plt.subplot(1,3,3)
plt.plot(rs['temp'], rs['z'], '.')
plt.plot(rs['dew_pt'], rs['z'], '.')


#%%
    ## not sure what the difference is between wnd3s, wnd3sa, and wnd3sr
#plt.plot(df.wnd3sa, 'r.')
#plt.plot(df.wnd3sr, 'b.')
#plt.plot(df.wnd3dr, 'r.') # this is direction



print([df.wnd3s.mean(), df.wnd3s.std()])
print([df.wnd3sa.mean(), df.wnd3sa.std()])
print([df.wnd3sr.mean(), df.wnd3sr.std()])
