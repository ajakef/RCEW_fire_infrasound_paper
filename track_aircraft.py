#import matplotlib.pyplot as plt
import xmltodict, obspy, glob
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth
#import json, glob, gzip
import pandas as pd
base_dir = '..'
#Replay map! This is really intuitive and easy. Click the aircraft and download its tracks. https://globe.adsbexchange.com/?replay=2023-10-06-20:19&lat=42.921&lon=-116.632&zoom=10.4
# example: the burn-day helicopter is reg N684H, A912C2, type code B407
# 20:23Z: went north of JD to refuel. Not seen N of Hill Rd (i.e., always south of quonset).

def read_adsb_kml(kml_files):
    times = []
    lat = []
    lon = []
    z = []
    for kml_file in kml_files:
        with open(kml_file, 'r') as file:
            xml_string=file.read()
        xml_dict=xmltodict.parse(xml_string)
        data = xml_dict['kml']['Folder']['Folder']['Placemark']['gx:Track']
        times += [obspy.UTCDateTime(t) for t in data['when']]
        for d in data['gx:coord']:
            s = d.split(' ')
            lat.append(float(s[1]))
            lon.append(float(s[0]))
            z.append(float(s[2]))
    return {'t':times, 'lat':lat, 'lon':lon, 'z':z}


fn = sorted(glob.glob(f'{base_dir}/data/helicopter_data/N684H-track*kml'))
fn = fn[5:] + fn[:5]
heli_track = read_adsb_kml(fn)


##################################

## find a mean location for every array
inv = obspy.read_inventory(f'{base_dir}/data/coordinates/RCEW_inventory.xml')[0]
station_list = ['TOP', 'JNA', 'JNB', 'JSA', 'JSB', 'QST', 'CHKB', 'VPT']
coords_dict = {'station':[], 'lat':[], 'lon':[], 'z':[]}
for station in station_list:
    coords_dict['station'].append(station)
    coords_dict['lat'].append(np.mean([s.latitude for s in inv if s.code[:len(station)] == station]))
    coords_dict['lon'].append(np.mean([s.longitude for s in inv if s.code[:len(station)] == station]))
    coords_dict['z'].append(np.mean([s.elevation for s in inv if s.code[:len(station)] == station]))

for key in coords_dict.keys():
    coords_dict[key] = np.array(coords_dict[key])
######################################
for station in station_list:
    lat_sta = coords_dict['lat'][coords_dict['station'] == station]
    lon_sta = coords_dict['lon'][coords_dict['station'] == station]
    z_sta = coords_dict['z'][coords_dict['station'] == station]
    az = []
    dist = []
    inc = []
    for i in range(len(heli_track['lat'])):
        r, th1, th2 = gps2dist_azimuth(lat_sta, lon_sta, heli_track['lat'][i], heli_track['lon'][i])
        az.append(th1)
        dist.append(r)
        inc.append(np.arctan2(r, heli_track['z'][i] - z_sta) * 180/np.pi)
    heli_track[f'{station}_az'] = az
    heli_track[f'{station}_dist'] = dist
    heli_track[f'{station}_inc'] = inc

heli_track = pd.DataFrame(heli_track)
heli_track.to_csv(f'{base_dir}/data/helicopter_data/helicopter_station_azimuths.csv', index = False)


