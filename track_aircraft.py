import matplotlib.pyplot as plt
import xmltodict, obspy, glob
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth

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


fn = sorted(glob.glob('helicopter_data/N684H-track*kml'))
fn = fn[5:] + fn[:5]
heli_track = read_adsb_kml(fn)


##################################

## find a mean location for every array
inv = obspy.read_inventory('coordinates/RCEW_inventory.xml')[0]
station_list = ['TOP', 'JDNA', 'JDNB', 'JDSA', 'JDSB', 'QST', 'CHKB', 'VPT']
coords_dict = {'station':[], 'lat':[], 'lon':[], 'z':[]}
for station in station_list:
    coords_dict['station'].append(station)
    coords_dict['lat'].append(np.mean([s.latitude for s in inv if s.code[:len(station)] == station]))
    coords_dict['lon'].append(np.mean([s.longitude for s in inv if s.code[:len(station)] == station]))
    coords_dict['z'].append(np.mean([s.elevation for s in inv if s.code[:len(station)] == station]))

######################################
for station in station_list:
    #station = 'TOP'
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
heli_track.to_csv('helicopter_data/heli_track.csv', index = False)





#####################################################
## code to read historical data files (readsb-hist and individual planes).
## this might not end up being useful. maybe best to download KMLs from the replay map page.
import json, glob, gzip
def download_many_snapshots(url, t1, t2):
    t = t1
    while t < t2:
        ts = t.isoformat()
        filename = f'{ts[11:13]}{ts[14:16]}{ts[17:19]}Z.json.gz'
        urlretrieve(url + filename)

def read_many_snapshots(t1, t2, lat_min = 43, lat_max = 43.3, lon_min = -117, lon_max = -116.6):
    t1 = obspy.UTCDateTime(t1)
    t2 = obspy.UTCDateTime(t2)
    t1 -= (t1.second % 5) + (t1.microsecond/1.0e6)
    hexcodes = set()
    t = t1
    while t < t2:
        print(t)
        ts = t.isoformat()
        url = f'https://samples.adsbexchange.com/readsb-hist/{ts[:4]}/{ts[5:7]}/{ts[8:10]}'
        filename = f'{ts[11:13]}{ts[14:16]}{ts[17:19]}Z.json.gz'
        with gzip.open(f'{url}/{filename}', 'rb') as f:
            data = json.load(f)
        for i, d in enumerate(data['aircraft']):
            good = False
            try:
                if (d['lat'] > lat_min) and (d['lat'] < lat_max) \
                   and (d['lon'] > lon_min) and (d['lon'] < lon_max):
                    good = True
                    print(1)
            except:
                try:
                    if (d['mlat'] > lat_min) and (d['mlat'] < lat_max) \
                       and (d['mlon'] > lon_min) and (d['mlon'] < lon_max):
                        good = True
                        print(2)
                except:
                    try:
                        d = d['lastPosition']
                        if (d['lat'] > lat_min) and (d['lat'] < lat_max) \
                           and (d['lon'] > lon_min) and (d['lon'] < lon_max):
                            good = True
                            print(3)
                    except:
                        pass
        if good:
            print(data['aircraft'][i]['hex'])
            hexcodes.add(data['aircraft'][i]['hex'])
        t += 5
    return hexcodes
    
hexcodes = set()
all_json_filenames = sorted(glob.glob('snapshot_json/*json'))
for json_filename in all_json_filenames:
    #json_filename = '2023-10-01-000000Z.json'
    #with gzip.open(gz_filename, 'rb') as f:
    with open(json_filename) as f:
        data = json.load(f)

    lat_min = 43
    lat_max = 43.3
    lon_min = -116.9
    lon_max = -116.7

    indices = []
    ## find all aircraft passing through region
    for i, d in enumerate(data['aircraft']):
    good = False
        try:
            if (d['lat'] > lat_min) and (d['lat'] < lat_max) \
               and (d['lon'] > lon_min) and (d['lon'] < lon_max):
                good = True
                print(1)
        except:
            try:
                if (d['mlat'] > lat_min) and (d['mlat'] < lat_max) \
                   and (d['mlon'] > lon_min) and (d['mlon'] < lon_max):
                    good = True
                    print(2)
            except:
                try:
                    d = d['lastPosition']
                    if (d['lat'] > lat_min) and (d['lat'] < lat_max) \
                       and (d['lon'] > lon_min) and (d['lon'] < lon_max):
                        good = True
                        print(3)
                except:
                    pass
        if good:
            print(data['aircraft'][i]['hex'])
            hexcodes.add(data['aircraft'][i]['hex'])

## track one aircraft            
json_filename = 'trace_full_a1311e.json' # plane passed RCEW at 00 UTC + 2000 sec, then stopped at Boise, Portland, etc
with open(json_filename) as json_file:
    data = json.load(json_file)

trace = data['trace']

t = [d[0] for d in trace]
lat = [d[1] for d in trace]
lon = [d[2] for d in trace]
z = [d[3] for d in trace]

plt.plot(lon, lat)
