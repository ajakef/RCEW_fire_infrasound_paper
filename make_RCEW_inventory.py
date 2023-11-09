import pandas as pd
import gemlog, obspy
import numpy as np
dtype = {key:str for key in ['station', 'location', 'network', 'SN', 'description','data issues']}
dtype.update({key:float for key in ['lat', 'lon', 'elevation']})
coords = pd.read_csv('coordinates/RCEW_all_coordinates.csv', dtype = dtype, keep_default_na=False)
response = gemlog.get_gem_response(gain = 'high')
inv = gemlog.make_gem_inventory('coordinates/RCEW_station_info.csv', coords, response) # station_info file must only have station/network/location/SN with optional elevation

#%% create a combined infraBSU/Datacube response file
infraBSU_response = obspy.read_inventory('https://service.iris.edu/irisws/nrl/1/combine?instconfig=sensor_JeffreyBJohnson_infraBSUv2_LP12_SG0.000047_STairPressure&format=resp')[0][0][0].response
datacube_response = obspy.read_inventory('https://service.iris.edu/irisws/nrl/1/combine?instconfig=datalogger_DiGOSOmnirecs_DataCube_PG64_FR400&format=resp')[0][0][0].response
datacube_response.response_stages[0].input_units = 'V'
datacube_response.response_stages[0].input_units_description = ''

## code copied and modified from obspy's nrl.get_response() to merge sensor and logger responses
response = datacube_response
sensor_resp = infraBSU_response

## incremment stage numbers for datacube stages, because the infraBSU will be stage 0
for i in range(3): response.response_stages[i].stage_sequence_number += 1

## insert the infraBSU as the first stage
sensor_stage0 = sensor_resp.response_stages[0]
response.response_stages.insert(0, sensor_stage0)

## change the input units of the response
response.instrument_sensitivity.input_units = sensor_stage0.input_units
response.instrument_sensitivity.input_units_description = sensor_stage0.input_units_description

## manually calculate the overall response--should be about 12000 counts/Pa
filter_gain_100 = np.abs(1-1/(1+100j*12.0)) # 12-second corner period
response.instrument_sensitivity.value = response.response_stages[0].stage_gain * response.response_stages[1].stage_gain * response.response_stages[2].stage_gain * filter_gain_100


#%% Correct the datacube stations' sensor and response info
for station in inv[0]:
    if station.code[:3] == 'VPT':
        for channel in station:
            channel.code = 'CDF'
            channel.response = response
            channel.sensor.description = 'infraBSU v2, 12s, 47 uV/Pa'
            channel.sensor.model = 'infraBSU v2'
            channel.sensor.serial_number = 'None'

#%% save the files
inv.write('coordinates/RCEW_inventory.xml', 'stationXML')
inv.write('coordinates/RCEW_inventory.kml', 'kml')
#%% just the OSIRIS-REx 
inv = inv.select(station = '[JT]*') # JD stations and TOP
inv = inv.remove(station='JDNA*')
#%%
inv.write('coordinates/RCEW_OSIRIS-REx_inventory.xml', 'stationXML')
inv.write('coordinates/RCEW_OSIRIS-REx_inventory.kml', 'kml')