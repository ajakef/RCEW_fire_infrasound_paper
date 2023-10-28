import pandas as pd
import gemlog
dtype = {key:str for key in ['station', 'location', 'network', 'SN', 'description','data issues']}
dtype.update({key:float for key in ['lat', 'lon', 'elevation']})
coords = pd.read_csv('RCEW_all_coordinates.csv', dtype = dtype, keep_default_na=False)
response = gemlog.get_gem_response(gain = 'high')
inv = gemlog.make_gem_inventory('RCEW_station_info.csv', coords, response) # station_info file must only have station/network/location/SN with optional elevation
inv.write('RCEW_inventory.xml', 'stationXML')
inv.write('RCEW_inventory.kml', 'kml')
