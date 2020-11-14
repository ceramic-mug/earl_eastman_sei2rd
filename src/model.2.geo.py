#########################################
# title: model.2.geo                    #
# author: Joshua Eastman                #
# description: Processing GIS Zip-code  #
#        data, subway data, and pop     #
#########################################


#########################################
#
# Zip code data extraction and prep    
#
#########################################

# import gis packages
import os
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon

# get dirname of the dir holding this file
dirname = os.path.dirname(os.path.abspath(__file__))

# get ZIPs shapes
fZip = fiona.open(dirname+'/../gis/ZIP_CODE_040114/ZIP_CODE_040114.shp')
# subway stops
fSubway = fiona.open(dirname+'/../gis/stops_nyc_subway_may2020/stops_nyc_subway_may2020.shp')

# parse out the stops
stops = [stop for stop in fSubway]

# begin writing to a new csv to hold the output
f = open(dirname+'/../in/zipSubwayData.out',"w+")
f.write(','.join(['zipcode','pop','borough','stop_id','trains']))

# for every Zip Code
for i, zip in enumerate(fZip):

    print('scanning '+zip['properties']['ZIPCODE'])

    # for every stop
    for j, s in enumerate(stops):

        # get the stop point
        stop = shape(s['geometry'])

        # if the stop is within the Zip Code
        if stop.within(shape(zip['geometry'])):

            print('   found stop '+s['properties']['stop_id'])

            try:
                f.write('\n'+','.join([zip['properties']['ZIPCODE'],'{:.0f}'.format(zip['properties']['POPULATION']),zip['properties']['COUNTY'],s['properties']['stop_id'],s['properties']['trains']]))
            except:
                f.write('\n'+','.join([zip['properties']['ZIPCODE'],'{:.0f}'.format(zip['properties']['POPULATION']),zip['properties']['COUNTY'],s['properties']['stop_id'],'CLOSED']))
f.close()

# creating subway enter/exit datafiles

import pandas as pd

turnstyle_counts = pd.read_csv('../gis/nyc-transit-data-turnstile_daily_counts_2020_-_2020-10-05-14-10-52.csv')