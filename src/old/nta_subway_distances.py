# Compute the distance from each NTA to its nearest subway stop
# output csv with nta, stop, colors of subway lines accessible by stop, and distance in ft to the stop

# import gis packages
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon
import numpy as np
import re

get_lines = re.compile('(\w*)\s*(\w*)\s*(\w*)\s*(\w*)')

# NTA shapes
ntas = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /nynta_19d/nynta.shp')

# Subway stops
subway = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /stops_nyc_subway_may2019/stops_nyc_subway_may2019.shp')

# subway line for matching colors
lines = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /routes_nyc_subway_may2019/routes_nyc_subway_may2019.shp')

line_colors = [(line['properties']['route_shor'],line['properties']['color']) for line in lines]

#[nta]: (stop ID, [colors], distance)
dict = {}

# stop: [colors]
stop_colors = {}

# assigning colors to stops
for stop in subway:

    trains_served = [i for i in get_lines.match(stop['properties']['trains']).groups() if i]

    for i in trains_served:

        for line in lines:

            if i == line['properties']['route_shor'] or i == line['properties']['route_id']:

                color = line['properties']['color']

                if stop['properties']['stop_id'] not in stop_colors.keys():
                    stop_colors[stop['properties']['stop_id']] = [color]
                else:
                    stop_colors[stop['properties']['stop_id']].append(color)

# begin writing to a new csv to hold the output
f = open("/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/nta_subway_stop_distances_to_closest.csv","w+")
f.write(','.join(['nta_code','stop_id','color1','color2','color3','color4','distance(ft)']))

# for every NTA
for i, nta in enumerate(ntas):

    n = shape(nta['geometry']).centroid

    print('NTA: '+nta['properties']['NTACode'])

    # distance,id,train
    shortest = [np.inf,'']

    # for every stop
    for j, stop in enumerate(subway):

        s = shape(stop['geometry'])

        d = n.distance(s)

        if d < shortest[0]:
            shortest[0]=d
            shortest[1]=stop['properties']['stop_id']

    print(' stop '+ str(shortest[1]) + ' closest')

    empties = ','*(4-len(stop_colors[shortest[1]]))

    f.write('\n'+','.join([str(nta['properties']['NTACode']),str(shortest[1])])+','+','.join(stop_colors[shortest[1]])+empties+','+str(shortest[0]))

f.close()
