# import gis packages
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon

# get NTA shapes
ntas = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /nynta_19d/nynta.shp')
# subway stops
subway = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /stops_nyc_subway_may2019/stops_nyc_subway_may2019.shp')

# parse out the stops
stops = [stop for stop in subway]

# begin writing to a new csv to hold the output
f = open("/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/nta_subway_stops.csv","w+")
f.write(','.join(['nta_code','stop_id','train']))

# for every NTA
for i, nta in enumerate(ntas):

    # for every stop
    for j, s in enumerate(stops):

        # get the stop point
        stop = shape(s['geometry'])

        # if the stop is within the NTA
        if stop.within(shape(nta['geometry'])):

            f.write('\n'+','.join([nta['properties']['NTACode'],s['properties']['stop_id'],s['properties']['trains']]))

f.close()
