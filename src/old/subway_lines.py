# import gis packages
import fiona
from shapely.geometry import shape,mapping, Point, Polygon, MultiPolygon
# get the borough shapefile
borough = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /boroughs/nybb.shp')
# get the subway line shapegile
subways = fiona.open('/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/gis /routes_nyc_subway_may2019/routes_nyc_subway_may2019.shp')

# this code altered from https://gis.stackexchange.com/questions/208546/check-if-a-point-falls-within-a-multipolygon-with-python
lines = ([line for line in subways])

f = open("/Users/joshua/Documents/School/Princeton/Sophomore Classes/Spring 2020/CEE302/term_project/SEI2RD-dt-model/dat/borough_subways.csv","w+")
f.write(','.join(['borough','route_id','color']))

for i, bo in enumerate(borough):
    for j, line in enumerate(lines):
        ln = shape(line['geometry'])
        if ln.within(shape(bo['geometry'])) or ln.intersects(shape(bo['geometry'])):
            colors.append(line['properties']['color'])
            f.write('\n'+','.join([bo['properties']['BoroName'],line['properties']['route_id'],line['properties']['color']]))

f.close()
