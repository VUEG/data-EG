#! /usr/bin/env pyhton3

# Run pip3 install -r requirements.txt first
import glob
import pprint
import rasterio


""" Check that the extent of all translated files matches.
"""

# Helper function to check equality of list items
def add_item(dic, key, value):
    if key in dic.keys():
        dic[key] = dic[key] + [value]
    else:
        dic[key] = [value]
    return dic

input_dir = "data/birds/"
bird_rasters = glob.glob(input_dir + "*.tif")

# Make lists to holds the bound values
minx = {}
miny = {}
maxx = {}
maxy = {}

for raster in bird_rasters:
    with rasterio.drivers():
        with rasterio.open(raster) as src:
            minx = add_item(minx, src.bounds[0], raster)
            miny = add_item(miny, src.bounds[1], raster)
            maxx = add_item(maxx, src.bounds[2], raster)
            maxy = add_item(maxy, src.bounds[3], raster)

print("{0} rasters found in {1}".format(len(bird_rasters), input_dir))

if len(minx.keys()) > 1:
    print("Multiple minx found")
    pprint.pprint(minx)
else:
    print("Only one minx found: {0}".format(list(minx)[0]))
if len(miny.keys()) > 1:
    print("Multiple miny found")
    pprint.pprint(miny)
else:
    print("Only one miny found: {0}".format(list(miny)[0]))
if len(maxx.keys()) > 1:
    print("Multiple maxx found")
    pprint.pprint(maxx)
else:
    print("Only one maxx found: {0}".format(list(maxx)[0]))
if len(maxy.keys()) > 1:
    print("Multiple maxy found")
    pprint.pprint(maxy)
else:
    print("Only one maxy found: {0}".format(list(maxy)[0]))
