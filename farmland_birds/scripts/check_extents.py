#! /usr/bin/env pyhton3

# Run pip3 install -r requirements.txt first
import glob
import logging
import pprint
import rasterio
import sys


""" Check that the extent of all translated files matches.
"""

# Helper function to check equality of list items
def add_item(dic, key, value):
    if key in dic.keys():
        dic[key] = dic[key] + [value]
    else:
        dic[key] = [value]
    return dic

# Get the inputs and outputs
input_dir = snakemake.input[0] if 'snakemake' in sys.modules and hasattr(snakemake, "input") else "data/birds"
log_file = snakemake.log[0] if 'snakemake' in sys.modules and hasattr(snakemake, "log") else "log/extents.log"

logger = logging.getLogger("check_extent")
logger.setLevel(logging.DEBUG)
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")

fileHandler = logging.FileHandler(log_file)
fileHandler.setFormatter(logFormatter)
logger.addHandler(fileHandler)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

bird_rasters = glob.glob(input_dir + "/*.tif")

if len(bird_rasters) == 0:
    logger.error("No GeoTIFFs found in the input folder")
    sys.exit(0)

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
    logger.warning("Multiple minx found")
    pprint.pprint(minx)
else:
    logger.info("Only one minx found: {0}".format(list(minx)[0]))
if len(miny.keys()) > 1:
    logger.warning("Multiple miny found")
    pprint.pprint(miny)
else:
    logger.info("Only one miny found: {0}".format(list(miny)[0]))
if len(maxx.keys()) > 1:
    logger.warning("Multiple maxx found")
    pprint.pprint(maxx)
else:
    logger.info("Only one maxx found: {0}".format(list(maxx)[0]))
if len(maxy.keys()) > 1:
    logger.warning("Multiple maxy found")
    pprint.pprint(maxy)
else:
    logger.info("Only one maxy found: {0}".format(list(maxy)[0]))
