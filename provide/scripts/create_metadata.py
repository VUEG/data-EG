#! /usr/bin/env pyhton3

# ONLY TO BE RUN FROM SNAKEMAKE

# Run pip3 install -r requirements.txt first
import datapackage
import json
import os
import sys
from collections import OrderedDict
from egpackager.datamanager import DataManager

# Create a DataManager instance to handle the data
dm = DataManager()

gdocs_uri = snakemake.params.gspread_uri
dm.add_datasource(type='gspread',
                  uri=gdocs_uri,
                  credentials_file=snakemake.params.gpsread_credentials,
                  spreadsheet_name=snakemake.params.gpsread_spreadsheet_name,
                  sheet=snakemake.params.gpsread_worksheet_name)

for i in range(0, len(snakemake.input)):
    input_raster = snakemake.input[i]
    # Get metadata directly from the raster
    dm.add_datasource(type='raster', uri=input_raster)
    # 'local_resource_name' is used to match the data from Google spreadsheet
    local_resource_name = os.path.basename(input_raster)
    # Find the dataset named based on 'local_resource_name'
    name = dm.find_name(local_resource_name)
    # If name is not found, report but don't stop
    if name is None:
        print("WARNING: Can't find dataset name for resource {0}".format(local_resource_name))
        next
    # Create a new DataPackage object. NOTE: 'name' is not the same as
    # 'resource_name' name (which is the file name).
    dp = datapackage.DataPackage(name=name)
    # Bump version to 1.0
    dp.bump_major_version()
    # Title
    dp.title = dm.get_metadata_value(name, 'original_name')
    # Description
    dp.description = dm.get_metadata_value(name, 'description')
    # Sources
    dp.sources = [{
        "name":"Environmental Geography group",
        "web":"https://www.falw.vu.nl/en/research/earth-sciences/earth-and-climate-cluster/research/environmental-geography/"}]
    # Publishers
    dp.publishers = [{
        "name":dm.get_metadata_value(name, "point_of_origin"),
        "email":dm.get_metadata_value(name, "point_of_origin_email")
    }]
    # Maintainers
    dp.maintainers = [{
        "name":dm.get_metadata_value(name, "added_by"),
        "email":dm.get_metadata_value(name, "added_by_email")
    }]
    # Resources
    dp.resources = [{
        "name":local_resource_name,
        "path":"{0}".format(local_resource_name),
        "raster_metadata":dm.resource_metadata
    }]
    # Make sure the order is correct
    attribute_order = ["name", "version", "title", "description", "sources",
                       "publishers", "maintainers", "resources"]
    dp_data = dp.as_dict()
    dp_ordered_dict = OrderedDict((k, dp_data[k]) for k in attribute_order)

    with open(snakemake.output[i], 'w') as outfile:
        outfile.write('{}\n'.format(json.dumps(dp_ordered_dict, indent=4)))
