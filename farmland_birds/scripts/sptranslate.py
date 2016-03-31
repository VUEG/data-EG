#! /usr/bin/env pyhton3

# Run pip3 install -r requirements.txt first

import logging
import os
import subprocess
import sys
import re
from pandas import read_csv


# Helper function to deal with whitespaces in the file name
def quote_item(item):
    return "\"" + item + "\""



# Get the inputs and outputs
input_dir = snakemake.input['input_dir'] if 'snakemake' in sys.modules and hasattr(snakemake, "input") else "data/birds/output"
sp_data_file = snakemake.input['sp_data'] if 'snakemake' in sys.modules and hasattr(snakemake, "input") else "data/birds/farmland_birds_sp.csv"
output_dir = snakemake.output[0] if 'snakemake' in sys.modules and hasattr(snakemake, "output") else "data/birds"
log_file = snakemake.log[0] if 'snakemake' in sys.modules and hasattr(snakemake, "output") else "log/sptranslate.log"

logger = logging.getLogger("")
logger.setLevel(logging.DEBUG)
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")

fileHandler = logging.FileHandler(log_file)
fileHandler.setFormatter(logFormatter)
logger.addHandler(fileHandler)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Filter the AIGs (these are dirs)
filtered_files = [f for f in os.listdir(input_dir) if re.match(r'^(b){1}.+(2000)$', f)]

# Load the additional CSV data
sp_data =  read_csv(sp_data_file)

n_files = len(filtered_files)
counter = 1

for aig in filtered_files:
    # Figure out the species. First, get rid of '2000' and 'a1' in the AIG name
    sp_name = aig.replace('a12000', '')
    sci_name = sp_data.loc[sp_data.id == sp_name, 'sci_name']
    if len(sci_name) == 0:
        logger.warning('No entry found for AIG {0}'.format(sp_name))
    else:
        sci_name = sci_name.values[0]
        input_file = os.path.join(input_dir, aig)
        output_file = sci_name.lower().replace(' ', '_') + '_a12000.tif'
        output_file = os.path.join(output_dir, output_file)
        # Define potentially platform-specific, GDAL-related variables
        # constants
        cmd = ['gdal_translate', '-of', 'GTiff', '-co', 'COMPRESS=DEFLATE']
        cmd.append(input_file)
        cmd.append(output_file)
        logger.info('[{0}/{1}] Translating {2} to {3}...'.format(counter,
                                                          n_files,
                                                          input_file,
                                                          output_file))
        subprocess.run(cmd)
        counter += 1
