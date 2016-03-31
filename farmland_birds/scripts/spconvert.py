#! /usr/bin/env pyhton3

# Run pip3 install -r requirements.txt first
#import snakemake
import sys
from pandas import concat, read_excel


# Load data
excel_file = snakemake.input[0] if 'snakemake' in sys.modules and hasattr(snakemake, "input") else "data/birds/ReadMe.xlsx"
# Skip the first row which has sum information for different sub-groups. Also,
# leave out first 2 cols; 1st is the same as 3rd, 2nd is NA.
bird_names = read_excel(excel_file, sheetname="BirdList", skiprows=1,
                        parse_cols = [2, 3, 4, 5])
# Fix column headers
bird_names.columns = ["id", "sci_name", "annex_I_spp", "farmland_spp"]

# Write output
output_file = snakemake.output[0] if 'snakemake' in sys.modules and hasattr(snakemake, "output") else "data/birds/farmland_birds_sp.csv"
bird_names.to_csv(output_file, index=False)
