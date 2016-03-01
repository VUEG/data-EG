#! /usr/bin/env pyhton3
# ONLY TO BE CALLED FROM SNAKEMAKE!!!
import errno
import os
import shutil

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

rasters = snakemake.input["src_rasters"]
rasters.sort()
renamed_rasters = [raster.replace(" ", "_").lower() for raster in rasters]
renamed_rasters = [raster.replace("org/", "") for raster in renamed_rasters]

# We'll need to manually create the final directory
mkdir_p(snakemake.output[0])

for raster in rasters:
    old_raster = raster
    # There might be a whitespace at the end of the basename
    new_raster = old_raster.replace(" .", "")
    # Replace the rest with underscores
    new_raster = new_raster.replace(" ", "_")
    # Finally, to lower case
    new_raster = new_raster.lower()
    # Set the location to be copied to
    new_raster = os.path.join(snakemake.output[0],
                              os.path.basename(new_raster))
    print("Renaming {0} -> {1}".format(old_raster, new_raster))
    shutil.copy2(old_raster, new_raster)
