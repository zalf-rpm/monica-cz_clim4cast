#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from datetime import datetime
import capnp
from collections import defaultdict
import copy
import csv
from datetime import date, timedelta
import json
import math
import numpy as np
import os
from pathlib import Path
from pyproj import CRS, Transformer
import sqlite3
import sqlite3 as cas_sq3
import sys
import time
import zmq
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from rasterio import features
import subprocess

import monica_io3
import cz_soil_io3
import monica_run_lib
import shared

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure"
PATH_TO_PYTHON_CODE = PATH_TO_MAS_INFRASTRUCTURE_REPO / "src/python"
if str(PATH_TO_PYTHON_CODE) not in sys.path:
    sys.path.insert(1, str(PATH_TO_PYTHON_CODE))

from pkgs.common import common
from pkgs.model import monica_io3

PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)

PATHS = {
    # adjust the local path to your environmen
    "ow-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
    },
    "mbm-local-remote": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
    },
    "mbm-local-local": {
        # "include-file-base-path": "/home/berg/GitHub/monica-parameters/", # path to monica-parameters
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
    },
    "hpc-local-remote": {
        # "path-to-climate-dir": "/beegfs/common/data/soil/global_soil_dataset_for_earth_system_modeling/",
        # mounted path to archive or hard drive with climate data
        "path-to-climate-dir": "/beegfs/common/data/climate/",  # mounted path to archive or hard drive with climate data
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        # mounted path to archive accessable by monica executable
        "path-to-data-dir": "./data/",  # mounted path to archive or hard drive with data
    },
}

DATA_SOIL_DB = "cz/cz_soil_500_woesten.sqlite"
SOIL_DB_URL = "https://github.com/zalf-rpm/monica-cz_clim4cast/raw/refs/heads/main/data/cz/cz_soil_500_woesten.sqlite"
DATA_GRID_HEIGHT = "cz/cz_dem_500_32633_etrs89-utm33n.asc"
DATA_GRID_SLOPE = "cz/cz_slope_500_32633_etrs89-utm33n.asc"
#DATA_GRID_LAND_USE = "germany/landuse_1000_31469_gk5.asc"
DATA_GRID_SOIL = "cz/cz_soil_500_32633_etrs89-utm33n.asc"
DATA_GRID_SOIL_OW = "germany/buek200_1000_25832_etrs89-utm32n_OW.asc"
DATA_GRID_CROPS = "cz/cz_crop-cw_500_32633_etrs89-utm33n.asc" # Added as a cropmap for winter wheat OW
# ORIGINAL DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/crops-all2017-2019_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/dwd-stations-pheno_1000_25832_etrs89-utm32n.asc"
# DATA_GRID_CROPS = "germany/germany-complete_1000_25832_etrs89-utm32n.asc"
TEMPLATE_PATH_LATLON = "{path_to_climate_dir}/latlon-to-rowcol.json"

# Additional data for masking the regions
NUTS3_REGIONS = "data/cz/cz_nuts3_32633.shp"

TEMPLATE_PATH_HARVEST = "{path_to_data_dir}/projects/monica-germany/ILR_SEED_HARVEST_doys_{crop_id}.csv"

gdf = gpd.read_file(NUTS3_REGIONS)


def run_producer(server={"server": None, "port": None}):
    context = zmq.Context()
    socket = context.socket(zmq.PUSH)  # pylint: disable=no-member
    # config_and_no_data_socket = context.socket(zmq.PUSH)

    config = {
        "mode": "mbm-local-remote",
        "server-port": server["port"] if server["port"] else "6666",
        "server": server["server"] if server["server"] else "login01.cluster.zalf.de",
        "start-row": "0",
        "end-row": "-1",
        "path_to_dem_grid": "",
        "sim.json": "sim_calibration.json",
        "crop.json": "crop_calibration.json",
        "site.json": "site.json",
        "setups-file": "sim_setups_calibration_OW.csv",
        "run-setups": "[1]",
        "reader_sr": None,
        "path_to_out": "out/",
        "only_nuts3_region_ids": "[]",  # "[10]",
    }


    common.update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    path_to_out_file = config["path_to_out"] + "/producer.out"
    if not os.path.exists(config["path_to_out"]):
        try:
            os.makedirs(config["path_to_out"])
        except OSError:
            print("run-calibration-producer.py: Couldn't create dir:", config["path_to_out"], "!")
    with open(path_to_out_file, "a") as _:
        _.write(f"config: {config}\n")

    with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        _.write(f"{datetime.now()} start producer in producer\n") 

    nuts3_region_ids = json.loads(config["only_nuts3_region_ids"])

     # select paths
    paths = PATHS[config["mode"]]

    soil_db_path = paths["path-to-data-dir"] + DATA_SOIL_DB
    subprocess.run(["wget", "-O", soil_db_path, SOIL_DB_URL], check=True)
    print("Downloaded soil db successfully.")

    # open soil db connection
    # soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB)
    soil_db_con = sqlite3.connect(soil_db_path)
    print("Connected to soil db successfully.")
    # soil_db_con = cas_sq3.connect(paths["path-to-data-dir"] + DATA_SOIL_DB) #CAS.
    # connect to monica proxy (if local, it will try to connect to a locally started monica)
    socket.connect("tcp://" + config["server"] + ":" + str(config["server-port"]))


    # read setup from csv file
    setups = monica_run_lib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    print("read sim setups: ", config["setups-file"])

    with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
        _.write(f"{datetime.now()} setup read\n") 

    # transforms geospatial coordinates from one coordinate reference system to another
    # transform wgs84 into gk5
    soil_crs_to_x_transformers = {}
    wgs84_crs = CRS.from_epsg(4326)
    utm32_crs = CRS.from_epsg(32633)
    # transformers[wgs84] = Transformer.from_crs(wgs84_crs, gk5_crs, always_xy=True)

    ilr_seed_harvest_data = defaultdict(
        lambda: {"interpolate": None, "data": defaultdict(dict), "is-winter-crop": None})

    # Load grids
    ## note numpy is able to load from a compressed file, ending with .gz or .bz2

    # soil data
    path_to_soil_grid = paths["path-to-data-dir"] + DATA_GRID_SOIL
    soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[3])
    soil_crs = CRS.from_epsg(soil_epsg_code)
    if wgs84_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[wgs84_crs] = Transformer.from_crs(soil_crs, wgs84_crs)
    soil_metadata, _ = Mrunlib.read_header(path_to_soil_grid)
    soil_grid_original = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    print("read: ", path_to_soil_grid)

    # height data
    path_to_dem_grid = paths["path-to-data-dir"] + DATA_GRID_HEIGHT
    dem_epsg_code = int(path_to_dem_grid.split("/")[-1].split("_")[3])
    dem_crs = CRS.from_epsg(dem_epsg_code)
    if dem_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[dem_crs] = Transformer.from_crs(soil_crs, dem_crs)
    dem_metadata, _ = Mrunlib.read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=6)
    dem_interpolate = Mrunlib.create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid)

    # slope data
    path_to_slope_grid = paths["path-to-data-dir"] + DATA_GRID_SLOPE
    slope_epsg_code = int(path_to_slope_grid.split("/")[-1].split("_")[3])
    slope_crs = CRS.from_epsg(slope_epsg_code)
    if slope_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[slope_crs] = Transformer.from_crs(soil_crs, slope_crs)
    slope_metadata, _ = Mrunlib.read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_interpolate = Mrunlib.create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid)

    # crop mask data
    path_to_crop_grid = paths["path-to-data-dir"] + DATA_GRID_CROPS
    crop_epsg_code = int(path_to_crop_grid.split("/")[-1].split("_")[3])
    crop_crs = CRS.from_epsg(crop_epsg_code)
    if crop_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[crop_crs] = Transformer.from_crs(soil_crs, crop_crs)
    crop_meta, _ = Mrunlib.read_header(path_to_crop_grid)
    crop_grid = np.loadtxt(path_to_crop_grid, dtype=int, skiprows=6)
    crop_interpolate = Mrunlib.create_ascii_grid_interpolator(crop_grid, crop_meta)
    print("read: ", path_to_crop_grid)

    # nuts3_regions
    path_to_nuts3_regions_grid = paths["path-to-data-dir"] + "cz/nuts3-regions_500_32633_etrs89-utm33n.asc"
    nuts3_regions_epsg_code = int(path_to_nuts3_regions_grid.split("/")[-1].split("_")[2])
    nuts3_regions_crs = CRS.from_epsg(nuts3_regions_epsg_code)
    if nuts3_regions_crs not in soil_crs_to_x_transformers:
        soil_crs_to_x_transformers[nuts3_regions_crs] = Transformer.from_crs(soil_crs, nuts3_regions_crs)
    nuts3_regions_metadata, _ = monica_run_lib.read_header(path_to_nuts3_regions_grid)
    nuts3_regions_grid = np.loadtxt(path_to_nuts3_regions_grid, dtype=float, skiprows=6)
    nuts3_regions_interpolate = monica_run_lib.create_ascii_grid_interpolator(nuts3_regions_grid, nuts3_regions_metadata)
    print("read: ", path_to_nuts3_regions_grid)



    # Create the function for the mask. This function will later use the additional column in a setup file!
    def create_mask_from_shapefile(NUTS3_REGIONS, region_name, path_to_soil_grid):
        regions_df = gpd.read_file(NUTS3_REGIONS)
        region = regions_df[regions_df["NUTS_ID"] == region_name]

        # This is needed to read the transformation data correctly from the file. With the original opening it does not work
        with rasterio.open(path_to_soil_grid) as dataset:
            soil_grid = dataset.read(1)
            transform = dataset.transform

        rows, cols = soil_grid.shape
        mask = rasterio.features.geometry_mask([region.geometry.values[0]], out_shape=(rows, cols), transform=transform,
                                               invert=True)

        return mask

    if len(run_setups) > 1 and run_setups[0] not in setups:
        return
    else:
        setup_id = run_setups[0]

    conman = common.ConnectionManager()
    reader = conman.try_connect(config["reader_sr"], cast_as=fbp_capnp.Channel.Reader, retry_secs=1)
    if reader:
        sent_env_count = 0
        while True:
            msg = reader.read().wait()
            # check for end of data from in port
            if msg.which() == "done":
                break

            with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
                _.write(f"{datetime.now()} connected\n") 

            env_template = None
            start_setup_time = None
            try:
                in_ip = msg.value.as_struct(fbp_capnp.IP)
                s: str = in_ip.content.as_text()
                params = json.loads(s)  # keys: MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
                if "only_nuts3_region_ids" in params:
                    nuts3_region_ids = params["only_nuts3_region_ids"]
                    del params["only_nuts3_region_ids"]

                start_setup_time = time.perf_counter()
            

                setup = setups[setup_id]
                crop_id = setup["crop-id"]
                #region_name = setup["region_name"]

                ## extract crop_id from crop-id name that has possible an extenstion
                crop_id_short = crop_id.split('_')[0]

                #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
                #    _.write(f"{datetime.now()} setup started producer\n") 

                #if region_name and len(region_name) > 0:
                    # Create the soil mask for the specific region
                    #path_to_soil_grid_ow = paths["path-to-data-dir"] + DATA_GRID_SOIL_OW
                    #mask = create_mask_from_shapefile(NUTS3_REGIONS, region_name, path_to_soil_grid_ow)

                    # Apply the soil mask to the soil grid
                    #soil_grid_copy = soil_grid.copy()
                    #soil_grid[mask == False] = -8888
                    #soil_grid[soil_grid_copy == -9999] = -9999

                # add crop id from setup file
                try:
                    # read seed/harvest dates for each crop_id
                    path_harvest = TEMPLATE_PATH_HARVEST.format(path_to_data_dir=paths["path-to-data-dir"],
                                                                crop_id=crop_id_short)
                    print("created seed harvest gk5 interpolator and read data: ", path_harvest)
                    monica_run_lib.create_seed_harvest_geoGrid_interpolator_and_read_data(path_harvest, wgs84_crs, utm32_crs,
                                                                                          ilr_seed_harvest_data)
                except IOError:
                    path_harvest = TEMPLATE_PATH_HARVEST.format(path_to_data_dir=paths["path-to-data-dir"],
                                                                crop_id=crop_id_short)
                    print("Couldn't read file:", path_harvest)
                    continue

                #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
                #    _.write(f"{datetime.now()} crop added producer\n") 

                cdict = {}
                # path to latlon-to-rowcol.json
                # path = TEMPLATE_PATH_LATLON.format(path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
                path = TEMPLATE_PATH_LATLON.format(
                    path_to_climate_dir=paths["path-to-climate-dir"] + setup["climate_path_to_latlon_file"] + "/")
                climate_data_interpolator = monica_run_lib.create_climate_geoGrid_interpolator_from_json_file(path, wgs84_crs,
                                                                                                              soil_crs, cdict)
                print("created climate_data to gk5 interpolator: ", path)

                #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
                #    _.write(f"{datetime.now()} climate data read producer\n") 

                # read template sim.json
                with open(setup.get("sim.json", config["sim.json"])) as _:
                    sim_json = json.load(_)
                # change start and end date according to setup
                if setup["start_date"]:
                    sim_json["climate.csv-options"]["start-date"] = str(setup["start_date"])
                if setup["end_date"]:
                    sim_json["climate.csv-options"]["end-date"] = str(setup["end_date"])
                    # sim_json["include-file-base-path"] = paths["include-file-base-path"]

                # read template site.json
                with open(setup.get("site.json", config["site.json"])) as _:
                    site_json = json.load(_)

                #with open(config["path_to_out"] + "/spot_setup.out", "a") as _:
                #    _.write(f"{datetime.now()} read site and sim json producer\n\n") 

                #site_json["EnvironmentParameters"]["rcp"] = scenario

                # read template crop.json
                with open(setup.get("crop.json", config["crop.json"])) as _:
                    crop_json = json.load(_)
                    crop_json["cropRotation"][2] = crop_id
                    real_crop_id = None
                    # set value of calibration params
                    for ws in crop_json["cropRotationTemplates"][crop_id][0]["worksteps"]:
                        if "Sowing" in ws["type"]:
                            real_crop_id = ws["crop"][2]
                    if real_crop_id:
                        ps = crop_json["crops"][real_crop_id]["cropParams"]
                        for pname, pval in params.items():
                            if pname == "DaylengthRequirement":
                                ps["cultivar"][pname][1] = pval
                                ps["cultivar"][pname][2] = pval
                                ps["cultivar"][pname][3] = pval
                            elif pname == "SpecificLeafArea":
                                ps["cultivar"][pname][0] *= pval
                                ps["cultivar"][pname][1] *= pval
                                ps["cultivar"][pname][2] *= pval
                                ps["cultivar"][pname][3] *= pval
                                ps["cultivar"][pname][4] *= pval
                                ps["cultivar"][pname][5] *= pval
                            elif pname == "BaseDaylength":
                                ps["cultivar"][pname][2] = pval
                                ps["cultivar"][pname][3] = pval
                            else:
                                pname_arr = pname.split("_")
                                i = None
                                if len(pname_arr) == 2:
                                    pname = pname_arr[0]
                                    i = int(pname_arr[1])
                                if pname in ps["species"]:
                                    if i:
                                        if len(ps["species"][pname]) < i:
                                            ps["species"][pname][i] = pval
                                    else:
                                        ps["species"][pname] = pval
                                elif pname in ps["cultivar"]:
                                    if i:
                                        if len(ps["cultivar"][pname]) > i:
                                            ps["cultivar"][pname][i] = pval
                                    else:
                                        ps["cultivar"][pname] = pval
                    else:
                        print("Error couldn't find sowing workstep in crop.json")
                        exit(1)

                crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup[
                    "use_vernalisation_fix"] if "use_vernalisation_fix" in setup else False

                # create environment template from json templates
                env_template = monica_io3.create_env_json_from_json_config({
                    "crop": crop_json,
                    "site": site_json,
                    "sim": sim_json,
                    "climate": ""
                })

                scols = int(soil_metadata["ncols"])
                srows = int(soil_metadata["nrows"])
                scellsize = int(soil_metadata["cellsize"])
                xllcorner = int(soil_metadata["xllcorner"])
                yllcorner = int(soil_metadata["yllcorner"])
                nodata_value = int(soil_metadata["nodata_value"])

                # unknown_soil_ids = set()
                soil_id_cache = {}
                print("All Rows x Cols: " + str(srows) + "x" + str(scols))
                for srow in range(0, srows):
                    print(srow, end=", ")

                    if srow < int(config["start-row"]):
                        continue
                    elif int(config["end-row"]) > 0 and srow > int(config["end-row"]):
                        break

                    for scol in range(0, scols):
                        soil_id = int(soil_grid[srow, scol])
                        if soil_id == nodata_value:
                            continue

                        # get coordinate of clostest climate element of real soil-cell
                        sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                        sr = xllcorner + (scellsize / 2) + scol * scellsize
                        # inter = crow/ccol encoded into integer
                        crow, ccol = climate_data_interpolator(sr, sh)

                        crop_grid_id = int(crop_grid[srow, scol])
                        # print(crop_grid_id)
                        if crop_grid_id != 1 or soil_id == -8888:
                            continue

                        tcoords = {}

                        if nuts3_regions_crs not in tcoords:
                            tcoords[nuts3_regions_crs] = soil_crs_to_x_transformers[nuts3_regions_crs].transform(sr, sh)
                        n3r, n3h = tcoords[nuts3_regions_crs]
                        nuts3_region_id = int(nuts3_regions_interpolate(n3r, n3h))
                        if not nuts3_region_id or (len(nuts3_region_ids) > 0 and nuts3_region_id not in nuts3_region_ids):
                            continue

                        if soil_id in soil_id_cache:
                            soil_profile = soil_id_cache[soil_id]
                        else:
                            soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
                            soil_id_cache[soil_id] = soil_profile
                        if not soil_profile or len(soil_profile) == 0:
                            continue

                        worksteps = env_template["cropRotation"][0]["worksteps"]
                        sowing_ws = next(filter(lambda ws: ws["type"][-6:] == "Sowing", worksteps))
                        harvest_ws = next(filter(lambda ws: ws["type"][-7:] == "Harvest", worksteps))

                        ilr_interpolate = ilr_seed_harvest_data[crop_id_short]["interpolate"]
                        seed_harvest_cs = ilr_interpolate(sr, sh) if ilr_interpolate else None

                        # set external seed/harvest dates
                        
                        if seed_harvest_cs:
                            seed_harvest_data = ilr_seed_harvest_data[crop_id_short]["data"][seed_harvest_cs]
                            if seed_harvest_data:
                                is_winter_crop = ilr_seed_harvest_data[crop_id_short]["is-winter-crop"]

                                if setup[
                                    "sowing-date"] == "fixed":  # fixed indicates that regionally fixed sowing dates will be used
                                    sowing_date = seed_harvest_data["sowing-date"]
                                elif setup[
                                    "sowing-date"] == "auto":  # auto indicates that automatic sowng dates will be used that vary between regions
                                    sowing_date = seed_harvest_data["latest-sowing-date"]
                                elif setup[
                                    "sowing-date"] == "fixed1":  # fixed1 indicates that a fixed sowing date will be used that is the same for entire germany
                                    sowing_date = sowing_ws["date"]

                                sds = [int(x) for x in sowing_date.split("-")]
                                sd = date(2001, sds[1], sds[2])
                                sdoy = sd.timetuple().tm_yday

                                if setup[
                                    "harvest-date"] == "fixed":  # fixed indicates that regionally fixed harvest dates will be used
                                    harvest_date = seed_harvest_data["harvest-date"]
                                elif setup[
                                    "harvest-date"] == "auto":  # auto indicates that automatic harvest dates will be used that vary between regions
                                    harvest_date = seed_harvest_data["latest-harvest-date"]
                                elif setup[
                                    "harvest-date"] == "auto1":  # fixed1 indicates that a fixed harvest date will be used that is the same for entire germany
                                    harvest_date = harvest_ws["latest-date"]

                                hds = [int(x) for x in harvest_date.split("-")]
                                hd = date(2001, hds[1], hds[2])
                                hdoy = hd.timetuple().tm_yday

                                esds = [int(x) for x in seed_harvest_data["earliest-sowing-date"].split("-")]
                                esd = date(2001, esds[1], esds[2])

                                # sowing after harvest should probably never occur in both fixed setup!
                                if setup["sowing-date"] == "fixed" and setup["harvest-date"] == "fixed":
                                    # calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy-1))
                                    if is_winter_crop:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                                    else:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                                    sowing_ws["date"] = seed_harvest_data["sowing-date"]
                                    harvest_ws["date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                       calc_harvest_date.day)
                                elif setup["sowing-date"] == "fixed" and setup["harvest-date"] == "auto":
                                    if is_winter_crop:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                                    else:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                                    sowing_ws["date"] = seed_harvest_data["sowing-date"]
                                    harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                              calc_harvest_date.day)
                                elif setup["sowing-date"] == "fixed" and setup["harvest-date"] == "auto1":
                                    if is_winter_crop:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                                    else:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                                    sowing_ws["date"] = seed_harvest_data["sowing-date"]
                                    harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], hds[1], hds[2])

                                elif setup["sowing-date"] == "auto" and setup["harvest-date"] == "fixed":
                                    sowing_ws["earliest-date"] = seed_harvest_data["earliest-sowing-date"] if esd > date(
                                        esd.year, 6, 20) else "{:04d}-{:02d}-{:02d}".format(sds[0], 6, 20)
                                    calc_sowing_date = date(2000, 12, 31) + timedelta(days=max(hdoy + 1, sdoy))
                                    sowing_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(sds[0], calc_sowing_date.month,
                                                                                             calc_sowing_date.day)
                                    harvest_ws["date"] = seed_harvest_data["harvest-date"]

                                elif setup["sowing-date"] == "auto" and setup["harvest-date"] == "auto":
                                    sowing_ws["earliest-date"] = seed_harvest_data["earliest-sowing-date"] if esd > date(
                                        esd.year, 6, 20) else "{:04d}-{:02d}-{:02d}".format(sds[0], 6, 20)
                                    if is_winter_crop:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                                    else:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                                    sowing_ws["latest-date"] = seed_harvest_data["latest-sowing-date"]
                                    harvest_ws["latest-date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                              calc_harvest_date.day)
                                elif setup["sowing-date"] == "fixed1" and setup["harvest-date"] == "fixed":
                                    # calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy-1))
                                    if is_winter_crop:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=min(hdoy, sdoy - 1))
                                    else:
                                        calc_harvest_date = date(2000, 12, 31) + timedelta(days=hdoy)
                                    sowing_ws["date"] = sowing_date
                                    # print(seed_harvest_data["sowing-date"])
                                    harvest_ws["date"] = "{:04d}-{:02d}-{:02d}".format(hds[0], calc_harvest_date.month,
                                                                                       calc_harvest_date.day)
                        # check if current grid cell is used for agriculture
                        if setup["landcover"]:
                            if landuse_crs not in tcoords:
                                tcoords[landuse_crs] = soil_crs_to_x_transformers[landuse_crs].transform(sr, sh)
                            lur, luh = tcoords[landuse_crs]
                            landuse_id = landuse_interpolate(lur, luh)
                            if landuse_id not in [2, 3, 4]:
                                continue

                        if dem_crs not in tcoords:
                            tcoords[dem_crs] = soil_crs_to_x_transformers[dem_crs].transform(sr, sh)
                        demr, demh = tcoords[dem_crs]
                        height_nn = dem_interpolate(demr, demh)

                        if slope_crs not in tcoords:
                            tcoords[slope_crs] = soil_crs_to_x_transformers[slope_crs].transform(sr, sh)
                        slr, slh = tcoords[slope_crs]
                        slope = slope_interpolate(slr, slh)

                        env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup[
                            "LeafExtensionModifier"]

                        # print("soil:", soil_profile)
                        env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

                        # setting groundwater level
                        if setup["groundwater-level"]:
                            groundwaterlevel = 20
                            layer_depth = 0
                            for layer in soil_profile:
                                if layer.get("is_in_groundwater", False):
                                    groundwaterlevel = layer_depth
                                    # print("setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m")
                                    break
                                layer_depth += monica_run_lib.get_value(layer["Thickness"])
                            env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
                            env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [
                                max(0, groundwaterlevel - 0.2), "m"]
                            env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [
                                groundwaterlevel + 0.2, "m"]

                        # setting impenetrable layer
                        if setup["impenetrable-layer"]:
                            impenetrable_layer_depth = monica_run_lib.get_value(
                                env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
                            layer_depth = 0
                            for layer in soil_profile:
                                if layer.get("is_impenetrable", False):
                                    impenetrable_layer_depth = layer_depth
                                    # print("setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m")
                                    break
                                layer_depth += monica_run_lib.get_value(layer["Thickness"])
                            env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = \
                                [impenetrable_layer_depth, "m"]
                            env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = \
                                [impenetrable_layer_depth, "m"]

                        if setup["elevation"]:
                            env_template["params"]["siteParameters"]["heightNN"] = float(height_nn)

                        if setup["slope"]:
                            env_template["params"]["siteParameters"]["slope"] = slope / 100.0

                        if setup["latitude"]:
                            clat, _ = cdict[(crow, ccol)]
                            env_template["params"]["siteParameters"]["Latitude"] = clat

                        if setup["CO2"]:
                            env_template["params"]["userEnvironmentParameters"]["AtmosphericCO2"] = float(setup["CO2"])

                        if setup["O3"]:
                            env_template["params"]["userEnvironmentParameters"]["AtmosphericO3"] = float(setup["O3"])

                        #if setup["FieldConditionModifier"]:
                            #env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["species"][
                                #"FieldConditionModifier"] = float(setup["FieldConditionModifier"])

                        if setup["StageTemperatureSum"]:
                            stage_ts = setup["StageTemperatureSum"].split('_')
                            stage_ts = [int(temp_sum) for temp_sum in stage_ts]
                            orig_stage_ts = env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                                "StageTemperatureSum"][0]
                            if len(stage_ts) != len(orig_stage_ts):
                                stage_ts = orig_stage_ts
                                print('The provided StageTemperatureSum array is not '
                                      'sufficiently long. Falling back to original StageTemperatureSum')

                            env_template["cropRotation"][0]["worksteps"][0]["crop"]["cropParams"]["cultivar"]["="][
                                "StageTemperatureSum"][0] = stage_ts

                        env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup[
                            "fertilization"]
                        env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

                        env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]

                env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]

                # crow = int(crow)
                # ccol = int(ccol)

                # subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario=scenario, ensmem=ensmem,
                #                                                   version=version, crow=str(int(crow)),
                #                                                   ccol=str(int(ccol)))
                # for _ in range(4):
                #     subpath_to_csv = subpath_to_csv.replace("//", "/")
                # env_template["pathToClimateCSV"] = [
                #     paths["monica-path-to-climate-dir"] + setup["climate_path_to_csvs"] + subpath_to_csv]

                # climate_csv_path = paths["monica-path-to-climate-dir"] + setup["climate_path_to_csvs"] + subpath_to_csv
                # env_template["pathToClimateCSV"] = [climate_csv_path]
                climate_csv_path = (paths["monica-path-to-climate-dir"] +
                                    f"czechglobe/hist_csv_1961-01-01_to_2023-01-01/row-{crow}/col-{ccol}.csv.gz")

                env_template["pathToClimateCSV"] = climate_csv_path
                print("pathToClimateCSV:", env_template["pathToClimateCSV"])

                # if setup["incl_hist"]:
                #
                #     if rcm[:3] == "UHO":
                #         hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm="CLMcom-CCLM4-8-17",
                #                                                                scenario="historical", ensmem=ensmem,
                #                                                                version=version, crow=str(crow),
                #                                                                ccol=str(ccol))
                #         for _ in range(4):
                #             hist_subpath_to_csv = hist_subpath_to_csv.replace("//", "/")
                #         env_template["pathToClimateCSV"].insert(0, paths["monica-path-to-climate-dir"] + setup[
                #             "climate_path_to_csvs"] + "/" + hist_subpath_to_csv)
                #
                #     elif rcm[:3] == "SMH":
                #         hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm="CLMcom-CCLM4-8-17",
                #                                                                scenario="historical", ensmem=ensmem,
                #                                                                version=version, crow=str(crow),
                #                                                                ccol=str(ccol))
                #         for _ in range(4):
                #             hist_subpath_to_csv = hist_subpath_to_csv.replace("//", "/")
                #         env_template["pathToClimateCSV"].insert(0, paths["monica-path-to-climate-dir"] + setup[
                #             "climate_path_to_csvs"] + "/" + hist_subpath_to_csv)
                #
                #     hist_subpath_to_csv = TEMPLATE_PATH_CLIMATE_CSV.format(gcm=gcm, rcm=rcm, scenario="historical",
                #                                                            ensmem=ensmem, version=version,
                #                                                            crow=str(crow), ccol=str(ccol))
                #     for _ in range(4):
                #         hist_subpath_to_csv = hist_subpath_to_csv.replace("//", "/")
                #     env_template["pathToClimateCSV"].insert(0, paths["monica-path-to-climate-dir"] + setup[
                #         "climate_path_to_csvs"] + "/" + hist_subpath_to_csv)
                # print("pathToClimateCSV:", env_template["pathToClimateCSV"])
                # if DEBUG_WRITE_CLIMATE:
                #     listOfClimateFiles.add(subpath_to_csv)

                env_template["customId"] = {
                    "setup_id": setup_id,
                    "srow": srow, "scol": scol,
                    "crow": int(crow), "ccol": int(ccol),
                    "soil_id": soil_id,
                    "env_id": sent_env_count,
                    "is_sensitivity_analysis": is_sensitivity_analysis,
                    "param_name": p_name,
                    "param_value": p_value,
                    "nodata": False
                }

                # print("Harvest type:", setup["harvest-date"])
                # print("Srow: ", env_template["customId"]["srow"], "Scol:", env_template["customId"]["scol"])
                # harvest_ws = next(
                #     filter(lambda ws: ws["type"][-7:] == "Harvest", env_template["cropRotation"][0]["worksteps"]))
                # if setup["harvest-date"] == "fixed":
                #     print("Harvest-date:", harvest_ws["date"])
                # elif setup["harvest-date"] == "auto":
                #     print("Harvest-date:", harvest_ws["latest-date"])

                if not DEBUG_DONOT_SEND:
                    socket.send_json(env_template)
                    print("sent env ", sent_env_count, " customId: ", env_template["customId"])

                sent_env_count += 1
                # write debug output, as json file
                if DEBUG_WRITE:
                    debug_write_folder = paths["path-debug-write-folder"]
                    if not os.path.exists(debug_write_folder):
                        os.makedirs(debug_write_folder)
                    if sent_env_count < DEBUG_ROWS:

                        path_to_debug_file = debug_write_folder + "/row_" + str(sent_env_count - 1) + "_" + str(
                            setup_id) + ".json"

                        if not os.path.isfile(path_to_debug_file):
                            with open(path_to_debug_file, "w") as _:
                                _.write(json.dumps(env_template))
                        else:
                            print("WARNING: Row ", (sent_env_count - 1), " already exists")
            # print("unknown_soil_ids:", unknown_soil_ids)
            except Exception as e:
                with open(path_to_out_file, "a") as _:
                    _.write(f"raised exception: {e}\n")
                print("Exception raised:", e)
                raise e            

        if env_template and is_sensitivity_analysis:
            env_template["pathToClimateCSV"] = ""
            env_template["customId"] = {
                "setup_id": setup_id,
                "no_of_sent_envs": sent_env_count,
                "is_sensitivity_analysis": is_sensitivity_analysis,
            }
            socket.send_json(env_template)

            # print("crows/cols:", crows_cols)
        # cs__.close()
        stop_setup_time = time.perf_counter()
        print("\nSetup ", sent_env_count, " envs took ", (stop_setup_time - start_setup_time), " seconds")
        sent_env_count = 0

    stop_time = time.perf_counter()

    # write summary of used json files
    if DEBUG_WRITE_CLIMATE:
        debug_write_folder = paths["path-debug-write-folder"]
        if not os.path.exists(debug_write_folder):
            os.makedirs(debug_write_folder)

        path_to_climate_summary = debug_write_folder + "/climate_file_list" + ".csv"
        with open(path_to_climate_summary, "w") as _:
            _.write('\n'.join(listOfClimateFiles))

    try:
        print("sending ", (sent_env_count - 1), " envs took ", (stop_time - start_time), " seconds")
        # print("ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
        print("exiting run_producer()")
    except Exception:
        raise


if __name__ == "__main__":
    run_producer()
