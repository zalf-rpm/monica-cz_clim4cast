#!/usr/bin/python
# -*- coding: UTF-8

"""
Collect the best parameter set and its RMSE of all regions into one table.

Reads every <region>_SCEUA_monica_results.csv or <region>_GRID_monica_results.csv in the given output folders (e.g. the output folders of
all region batch jobs), takes the run with the lowest like1 (= RMSE in kg/ha) per region and writes
one csv row per region.

Usage:
    python collect-calibration-results_de.py OUT_FOLDER [OUT_FOLDER ...] [--out best_parameters.csv]
                                             [--obs calibration_WW_de.csv]

    --obs adds the NUTS id of each region from the observation file.
    Folders may contain wildcards, e.g. "/beegfs/rpm/projects/monica/out/palka_*".
"""

import argparse
import csv
import glob
import math
import os
import re
import sys

RESULTS_FILE_PATTERN = re.compile(r"^(?P<region>.+)_(SCEUA|GRID)_monica_results\.csv$")


def best_run(path_to_results_csv):
    """Return (no_of_runs, best_row_index, best_like1, {param_name: value}) of a spotpy results csv."""
    with open(path_to_results_csv, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        like_col = header.index("like1")
        par_cols = [(i, name[3:]) for i, name in enumerate(header) if name.startswith("par")]
        no_of_runs = 0
        best = None
        for row_index, row in enumerate(reader):
            if len(row) < len(header):
                continue  # incomplete last line of a still running job
            no_of_runs += 1
            try:
                like1 = float(row[like_col])
            except ValueError:
                continue
            if math.isnan(like1) or math.isinf(like1):
                continue
            if best is None or like1 < best[1]:
                best = (row_index, like1, {name: row[i] for i, name in par_cols})
    if best is None:
        return no_of_runs, None, None, {}
    return (no_of_runs,) + best


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folders", nargs="+", help="output folders of the calibration jobs (wildcards allowed)")
    parser.add_argument("--out", default="best_parameters.csv", help="output csv (default: best_parameters.csv)")
    parser.add_argument("--obs", help="observation csv (calibration_<crop>_de.csv) to add the NUTS ids")
    args = parser.parse_args()

    region_id_to_nuts_id = {}
    if args.obs:
        with open(args.obs, newline="") as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                region_id_to_nuts_id[row[1].strip()] = row[0].strip()

    folders = sorted({f for pattern in args.folders for f in glob.glob(pattern) if os.path.isdir(f)})
    if not folders:
        sys.exit(f"no folders found for {args.folders}")

    region_to_result = {}
    param_names = []
    for folder in folders:
        for file_name in sorted(os.listdir(folder)):
            m = RESULTS_FILE_PATTERN.match(file_name)
            if not m:
                continue
            region = m.group("region")
            no_of_runs, best_index, like1, params = best_run(os.path.join(folder, file_name))
            for name in params:
                if name not in param_names:
                    param_names.append(name)
            result = {"region_id": region, "nuts_id": region_id_to_nuts_id.get(region, ""),
                      "rmse_kg_ha": like1, "no_of_runs": no_of_runs, "best_run": best_index,
                      "folder": folder, **params}
            if region in region_to_result:
                prev = region_to_result[region]
                print(f"WARNING: region {region} found in {prev['folder']} and {folder}; "
                      f"keeping the one with the lower RMSE")
                if like1 is None or (prev["rmse_kg_ha"] is not None and prev["rmse_kg_ha"] <= like1):
                    continue
            region_to_result[region] = result

    if not region_to_result:
        sys.exit(f"no *_SCEUA_monica_results.csv or *_GRID_monica_results.csv found in {len(folders)} folder(s)")

    def sort_key(region):
        return (0, int(region)) if region.isdigit() else (1, region)

    columns = ["region_id", "nuts_id", "rmse_kg_ha", "no_of_runs", "best_run"] + param_names + ["folder"]
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for region in sorted(region_to_result, key=sort_key):
            writer.writerow(region_to_result[region])

    no_result = [r for r, res in region_to_result.items() if res["rmse_kg_ha"] is None]
    print(f"{len(region_to_result)} regions from {len(folders)} folder(s) written to {args.out}")
    if no_result:
        print(f"WARNING: no valid run for regions {sorted(no_result, key=sort_key)}")


if __name__ == "__main__":
    main()
