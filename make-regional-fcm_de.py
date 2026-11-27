#!/usr/bin/python
# -*- coding: UTF-8

"""
Build the FieldConditionModifier lookup table per NUTS3 region for the production simulations.

Every region of the weights file gets a value:
  regional  the calibrated value of the region (from collect-calibration-results_de.py)
  national  the national value, for regions without a reliable regional value: no valid calibration run,
            fewer than --min-obs observed years in the calibration period, or no pixels of the crop

Usage:
    python make-regional-fcm_de.py --best regional_best_WW.csv --obs calibration_WW_de.csv
                                   --weights Weights_WW_de.csv --national 0.93 [--min-obs 2]
                                   [--first-year 1995] [--out regional_fcm_WW_de.csv]
"""

import argparse
import csv
import sys

PARAM = "FieldConditionModifier"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--best", required=True, help="output of collect-calibration-results_de.py")
    parser.add_argument("--obs", required=True, help="observation csv, e.g. calibration_WW_de.csv")
    parser.add_argument("--weights", required=True, help="weights csv, e.g. Weights_WW_de.csv")
    parser.add_argument("--national", required=True, type=float, help="national FieldConditionModifier")
    parser.add_argument("--min-obs", type=int, default=2,
                        help="minimum number of observed years for a regional value (default 2)")
    parser.add_argument("--first-year", type=int, default=1995, help="first calibration year (default 1995)")
    parser.add_argument("--out", default="regional_fcm.csv", help="output csv (default regional_fcm.csv)")
    args = parser.parse_args()

    # observed years per region in the calibration period
    region_to_no_of_obs = {}
    with open(args.obs, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        year_cols = [i for i, h in enumerate(header) if h.strip().isdigit() and int(h) >= args.first_year]
        for row in reader:
            region_to_no_of_obs[row[1].strip()] = sum(1 for i in year_cols if float(row[i]) >= 0)

    regions = []
    with open(args.weights, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            regions.append((row[2].strip(), row[0].strip(), int(float(row[3]))))

    best = {}
    with open(args.best, newline="") as f:
        for row in csv.DictReader(f):
            best[row["region_id"]] = row
    if best and PARAM not in next(iter(best.values())):
        sys.exit(f"{args.best} has no column {PARAM}")

    rows = []
    for region_id, nuts_id, pixels in regions:
        no_of_obs = region_to_no_of_obs.get(region_id, 0)
        b = best.get(region_id)
        value, source, reason, rmse = args.national, "national", "", ""
        if pixels == 0:
            reason = "no crop pixels"
        elif b is None:
            reason = "not calibrated"
        elif not b["rmse_kg_ha"] or not b.get(PARAM):
            reason = "no valid calibration run"
        elif no_of_obs < args.min_obs:
            reason = f"only {no_of_obs} observed year(s)"
        else:
            value, source, rmse = float(b[PARAM]), "regional", b["rmse_kg_ha"]
        rows.append({"region_id": region_id, "nuts_id": nuts_id, PARAM: value, "source": source,
                     "reason": reason, "no_of_obs": no_of_obs, "crop_pixels": pixels, "rmse_kg_ha": rmse})

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: int(r["region_id"])))

    regional = [r for r in rows if r["source"] == "regional"]
    with_pixels = [r for r in rows if r["crop_pixels"] > 0]
    print(f"{len(rows)} regions written to {args.out}: {len(regional)} regional values, "
          f"{len(rows) - len(regional)} national ({args.national}), "
          f"of which {len(with_pixels) - len(regional)} with crop pixels")
    for r in with_pixels:
        if r["source"] == "national":
            print(f"  national value for {r['region_id']} ({r['nuts_id']}): {r['reason']}")


if __name__ == "__main__":
    main()
