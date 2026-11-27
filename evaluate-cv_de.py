#!/usr/bin/python
# -*- coding: UTF-8

"""
Evaluate leave-one-year-out cross-validation of the German calibration.

Reads every calibration in the given output folders, i.e. each <prefix>_observations.csv together with
<prefix>_SCEUA_monica_results.csv or <prefix>_GRID_monica_results.csv, and predicts held-out years:

  fold       calibrations run with leave_out_years (e.g. national DE_loyo2003_...): the best run of the
             calibration (lowest like1, computed without the held-out years) predicts the held-out year(s)
  grid_loyo  grid searches on all years (e.g. regional 150_GRID_...): for every year, the grid value with
             the lowest RMSE on the other years is chosen and predicts that year; no new simulations needed

Predictions are national means (calibration_level national, weighted by region weight over the regions
with observation and simulation in that year) or region yields (nuts3).
SCE-UA calibrations without held-out years have nothing to predict and are skipped.

Usage:
    python evaluate-cv_de.py OUT_FOLDER [OUT_FOLDER ...] [--out cv_predictions.csv] [--summary cv_summary.csv]
    Folders may contain wildcards.
"""

import argparse
import csv
import glob
import math
import os
import sys

RESULTS_TYPES = [("_SCEUA_monica_results.csv", "sceua"), ("_GRID_monica_results.csv", "grid")]


def is_nan(x):
    return x != x


def read_observations(path):
    obs = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            obs.append({"region_id": row["region_id"], "year": int(row["year"]), "obs": float(row["obs_kg_ha"]),
                        "weight": float(row["weight"]), "held_out": row["held_out"] == "1",
                        "level": row["calibration_level"]})
    return obs


def read_runs(path, no_of_obs):
    """Return a list of (like1, {param: value}, [simulations]) of a results csv (incomplete lines skipped)."""
    runs = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        par_cols = [(i, name[3:]) for i, name in enumerate(header) if name.startswith("par")]
        sim_start = header.index("simulation_0")
        if sum(1 for name in header if name.startswith("simulation_")) != no_of_obs:
            raise ValueError(f"{path}: number of simulations does not match the observations file")
        for row in reader:
            if len(row) < len(header):
                continue
            runs.append((float(row[0]), {name: float(row[i]) for i, name in par_cols},
                         [float(x) for x in row[sim_start:sim_start + no_of_obs]]))
    return runs


def no_of_complete_runs(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader, [])
        return sum(1 for row in reader if len(row) >= len(header))


def valid_pairs(obs, sims, use):
    """Indices with observation and simulation available and use(observation) true."""
    return [i for i, o in enumerate(obs) if use(o) and not is_nan(o["obs"]) and not is_nan(sims[i])]


def national_means(obs, sims, use):
    """{year: (obs_mean, sim_mean)} weighted over the regions having both values in that year."""
    sums = {}
    for i in valid_pairs(obs, sims, use):
        o = obs[i]
        s = sums.setdefault(o["year"], [0.0, 0.0, 0.0])
        s[0] += o["weight"] * o["obs"]
        s[1] += o["weight"] * sims[i]
        s[2] += o["weight"]
    return {y: (s[0] / s[2], s[1] / s[2]) for y, s in sums.items() if s[2] > 0}


def objective(obs, sims, level, use):
    """Calibration objective restricted to the observations with use(o) true (mirrors the calibration)."""
    if level == "national":
        means = national_means(obs, sims, use)
        if not means:
            return math.nan
        return math.sqrt(sum((s - o) ** 2 for o, s in means.values()) / len(means))
    idx = valid_pairs(obs, sims, use)
    w = sum(obs[i]["weight"] for i in idx)
    if not idx or w <= 0:
        return math.nan
    return math.sqrt(sum(obs[i]["weight"] * (sims[i] - obs[i]["obs"]) ** 2 for i in idx) / w)


def predictions(obs, sims, level, years):
    """Observed and simulated values of the given years: national means or region yields."""
    preds = []
    if level == "national":
        means = national_means(obs, sims, lambda o: o["year"] in years)
        for y in sorted(means):
            preds.append(("DE", y, means[y][0], means[y][1]))
    else:
        for i in valid_pairs(obs, sims, lambda o: o["year"] in years):
            preds.append((obs[i]["region_id"], obs[i]["year"], obs[i]["obs"], sims[i]))
    return preds


def best_run(runs, scores):
    scored = [(s, k) for k, s in enumerate(scores) if not is_nan(s)]
    return min(scored)[1] if scored else None


def evaluate(prefix, results_path, sampler, records):
    obs = read_observations(prefix + "_observations.csv")
    runs = read_runs(results_path, len(obs))
    if not runs:
        return "no complete run"
    level = obs[0]["level"]
    source = os.path.basename(prefix)
    held_out_years = sorted({o["year"] for o in obs if o["held_out"]})

    if held_out_years:
        # fold: like1 was computed without the held-out years
        k = best_run(runs, [r[0] for r in runs])
        if k is None:
            return "no valid run"
        like1, params, sims = runs[k]
        for region_id, year, o, s in predictions(obs, sims, level, set(held_out_years)):
            records.append({"source": source, "method": "fold", "calibration_level": level,
                            "region_id": region_id, "year": year, "obs_kg_ha": o, "sim_kg_ha": s,
                            "train_rmse_kg_ha": like1, **params})
        return None

    if sampler == "grid":
        # internal leave-one-year-out over the stored grid simulations
        for year in sorted({o["year"] for o in obs if not is_nan(o["obs"])}):
            scores = [objective(obs, sims, level, lambda o, y=year: o["year"] != y) for _, _, sims in runs]
            k = best_run(runs, scores)
            if k is None:
                continue
            for region_id, y, o, s in predictions(obs, runs[k][2], level, {year}):
                records.append({"source": source, "method": "grid_loyo", "calibration_level": level,
                                "region_id": region_id, "year": y, "obs_kg_ha": o, "sim_kg_ha": s,
                                "train_rmse_kg_ha": scores[k], **runs[k][1]})
        return None

    return "SCE-UA calibration without held-out years (nothing to predict)"


def stats(rs):
    d = [r["sim_kg_ha"] - r["obs_kg_ha"] for r in rs]
    return {"n": len(d), "rmse_kg_ha": math.sqrt(sum(x * x for x in d) / len(d)), "bias_kg_ha": sum(d) / len(d),
            "mean_obs_kg_ha": sum(r["obs_kg_ha"] for r in rs) / len(rs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folders", nargs="+", help="output folders of the calibration jobs (wildcards allowed)")
    parser.add_argument("--out", default="cv_predictions.csv", help="one row per prediction")
    parser.add_argument("--summary", default="cv_summary.csv", help="RMSE and bias per level and region")
    args = parser.parse_args()

    folders = sorted({f for pattern in args.folders for f in glob.glob(pattern) if os.path.isdir(f)})
    if not folders:
        sys.exit(f"no folders found for {args.folders}")

    # the same calibration (e.g. region 150, or DE_loyo2003) may be in several folders, e.g. after a re-run
    # with a wider grid: use the results file with the most complete runs, with equal runs the newer one
    candidates = {}
    for folder in folders:
        for file_name in sorted(os.listdir(folder)):
            if not file_name.endswith("_observations.csv"):
                continue
            name = file_name[:-len("_observations.csv")]
            prefix = os.path.join(folder, name)
            found = [(prefix + suffix, sampler) for suffix, sampler in RESULTS_TYPES if os.path.exists(prefix + suffix)]
            if not found:
                print(f"skipped {prefix}: no results file")
                continue
            for results_path, sampler in found:
                key = (no_of_complete_runs(results_path), os.path.getmtime(results_path))
                candidates.setdefault((name, sampler), []).append((key, prefix, results_path))

    records = []
    for (name, sampler), files in sorted(candidates.items()):
        files.sort(reverse=True)
        _, prefix, results_path = files[0]
        if len(files) > 1:
            print(f"NOTE: {name} found in {len(files)} folders; using {os.path.dirname(results_path)} "
                  f"({files[0][0][0]} runs)")
        problem = evaluate(prefix, results_path, sampler, records)
        if problem:
            print(f"skipped {os.path.basename(results_path)}: {problem}")

    if not records:
        sys.exit("no cross-validation predictions found")

    seen = set()
    for r in records:
        key = (r["calibration_level"], r["region_id"], r["year"])
        if key in seen:
            print(f"WARNING: held-out year {r['year']} of {r['region_id']} predicted more than once")
        seen.add(key)
        r["diff_kg_ha"] = r["sim_kg_ha"] - r["obs_kg_ha"]

    param_names = sorted({k for r in records for k in r} - {
        "source", "method", "calibration_level", "region_id", "year", "obs_kg_ha", "sim_kg_ha", "diff_kg_ha",
        "train_rmse_kg_ha"})
    columns = ["source", "method", "calibration_level", "region_id", "year", "obs_kg_ha", "sim_kg_ha",
               "diff_kg_ha", "train_rmse_kg_ha"] + param_names
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda r: (r["calibration_level"], r["region_id"], r["year"])))

    groups = {}
    for r in records:
        groups.setdefault((r["calibration_level"], "ALL"), []).append(r)
        if r["calibration_level"] == "nuts3":
            groups.setdefault(("nuts3", r["region_id"]), []).append(r)
    with open(args.summary, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["calibration_level", "region_id", "n", "rmse_kg_ha", "bias_kg_ha", "mean_obs_kg_ha"])
        for (level, region_id), rs in sorted(groups.items()):
            s = stats(rs)
            writer.writerow([level, region_id, s["n"], s["rmse_kg_ha"], s["bias_kg_ha"], s["mean_obs_kg_ha"]])
            if region_id == "ALL":
                print(f"{level}: {s['n']} predictions, cross-validated RMSE {s['rmse_kg_ha']:.0f} kg/ha, "
                      f"bias {s['bias_kg_ha']:+.0f} kg/ha (mean observed {s['mean_obs_kg_ha']:.0f} kg/ha)")
    print(f"predictions written to {args.out}, summary to {args.summary}")


if __name__ == "__main__":
    main()
