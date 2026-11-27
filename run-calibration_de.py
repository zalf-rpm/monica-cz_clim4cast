
import atexit
from datetime import datetime
import capnp
from collections import defaultdict
import json
import csv
import itertools
import matplotlib.pyplot as plt
import monica_run_lib
import numpy as np
import os
from pathlib import Path
import spotpy
import subprocess as sp
import sys
import time
import uuid

import calibration_spotpy_setup_MONICA

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure"
PATH_TO_PYTHON_CODE = PATH_TO_MAS_INFRASTRUCTURE_REPO / "src/python"
if str(PATH_TO_PYTHON_CODE) not in sys.path:
    sys.path.insert(1, str(PATH_TO_PYTHON_CODE))

from pkgs.common import common

PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)


def get_reader_writer_srs_from_channel(path_to_channel_binary, chan_name=None):
    chan = sp.Popen([
        path_to_channel_binary,
        "--name=chan_{}".format(chan_name if chan_name else str(uuid.uuid4())),
        "--output_srs",
    ], stdout=sp.PIPE, text=True)
    reader_sr = None
    writer_sr = None
    while True:
        s = chan.stdout.readline().split("=", maxsplit=1)
        id, sr = s if len(s) == 2 else (None, None)
        if id and id == "readerSR":
            reader_sr = sr.strip()
        elif id and id == "writerSR":
            writer_sr = sr.strip()
        if reader_sr and writer_sr:
            break
    return {"chan": chan, "reader_sr": reader_sr, "writer_sr": writer_sr}


def balanced_region_batches(region_ids, region_pixels, no_of_batches, overhead_pixels=100):
    """
    Split regions into no_of_batches batches with similar expected run time.

    The cost of a region is its number of crop pixels (MONICA runs per iteration) plus a fixed
    overhead_pixels for the per-iteration costs that don't depend on the region size. Regions are
    assigned largest first to the currently cheapest batch. The result is deterministic, so every
    HPC job computes the same split and picks its own batch.
    """
    batches = [[] for _ in range(no_of_batches)]
    loads = [0] * no_of_batches
    for id_ in sorted(region_ids, key=lambda i: (-region_pixels[i], i)):
        b = min(range(no_of_batches), key=lambda k: (loads[k], k))
        batches[b].append(id_)
        loads[b] += region_pixels[id_] + overhead_pixels
    return [sorted(batch) for batch in batches]


def write_observations_file(path, observations, weights, calibration_level):
    """
    Write which observation belongs to which simulation_<i> column of the results file: index, region id,
    year, observed yield (kg/ha, also for held-out years), region weight, held out (1/0), calibration level.
    """
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "region_id", "year", "obs_kg_ha", "weight", "held_out", "calibration_level"])
        for i, o in enumerate(observations):
            writer.writerow([i, o["id"], o["year"], o["observed"], weights[o["id"]], int(o["held_out"]),
                             calibration_level])


class GridParameterVector(list):
    """Parameter values plus their names, as spot_setup.simulation expects from spotpy."""
    def __init__(self, names, values):
        super().__init__(values)
        self.name = names


def grid_search(spot_setup, params, grid_points, path_to_results_csv, path_to_best_out, path_to_plot):
    """
    Evaluate all combinations of grid_points equally spaced values between low and high of each calibrated
    parameter. Writes a results csv in the same layout as spotpy's (like1, par<name>..., simulation_<i>...,
    chain), the best parameter set and a plot of the objective function.
    """
    free_params = [p for p in params if "derive_function" not in p]
    names = [p["name"] for p in free_params]
    axes = [[round(float(v), 10) for v in np.linspace(p["low"], p["high"], grid_points)] for p in free_params]
    no_of_runs = grid_points ** len(free_params)
    if no_of_runs > 10000:
        raise ValueError(f"grid search with {len(free_params)} parameters and {grid_points} points "
                         f"would need {no_of_runs} runs")

    evaluation = spot_setup.evaluation()
    runs = []
    with open(path_to_results_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["like1"] + [f"par{n}" for n in names] +
                        [f"simulation_{i}" for i in range(len(evaluation))] + ["chain"])
        for values in itertools.product(*axes):
            simulation = spot_setup.simulation(GridParameterVector(names, list(values)))
            like1 = spot_setup.objectivefunction(simulation, evaluation)
            writer.writerow([like1] + list(values) + list(simulation) + [0])
            f.flush()
            runs.append((like1, values))

    best_like1, best_values = min(runs, key=lambda r: r[0])
    with open(path_to_best_out, "a") as f:
        f.write(f"*** Grid search summary ***\nTotal runs: {len(runs)}\n"
                f"Minimal objective value: {best_like1:g}\nCorresponding parameter setting:\n")
        for name, value in zip(names, best_values):
            f.write(f"par{name}: {value:g}\n")

    fig = plt.figure(1, figsize=(9, 6))
    if len(names) == 1:
        plt.plot([r[1][0] for r in runs], [r[0] for r in runs], "r+-")
        plt.xlabel(names[0])
    else:
        plt.plot([r[0] for r in runs], "r+")
        plt.xlabel("Run")
    plt.ylabel("RMSE (kg/ha)")
    fig.savefig(path_to_plot, dpi=150)
    plt.close(fig)
    return best_like1, dict(zip(names, best_values))


local_run = False


def run_calibration(server=None, prod_port=None, cons_port=None):
    config = {
        "mode": "mbm-local-remote",
        "prod-port": prod_port if prod_port else "6666",  # local: 6667, remote 6666
        "cons-port": cons_port if cons_port else "7777",  # local: 6667, remote 6666
        "server": server if server else "login01.cluster.zalf.de",
        "sim.json": "sim.json",
        "crop.json": "crop.json",
        "site.json": "site.json",
        "setups-file": "sim_setups_calibration_de.csv",
        "path_to_out": "out/",
        "run-setups": "[1]", ## Define on rundeck ##
        "path_to_channel": "/home/berg/GitHub/mas-infrastructure/src/cpp/common/_cmake_debug/channel" if local_run else
        "/home/rpm/start_manual_test_services/GitHub/mas-infrastructure/src/cpp/common/_cmake_release/channel",
        "path_to_python": "python" if local_run else "/home/rpm/.conda/envs/clim4cast/bin/python",
        "repetitions": "500", ## Define on rundeck ## SCE-UA needs >> ngs * (2 * no_of_params + 1) to evolve
        "ngs": "4",  # number of SCE-UA complexes; ~no_of_params is usually enough
        # simulations start with the setup's start_date (1990); the years before are spin-up and not calibrated
        "first_calibration_year": "1995",
        # split the regions of a region-by-region run into balanced batches, one HPC job per batch:
        # e.g. no_of_region_batches=20 and region_batch=0..19 (requires all_nuts3_regions_one_by_one=true)
        "no_of_region_batches": "1",
        "region_batch": "0",
        # nuts3: calibrate against regional yields (region by region or pooled)
        # national: one parameter set for Germany, objective = RMSE of the yearly weighted national means
        #   of simulations and observations over the regions having both values in a year (option B)
        "calibration_level": "nuts3",
        # share of each region's crop cells that is simulated (fixed random sample, at least
        # min_cells_per_region cells); "auto" = 0.1 for national, 1.0 (all cells) for nuts3
        "cell_sample_fraction": "auto",
        "min_cells_per_region": "3",
        "cell_sample_seed": "1",
        # parameter file; empty = calibratethese_<crop>_de.csv, otherwise calibratethese_<crop>.csv
        "params_file": "",
        # sceua: SCE-UA optimisation (repetitions, ngs)
        # grid: evaluate grid_points equally spaced values between low and high of each parameter
        #   (all combinations), e.g. 51 points = steps of 0.01 between 0.5 and 1.0
        "sampler": "sceua",
        "grid_points": "51",
        # cross-validation: years excluded from the objective, e.g. [2003]; their simulations are still
        # stored, output files get the suffix _loyo<years>
        "leave_out_years": "[]",
        "test_mode": "false",
        "all_nuts3_regions_one_by_one": False, ## Define on rundeck ##
        "only_nuts3_region_ids": "[]",  ## Define on rundeck ##
    }

    common.update_config(config, sys.argv, print_config=True, allow_new_keys=False)

    path_to_out_folder = config['path_to_out']
    if not os.path.exists(path_to_out_folder):
        try:
            os.makedirs(path_to_out_folder)
        except OSError:
            print("run-calibration.py: Couldn't create dir:", path_to_out_folder, "!")
    path_to_out_file = path_to_out_folder + "/run-calibration.out"
    with open(path_to_out_file, "a") as _:
        _.write(f"{datetime.now()} config: {config}\n")

    # read and validate all inputs before starting any subprocess, so that input errors fail fast
    # and don't leave channels/producer/consumer running on the HPC
    setups = monica_run_lib.read_sim_setups(config["setups-file"])
    run_setups = json.loads(config["run-setups"])
    if len(run_setups) < 1:
        print("run-calibration_de.py: run-setups is empty, nothing to calibrate.")
        return
    setup_id = run_setups[0]
    setup = setups[setup_id]

    # the crop is taken from the setup's crop-id (e.g. "WW_sfix_hauto" -> "WW"), so the input files
    # always match the crop that the producer simulates
    crop = setup["crop-id"].split("_")[0]
    path_to_obs_csv = f"calibration_{crop}_de.csv"
    path_to_weights_csv = f"Weights_{crop}_de.csv"
    # parameter file: as configured, otherwise prefer a German parameter file if there is one
    if config["params_file"]:
        path_to_params_csv = config["params_file"] if os.path.exists(config["params_file"]) else None
        params_file_candidates = config["params_file"]
    else:
        path_to_params_csv = next((p for p in [f"calibratethese_{crop}_de.csv", f"calibratethese_{crop}.csv"]
                                   if os.path.exists(p)), None)
        params_file_candidates = f"calibratethese_{crop}_de.csv or calibratethese_{crop}.csv"
    missing_files = [p for p in [path_to_obs_csv, path_to_weights_csv] if not os.path.exists(p)]
    if path_to_params_csv is None:
        missing_files.append(params_file_candidates)

    calibration_level = config["calibration_level"]
    if calibration_level not in ["nuts3", "national"]:
        raise ValueError(f"run-calibration_de.py: calibration_level must be nuts3 or national, got {calibration_level}")
    if calibration_level == "national" and config["all_nuts3_regions_one_by_one"]:
        raise ValueError("run-calibration_de.py: calibration_level=national requires all_nuts3_regions_one_by_one=false")
    cell_sample_fraction = (0.1 if calibration_level == "national" else 1.0) \
        if config["cell_sample_fraction"] == "auto" else float(config["cell_sample_fraction"])
    if not 0.0 < cell_sample_fraction <= 1.0:
        raise ValueError(f"run-calibration_de.py: cell_sample_fraction must be in (0, 1], got {cell_sample_fraction}")
    if config["sampler"] not in ["sceua", "grid"]:
        raise ValueError(f"run-calibration_de.py: sampler must be sceua or grid, got {config['sampler']}")
    if config["sampler"] == "grid" and int(config["grid_points"]) < 2:
        raise ValueError("run-calibration_de.py: grid_points must be at least 2")
    if missing_files:
        raise FileNotFoundError(f"run-calibration_de.py: missing input files for crop {crop}: {missing_files}")
    with open(path_to_out_file, "a") as _:
        _.write(f"{datetime.now()} crop: {crop}, observations: {path_to_obs_csv}, "
                f"parameters: {path_to_params_csv}, weights: {path_to_weights_csv}, "
                f"calibration level: {calibration_level}, cell sample fraction: {cell_sample_fraction}\n")

    # read observations (t/ha -> kg/ha, negative values = missing); years are taken from the header
    observations = []
    nuts3_region_id_to_name = {}
    with open(path_to_obs_csv) as file:
        dialect = csv.Sniffer().sniff(file.read(), delimiters=';,\t')
        file.seek(0)
        reader = csv.reader(file, dialect)
        years = [int(y) for y in next(reader)[2:]]
        for row in reader:
            id = int(row[1])
            name = row[0].strip()
            nuts3_region_id_to_name[id] = name
            for year, value in zip(years, row[2:]):
                yield_t = float(value)
                observations.append({
                    "id": id,
                    "year": year,
                    "value": np.nan if yield_t < 0.0 else yield_t * 1000.0
                })

    # spin-up: only observations from first_calibration_year on are used for calibration
    first_calibration_year = int(config["first_calibration_year"])
    observations = [o for o in observations if o["year"] >= first_calibration_year]

    # cross-validation: left-out years stay in the observation list (so their simulations are stored in
    # the results file), but with value NaN they are ignored by the objective
    leave_out_years = sorted(json.loads(config["leave_out_years"]))
    unknown_years = [y for y in leave_out_years if y < first_calibration_year or y > max(years)]
    if unknown_years:
        raise ValueError(f"run-calibration_de.py: leave_out_years {unknown_years} outside of the calibration "
                         f"years {first_calibration_year}-{max(years)}")
    for o in observations:
        o["observed"] = o["value"]
        o["held_out"] = o["year"] in leave_out_years
        if o["held_out"]:
            o["value"] = np.nan
    loyo_suffix = ("_loyo" + "-".join(map(str, leave_out_years))) if leave_out_years else ""

    with open(path_to_out_file, "a") as _:
        _.write(f"{datetime.now()} calibrating years {first_calibration_year}-{max(years)} "
                f"(simulation start {setup['start_date']}), left out: {leave_out_years}\n")

    # order obs list by id to avoid mismatch between observation/evaluation lists
    observations.sort(key=lambda r: [r["id"], r["year"]])

    # read parameters which are to be calibrated
    params = []
    with open(path_to_params_csv) as params_csv:
        dialect = csv.Sniffer().sniff(params_csv.read(), delimiters=';,\t')
        params_csv.seek(0)
        reader = csv.reader(params_csv, dialect)
        next(reader, None)  # skip the header
        for row in reader:
            p = {"name": row[0]}
            if len(row[1]) > 0:
                p["array"] = int(row[1])
            for n, i in [("low", 2), ("high", 3), ("step", 4), ("optguess", 5), ("minbound", 6), ("maxbound", 7)]:
                if len(row[i]) > 0:
                    p[n] = float(row[i])
            if len(row) == 9 and len(row[8]) > 0:
                # bind the expression now (default argument); a plain closure over row would
                # make every derived parameter use the expression of the last csv row
                p["derive_function"] = lambda _, _2, expr=row[8]: eval(expr)
            params.append(p)

    # read weights and number of crop pixels per region
    weights = {}
    region_pixels = {}
    with open(path_to_weights_csv) as weights_csv:
        dialect = csv.Sniffer().sniff(weights_csv.read(), delimiters=';,\t')
        weights_csv.seek(0)
        reader = csv.reader(weights_csv, dialect)
        next(reader, None)  # skip the header
        for row in reader:
            region_pixels[int(row[2])] = int(float(row[3]))
    # weight = share of the region's crop pixels, computed from the pixel counts: the weight column of the
    # file is rounded to 4 decimals, which is 0.0000 for small regions (e.g. 28 WW regions)
    total_pixels = sum(region_pixels.values())
    weights = {id_: px / total_pixels for id_, px in region_pixels.items()}
    #print("weights:", weights)

    # configure MONICA setup for spotpy
    only_nuts3_region_ids = json.loads(config["only_nuts3_region_ids"])

    # every requested region must have observations and a weight, otherwise the run would
    # either evaluate nothing or crash late with a KeyError
    selected_region_ids = only_nuts3_region_ids if len(only_nuts3_region_ids) > 0 \
        else sorted(nuts3_region_id_to_name.keys())
    missing_obs = sorted(set(selected_region_ids) - set(nuts3_region_id_to_name.keys()))
    missing_weights = sorted(set(selected_region_ids) - set(weights.keys()))
    if missing_obs or missing_weights:
        raise ValueError(f"run-calibration_de.py: regions without observations: {missing_obs}, "
                         f"regions without weights: {missing_weights}")

    no_of_region_batches = int(config["no_of_region_batches"])
    region_batch = int(config["region_batch"])
    if not 0 <= region_batch < no_of_region_batches:
        raise ValueError(f"run-calibration_de.py: region_batch must be in 0..{no_of_region_batches - 1}, "
                         f"got {region_batch}")
    if no_of_region_batches > 1 and not config["all_nuts3_regions_one_by_one"]:
        raise ValueError("run-calibration_de.py: no_of_region_batches > 1 requires all_nuts3_regions_one_by_one=true")

    to_be_run_only_nuts3_region_ids = []
    if config["all_nuts3_regions_one_by_one"]:
        region_ids = sorted(selected_region_ids)
        # regions without any pixel of this crop have nothing to simulate
        no_crop_ids = [id_ for id_ in region_ids if region_pixels[id_] == 0]
        region_ids = [id_ for id_ in region_ids if region_pixels[id_] > 0]
        if no_of_region_batches > 1:
            region_ids = balanced_region_batches(region_ids, region_pixels, no_of_region_batches)[region_batch]
        to_be_run_only_nuts3_region_ids = [[id_] for id_ in region_ids]
        with open(path_to_out_file, "a") as _:
            _.write(f"{datetime.now()} region batch {region_batch} of {no_of_region_batches}: "
                    f"{len(region_ids)} regions, {sum(region_pixels[i] for i in region_ids)} crop pixels: "
                    f"{region_ids}; regions without {crop} pixels (not run): {no_crop_ids}\n")
    else:
        to_be_run_only_nuts3_region_ids = [only_nuts3_region_ids]

    # all inputs are fine -> start channels, producer and consumer
    procs = []
    # make sure channels, producer and consumer are stopped even if the calibration raises an error
    atexit.register(lambda: [proc.terminate() for proc in procs if proc.poll() is None])

    prod_chan_data = get_reader_writer_srs_from_channel(config["path_to_channel"], "prod_chan")
    procs.append(prod_chan_data["chan"])
    cons_chan_data = get_reader_writer_srs_from_channel(config["path_to_channel"], "cons_chan")
    procs.append(cons_chan_data["chan"])

    procs.append(sp.Popen([
        config["path_to_python"],
        "run-producer_calibration_de.py",
        "mode=mbm-local-remote" if local_run else "mode=hpc-local-remote",
        f"server={config['server']}",
        f"port={config['prod-port']}",
        f"setups-file={config['setups-file']}",
        f"run-setups={config['run-setups']}",
        f"reader_sr={prod_chan_data['reader_sr']}",
        f"test_mode={config['test_mode']}",
        f"path_to_out={config['path_to_out']}",
        f"cell_sample_fraction={cell_sample_fraction}",
        f"min_cells_per_region={config['min_cells_per_region']}",
        f"cell_sample_seed={config['cell_sample_seed']}",
    ]))

    procs.append(sp.Popen([
        config["path_to_python"],
        "run-consumer_calibration_de.py",
        # "mode=remoteConsumer-remoteMonica",
        f"server={config['server']}",
        f"port={config['cons-port']}",
        f"run-setups={config['run-setups']}",
        f"writer_sr={cons_chan_data['writer_sr']}",
        f"path_to_out={config['path_to_out']}",
    ]))

    con_man = common.ConnectionManager()
    cons_reader = con_man.try_connect(cons_chan_data["reader_sr"], cast_as=fbp_capnp.Channel.Reader, retry_secs=1)
    prod_writer = con_man.try_connect(prod_chan_data["writer_sr"], cast_as=fbp_capnp.Channel.Writer, retry_secs=1)

    spot_setup = None
    for current_only_nuts3_region_ids in to_be_run_only_nuts3_region_ids:
        # start timer 
        start_time = time.time()

        # output file prefix: region id(s), or DE if all regions are calibrated together, plus left-out years
        nuts3_region_folder_name = ("-".join(map(str, current_only_nuts3_region_ids)) or "DE") + loyo_suffix
        # regions without any pixel of the crop can't be simulated
        filtered_observations = [o for o in observations if region_pixels[o["id"]] > 0]
        if len(current_only_nuts3_region_ids) > 0:
            filtered_observations = [o for o in filtered_observations if o["id"] in current_only_nuts3_region_ids]
        if len(filtered_observations) == 0:
            continue
        # regions without any observation in the calibration years (or only held-out ones) have nothing to
        # calibrate against; the objective would fail and stop all following regions of the job
        if all(np.isnan(o["value"]) for o in filtered_observations):
            if len(to_be_run_only_nuts3_region_ids) > 1:
                with open(path_to_out_file, "a") as _:
                    _.write(f"{datetime.now()} skipped region(s) {current_only_nuts3_region_ids}: no observations "
                            f"in the calibration years\n")
                continue
            raise ValueError(f"run-calibration_de.py: no observations in the calibration years for regions "
                             f"{current_only_nuts3_region_ids or 'DE'}")
        if spot_setup:
            del spot_setup

        # weight of every observation = crop pixel share of its region
        # nuts3: weighted RMSE over all region-years (= plain RMSE for a single region)
        # national: RMSE of the yearly weighted national means (see calibration_level)
        spot_setup = calibration_spotpy_setup_MONICA.spot_setup(
            params, filtered_observations, prod_writer, cons_reader, path_to_out_folder,
            current_only_nuts3_region_ids,
            objective="national_mean" if calibration_level == "national" else "rmse",
            obs_weights=[weights[o["id"]] for o in filtered_observations])

        # which observation (region, year, observed value, weight, held out) belongs to which
        # simulation_<i> column of the results file; used by the cross-validation evaluation
        write_observations_file(f"{path_to_out_folder}/{nuts3_region_folder_name}_observations.csv",
                                filtered_observations, weights, calibration_level)

        if config["sampler"] == "grid":
            try:
                grid_search(spot_setup, params, int(config["grid_points"]),
                            f"{path_to_out_folder}/{nuts3_region_folder_name}_GRID_monica_results.csv",
                            f"{path_to_out_folder}/{nuts3_region_folder_name}_best.out",
                            f"{path_to_out_folder}/{nuts3_region_folder_name}_GRID_objectivefunction_MONICA.png")
            except calibration_spotpy_setup_MONICA.NoSimulationResultsError as e:
                if len(to_be_run_only_nuts3_region_ids) > 1:
                    with open(path_to_out_file, "a") as _:
                        _.write(f"{datetime.now()} skipped region(s) {current_only_nuts3_region_ids}: {e}\n")
                    continue
                raise
            with open(path_to_out_folder + "/spot_setup.out", "a") as _:
                _.write(f"{datetime.now()} grid search for {nuts3_region_folder_name} took "
                        f"{time.time() - start_time:.2f} seconds\n")
            continue

        rep = int(config["repetitions"]) #initial number was 10
        results = []
        #Set up the sampler with the model above
        sampler = spotpy.algorithms.sceua(spot_setup, dbname=f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_monica_results", dbformat="csv")
        # sampler = spotpy.algorithms.dream(spot_setup, dbname=f"{path_to_out_folder}/{nuts3_region_folder_name}_DREAM_monica_results", dbformat="csv")
        #Run the sampler to produce the paranmeter distribution
        #and identify optimal parameters based on objective function
        #ngs = number of complexes
        #kstop = max number of evolution loops before convergence
        #peps = convergence criterion
        #pcento = percent change allowed in kstop loops before convergence
        with open(path_to_out_folder + "/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} sampler starts run-cal\n")

        try:
            sampler.sample(rep, ngs=int(config["ngs"]), kstop=100, peps=0.0001, pcento=0.0001)
        except calibration_spotpy_setup_MONICA.NoSimulationResultsError as e:
            # in region-by-region runs, skip regions without any cell of the crop mask instead of
            # stopping the whole run; a single selected region (or the national set) still fails
            if len(to_be_run_only_nuts3_region_ids) > 1:
                with open(path_to_out_file, "a") as _:
                    _.write(f"{datetime.now()} skipped region(s) {current_only_nuts3_region_ids}: {e}\n")
                continue
            raise


        # sampler.sample(rep, nChains = 20, nCr = 3, eps = (10e-6), convergence_limit=1.0)

        with open(path_to_out_folder + "/spot_setup.out", "a") as _:
            _.write(f"{datetime.now()} sampler ends run-cal\n")
        # end timer
        end_time = time.time()
        time_taken = end_time - start_time
        if time_taken > 10:
            with open(path_to_out_folder + "/spot_setup.out", "a") as _:
                _.write(f"{datetime.now()} Time taken to calibrate: {time_taken:.2f} seconds\n")
            #print(f"Time taken to calibrate: {time_taken:.2f} seconds")


        def print_status_final(self, stream):
            print("\n*** Final SPOTPY summary ***")
            print(
                "Total Duration: "
                + str(round((time.time() - self.starttime), 2))
                + " seconds"
            , file=stream)
            print("Total Repetitions:", self.rep, file=stream)

            if self.optimization_direction == "minimize":
                print("Minimal objective value: %g" % (self.objectivefunction_min), file=stream)
                print("Corresponding parameter setting:", file=stream)
                for i in range(self.parameters):
                    text = "%s: %g" % (self.parnames[i], self.params_min[i])
                    print(text, file=stream)

            if self.optimization_direction == "maximize":
                print("Maximal objective value: %g" % (self.objectivefunction_max), file=stream)
                print("Corresponding parameter setting:", file=stream)
                for i in range(self.parameters):
                    text = "%s: %g" % (self.parnames[i], self.params_max[i])
                    print(text, file=stream)

            if self.optimization_direction == "grid":
                print("Minimal objective value: %g" % (self.objectivefunction_min), file=stream)
                print("Corresponding parameter setting:", file=stream)
                for i in range(self.parameters):
                    text = "%s: %g" % (self.parnames[i], self.params_min[i])
                    print(text, file=stream)

                print("Maximal objective value: %g" % (self.objectivefunction_max), file=stream)
                print("Corresponding parameter setting:", file=stream)
                for i in range(self.parameters):
                    text = "%s: %g" % (self.parnames[i], self.params_max[i])
                    print(text, file=stream)

            print("******************************\n", file=stream)


        path_to_best_out_file = f"{path_to_out_folder}/{nuts3_region_folder_name}_best.out"
        with open(path_to_best_out_file, "a") as _:
            print_status_final(sampler.status, _)

        #with open(path_to_out_folder + "/spot_setup.out", "a") as _:
        #    _.write(f"{datetime.now()} results written run-cal\n\n")

        #Extract the parameter samples from distribution
        results = spotpy.analyser.load_csv_results(f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_monica_results")

        # Plot how the objective function was minimized during sampling
        #font = {"family": "calibri",
        #        "weight": "normal",
        #        "size": 18}
        fig = plt.figure(1, figsize=(9, 6))
        #plt.plot(results["like1"],  marker='o')
        plt.plot(results["like1"], "r+")
        plt.show()
        plt.ylabel("RMSE (kg/ha)")
        plt.xlabel("Iteration")
        fig.savefig(f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_objectivefunctiontrace_MONICA.png", dpi=150)
        plt.close(fig)

        # OW addition
        #df = pd.read_csv (f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_monica_results.csv")
        #columns_of_interest = ['like1','parSpecificLeafArea', 'parMaxAssimilationRate', 'parDaylengthRequirement', 'parBaseDaylength', 'parCropSpecificMaxRootingDepth']
        #df_selected = df[columns_of_interest]
        #lowest_like1_values = df_selected.nsmallest(100, 'like1')['like1']
        #df_lowest_like1 = df_selected[df_selected['like1'].isin(lowest_like1_values)]

        # Drop the 'like1' column from the DataFrame as it's no longer needed for plotting
        #df_lowest_like1 = df_lowest_like1.drop(columns=['like1'])

        #Drop any non-numeric columns (like 'chain') before creating the pair plot
        #df_lowest_like1_numeric = df_lowest_like1.select_dtypes(include='number')
        #fig1 = sns.pairplot(df_lowest_like1_numeric)
        #fig1.savefig(f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_pair_MONICA.png", dpi=150)
        #plt.close(fig1.fig)



        # Plot the percentage differences
        #fig = plt.figure(1, figsize=(9, 6))
        #plt.plot(results["like1"], "r+")
        #plt.ylabel("Percentage Difference (%)")
        #plt.xlabel("Iteration")
        #plt.show()

        # Save the plot
        #fig.savefig(f"{path_to_out_folder}/{nuts3_region_folder_name}_SCEUA_percentage_difference_MONICA.png", dpi=150)
        #plt.close(fig)


        del results
    # kill the two channels and the producer and consumer
    for proc in procs:
        proc.terminate()

    
    print("sampler_MONICA.py finished")

if __name__ == "__main__":
    run_calibration()


