# German MONICA calibration – step-by-step guide

The calibration has three steps:

1. **Experimental sites** – own scripts (not part of this guide). Results go into `crop_calibration_de.json`.
2. **National** – 4 parameters, one set for Germany, calibrated against yearly weighted national mean yields.
3. **Regional** – `FieldConditionModifier` per NUTS3 region, by grid search.

Results are passed from one step to the next **by hand** via `crop_calibration_de.json`.
Both the national and the regional step are cross-validated with leave-one-year-out (LOYO).

The examples use winter wheat (`run-setups=[3]`, crop `WW`). For other crops, use the matching setup id
from `sim_setups_calibration_de.csv` and the crop's files.

---

## Before you start

- **Every Rundeck job checks out the repository when it starts.** Changes (code, `crop_calibration_de.json`,
  parameter files) reach the HPC only after `git push`, and only for jobs started afterwards.
- Each Rundeck job starts its own MONICA instance, so jobs can run in parallel.
- Output of a job: `/beegfs/rpm/projects/monica/out/palka_<job-id>_<date>/`.
  Checkout of a job: `/beegfs/rpm/projects/monica/projectcheckouts/palka_<job-id>_<date>/monica-cz_clim4cast/`.
- Input files per crop: `calibration_<crop>_de.csv` (observations, dry matter, t/ha),
  `Weights_<crop>_de.csv` (crop pixels per region), `calibratethese_<crop>.csv` (national parameters).
  `calibratethese_WB.csv` and `calibratethese_PO.csv` do not exist yet.

### Running the evaluation scripts on the HPC

`collect-calibration-results_de.py` and `evaluate-cv_de.py` need only the Python standard library
(Python 3.6 or newer). They only **read** the job output folders; their own output goes to the paths
given with `--out` / `--summary`, e.g. in `/beegfs/palka/`.

One-time setup on the HPC (own copy of the scripts):
```
cd /beegfs/palka
git clone https://github.com/zalf-rpm/monica-cz_clim4cast.git
mkdir -p /beegfs/palka/calibration_test
```
Before each use, update the copy:
```
cd /beegfs/palka/monica-cz_clim4cast && git pull
```
Run the scripts from there with `python3` (check `python3 --version`; if it is older than 3.6, use
`/home/rpm/.conda/envs/clim4cast/bin/python` instead):
```
cd /beegfs/palka/monica-cz_clim4cast
python3 <script> <output folders> ... --out /beegfs/palka/calibration_test/<file>.csv
```

List the job output folders with `ls -d /beegfs/rpm/projects/monica/out/palka_*`. Output folders can be
given with wildcards, e.g. `"/beegfs/rpm/projects/monica/out/palka_2216*"`, but make sure a pattern does
not match folders of other (old) runs.

---

## Test chain (WW, small settings)

Test settings: few repetitions, few cells, coarse grid. Durations are rough estimates.

### Step 0 – Preparation (local)

1. The experimental results for WW are in `crop_calibration_de.json`.
2. Create `calibratethese_WW_regional.csv` (only `FieldConditionModifier`, range 0.5–1.0):
   ```
   name,position in array,low,high,step,optguess,minbound,maxbound,derived params
   FieldConditionModifier,,0.5,1.0,,,,,
   ```
3. Commit and push.

### Step 1 – National calibration on all years (1 job, ~5 min)

Rundeck parameters:
```
calibration_level=national run-setups=[3] repetitions=10 cell_sample_fraction=0.02
```
Check in the output folder:
- `run-calibration.out`: `calibration level: national, cell sample fraction: 0.02`,
  `calibrating years 1995-2024 ..., left out: []`
- `producer.out`: `collected ... cells in 380 regions`, `sampled ... cells (fraction 0.02, ...)`
- `consumer.out`: one line per iteration with `0 with MONICA errors`
- `DE_SCEUA_monica_results.csv`, `DE_best.out`, `DE_observations.csv`, `DE_SCEUA_objectivefunctiontrace_MONICA.png`

### Step 2 – National cross-validation folds (2 jobs for the test, ~5 min each)

One job per left-out year; for the test e.g. 2003 and 2018:
```
calibration_level=national run-setups=[3] repetitions=10 cell_sample_fraction=0.02 leave_out_years=[2003]
calibration_level=national run-setups=[3] repetitions=10 cell_sample_fraction=0.02 leave_out_years=[2018]
```
Check: `run-calibration.out` shows `left out: [2003]`; files are named `DE_loyo2003_...`.

### Step 3 – Evaluate the national results (HPC, see "Running the evaluation scripts")

Best national parameters (from step 1):
```
python3 collect-calibration-results_de.py /beegfs/rpm/projects/monica/out/palka_<id-step1>_<date> --out /beegfs/palka/calibration_test/national_best_WW.csv
```
National cross-validation (from step 2):
```
python3 evaluate-cv_de.py /beegfs/rpm/projects/monica/out/palka_<id-fold2003>_<date> /beegfs/rpm/projects/monica/out/palka_<id-fold2018>_<date> \
        --out /beegfs/palka/calibration_test/national_cv_WW.csv --summary /beegfs/palka/calibration_test/national_cv_summary_WW.csv
```
Expected: `national: 2 predictions, cross-validated RMSE ... kg/ha, bias ... kg/ha`.

### Step 4 – Transfer the national parameters (local)

Take the best parameters from step 3 (or from `DE_best.out`) and enter them for `WW` in
`crop_calibration_de.json`:

| Parameter | Where in `crop_calibration_de.json` | How |
|---|---|---|
| `SpecificLeafArea` | `cultivar` → `SpecificLeafArea` (6 values) | **multiply** all 6 values by the calibrated factor |
| `MaxAssimilationRate` | `cultivar` → `MaxAssimilationRate` | replace the value |
| `CropSpecificMaxRootingDepth` | `cultivar` → `CropSpecificMaxRootingDepth` | replace the value |
| `FieldConditionModifier` | `species` → `FieldConditionModifier` | replace the value (it is re-calibrated per region in step 3) |

For the test, commit with a clear message (e.g. `TEST: national WW parameters`) so it can be reverted
afterwards with `git revert <commit>`. Push.

### Step 5 – Regional grid search (1 job for the test, ~15 min)

One region batch with a coarse grid (6 values: 0.5, 0.6, ..., 1.0):
```
all_nuts3_regions_one_by_one=true no_of_region_batches=10 region_batch=0 run-setups=[3] sampler=grid grid_points=6 params_file=calibratethese_WW_regional.csv
```
Check:
- `run-calibration.out`: `parameters: calibratethese_WW_regional.csv`,
  `region batch 0 of 10: 37 regions, ...`
- per region: `<region>_GRID_monica_results.csv` (6 runs), `<region>_best.out`, `<region>_observations.csv`,
  `<region>_GRID_objectivefunction_MONICA.png` (RMSE against `FieldConditionModifier`)
- regions without simulation cells are listed as `skipped region(s)` in `run-calibration.out`

### Step 6 – Evaluate the regional results (HPC)

Best `FieldConditionModifier` per region (fitted on all years):
```
python3 collect-calibration-results_de.py /beegfs/rpm/projects/monica/out/palka_<id-step5>_<date> --obs calibration_WW_de.csv --out /beegfs/palka/calibration_test/regional_best_WW.csv
```
Regional cross-validation (all LOYO folds are computed from the grid results, no extra runs):
```
python3 evaluate-cv_de.py /beegfs/rpm/projects/monica/out/palka_<id-step5>_<date> --out /beegfs/palka/calibration_test/regional_cv_WW.csv --summary /beegfs/palka/calibration_test/regional_cv_summary_WW.csv
```
Expected: `nuts3: ... predictions, cross-validated RMSE ... kg/ha, ...`; the summary has one row per region.

### Step 7 – Clean up after the test

Revert the test parameters in `crop_calibration_de.json` (`git revert <commit of step 4>`) and push.

---

## Production settings

| Run | Rundeck parameters | Jobs | Rough duration (WW) |
|---|---|---|---|
| National, all years | `calibration_level=national run-setups=[3] repetitions=500` | 1 | ~13-14 h |
| National LOYO folds | `... repetitions=500 leave_out_years=[<year>]` for 1995 ... 2024 | 30 | ~13-14 h each |
| Regional grid (incl. LOYO) | `all_nuts3_regions_one_by_one=true no_of_region_batches=10 region_batch=<0..9> run-setups=[3] sampler=grid grid_points=<n> params_file=calibratethese_WW_regional.csv` | 10 | ~45-60 min each (20 points) |

- National runs use 10% of the cells by default (`cell_sample_fraction`).
- Regional range: national `FieldConditionModifier` -20% to +20%, rounded outwards to the grid step, set as
  `low`/`high` in `calibratethese_<crop>_regional.csv`. `grid_points = (high - low) / step + 1`, e.g. national 0.93,
  step 0.02: 0.74-1.12, `grid_points=20`. Without `grid_points`, 51 values are used.
- Regions without observations in the calibration years are skipped (listed in `run-calibration.out`).
- To re-run some regions (e.g. with a wider range): `only_nuts3_region_ids=[...]` with
  `all_nuts3_regions_one_by_one=true`, optionally split with `no_of_region_batches`. Give the old and the new
  folders to the evaluation scripts: for a region found more than once, the results file with the most runs is used.
- Lookup table for the production simulations (regional values, national value where there is no reliable
  regional one):
  ```
  python3 make-regional-fcm_de.py --best regional_best_WW.csv --obs calibration_WW_de.csv           --weights Weights_WW_de.csv --national 0.93 --out regional_fcm_WW_de.csv
  ```
  The national value is used for regions without a valid run, with fewer than `--min-obs` (default 2) observed
  years, or without crop pixels.
- All jobs of one crop must use the same `no_of_region_batches`.
- Run the regional step only after the national parameters of all years are in `crop_calibration_de.json`.

---

## Configuration reference (`run-calibration_de.py`)

| Key | Default | Meaning |
|---|---|---|
| `run-setups` | `[1]` | setup id in `sim_setups_calibration_de.csv`; determines the crop and its input files |
| `calibration_level` | `nuts3` | `nuts3`: regional yields; `national`: yearly weighted national means (all regions, one job) |
| `all_nuts3_regions_one_by_one` | `false` | calibrate each region separately (regional step) |
| `only_nuts3_region_ids` | `[]` | restrict to these regions (empty = all) |
| `no_of_region_batches`, `region_batch` | `1`, `0` | split the regions of a one-by-one run into balanced batches (one job per batch) |
| `params_file` | (empty) | parameter file; empty = `calibratethese_<crop>_de.csv`, else `calibratethese_<crop>.csv` |
| `sampler` | `sceua` | `sceua`: SCE-UA optimisation; `grid`: grid search |
| `repetitions`, `ngs` | `500`, `4` | SCE-UA runs and complexes (repetitions should be well above `ngs * (2 * no_of_params + 1)`) |
| `grid_points` | `51` | grid values per parameter between `low` and `high` of the parameter file |
| `first_calibration_year` | `1995` | years before are spin-up (simulations start 1990-01-01, see setup file) |
| `leave_out_years` | `[]` | years excluded from the objective (cross-validation); files get `_loyo<years>` |
| `cell_sample_fraction` | `auto` | share of each region's cells simulated; auto = 0.1 national, 1.0 regional |
| `min_cells_per_region`, `cell_sample_seed` | `3`, `1` | minimum cells per region and seed of the fixed cell sample |

## Output files per calibration (`<prefix>` = region id or `DE`, plus `_loyo<years>`)

| File | Content |
|---|---|
| `<prefix>_SCEUA_monica_results.csv` / `<prefix>_GRID_monica_results.csv` | one row per run: `like1` (RMSE, kg/ha), parameters, simulated yield per region-year |
| `<prefix>_observations.csv` | observation (region, year, observed yield, weight, held out) of each `simulation_<i>` column |
| `<prefix>_best.out` | best parameter set |
| `*.png` | objective function trace (SCE-UA) or RMSE against parameter value (grid) |
| `run-calibration.out`, `producer.out`, `consumer.out`, `spot_setup.out` | logs of the job |

## Notes

- The regional cross-validation uses national parameters calibrated on all years, including the held-out
  year, so it is slightly optimistic compared to a fully nested cross-validation.
- Region-specific `FieldConditionModifier` values from step 3 cannot be stored in `crop_calibration_de.json`
  (one value per crop); production simulations use the lookup table of `make-regional-fcm_de.py`.
