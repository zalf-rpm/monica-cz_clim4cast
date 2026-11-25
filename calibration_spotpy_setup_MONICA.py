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

import capnp
from datetime import datetime
import json
import logging
import os
from pathlib import Path

import numpy as np
import spotpy
import re

PATH_TO_REPO = Path(os.path.realpath(__file__)).parent
PATH_TO_MAS_INFRASTRUCTURE_REPO = PATH_TO_REPO / "../mas-infrastructure" #/nb(*T)
PATH_TO_CAPNP_SCHEMAS = (PATH_TO_MAS_INFRASTRUCTURE_REPO / "capnproto_schemas").resolve()
abs_imports = [str(PATH_TO_CAPNP_SCHEMAS)]
fbp_capnp = capnp.load(str(PATH_TO_CAPNP_SCHEMAS / "fbp.capnp"), imports=abs_imports)

class NoSimulationResultsError(RuntimeError):
    """MONICA returned no results at all, e.g. because the selected region has no cells of the crop mask."""


class spot_setup(object):
    """
    objective:
        "rmse"          weighted RMSE over all obs/sim pairs (obs_weights, default equal weights; like1 in kg/ha)
        "national_mean" RMSE over the years of the weighted means of obs and sim per year, using only the
                        regions with both values in that year (obs_weights = region weight per observation)
    weight_per_region is kept for compatibility with run-calibration.py (CZ) and not used: a single
    constant weight scales the objective without changing the optimum.
    """
    def __init__(self, user_params, observations, prod_writer, cons_reader, path_to_out, only_nuts3_region_ids,
                 weight_per_region=None, objective="rmse", obs_weights=None):
        self.user_params = user_params
        self.params = []
        self.observations = observations
        self.obs_flat_list = list(map(lambda d: d["value"], observations))
        self.obs_years = [d["year"] for d in observations]
        self.prod_writer = prod_writer
        self.cons_reader = cons_reader
        self.path_to_out_file = path_to_out + "/spot_setup.out"
        self.only_nuts3_region_ids = only_nuts3_region_ids
        self.weight_per_region = weight_per_region
        if objective not in ["rmse", "national_mean"]:
            raise ValueError(f"unknown objective {objective}")
        if obs_weights is not None and len(obs_weights) != len(observations):
            raise ValueError("obs_weights must have one weight per observation")
        if objective == "national_mean" and obs_weights is None:
            raise ValueError("objective national_mean needs obs_weights")
        self.objective = objective
        self.obs_weights = obs_weights

        if not os.path.exists(path_to_out):
            try:
                os.makedirs(path_to_out)
            except OSError:
                print("spot_setup.__init__: Couldn't create dir:", path_to_out, "!")

        with open(self.path_to_out_file, "a") as _:
            _.write(f"observations: {self.observations}\n")
            _.write(f"obs_flat_list: {self.obs_flat_list}\n")

        for par in user_params:
            par_name = par["name"]
            if "array" in par:
                par["name"] = f"{par_name}_{par['array']}"  # spotpy does not allow two parameters to have the same name
                del par["array"]
            if "derive_function" not in par:  # spotpy does not care about derived params
                self.params.append(spotpy.parameter.Uniform(**par))

    def parameters(self):
        return spotpy.parameter.generate(self.params)

    def simulation(self, vector):
        # vector = MaxAssimilationRate, AssimilateReallocation, RootPenetrationRate
        msg_content = dict(zip(vector.name, vector))
        msg_content["only_nuts3_region_ids"] = self.only_nuts3_region_ids
        out_ip = fbp_capnp.IP.new_message(content=json.dumps(msg_content))
        self.prod_writer.write(value=out_ip).wait()
        with open(self.path_to_out_file, "a") as _:
            _.write(f"{datetime.now()} sent params to monica setup: {vector}\n")
        print("sent params to monica setup:", vector, flush=True)

        msg = self.cons_reader.read().wait()
        # check for end of data from in port
        if msg.which() == "done":
            return

        in_ip = msg.value.as_struct(fbp_capnp.IP)
        s: str = in_ip.content.as_text()
        nuts3_region_id_and_year_to_avg_yield = json.loads(s)
        if len(nuts3_region_id_and_year_to_avg_yield) == 0:
            # nothing was simulated (e.g. region ids unknown to the producer's region grid);
            # stop instead of silently calibrating against nothing
            raise NoSimulationResultsError(f"No simulation results received for regions {self.only_nuts3_region_ids}. "
                                           f"Check producer.out (number of sent envs) and the region ids.")

        #with open(self.path_to_out_file, "a") as _:
        #    _.write(f"{datetime.now()} jsons loaded cal-sp-set-M\n")
        # print("received monica results:", country_id_and_year_to_avg_yield, flush=True)

        # remove all simulation results which are not in the observed list
        sim_list = []
        for d in self.observations:
            key = f"{d['id']}|{d['year']}"
            if key in nuts3_region_id_and_year_to_avg_yield:
                #if np.isnan(d["value"]):
                #    sim_list.append(np.nan)
                #else:
                sim_list.append(nuts3_region_id_and_year_to_avg_yield[key])
            else:
                sim_list.append(np.nan)

        #with open(self.path_to_out_file, "a") as _:
        #    _.write(f"{datetime.now()} simulation and observation matchedcal-sp-set-M\n\n")

        print("len(sim_list):", len(sim_list), "== len(self.obs_list):", len(self.obs_flat_list), flush=True)
        with open(self.path_to_out_file, "a") as _:
            #_.write(f"received monica results: {country_id_and_year_to_avg_yield}\n")
            _.write(f"{datetime.now()}  len(sim_list): {len(sim_list)} == len(self.obs_list): {len(self.obs_flat_list)}\n")
            #_.write(f"sim_list: {sim_list}\n")
            #_.write(f"obs_list: {self.obs_flat_list}\n")
        # besides the order the length of observation results and simulation results should be the same
        assert len(sim_list) == len(self.obs_flat_list)
        return sim_list if len(sim_list) > 0 else None

    def evaluation(self):
        return self.obs_flat_list


    def objectivefunction(self, simulation, evaluation):

        #return unbiased_rmse_RB(evaluation, simulation)
        #return spotpy.objectivefunctions.rmse(evaluation, simulation)
        #return calculate_percentage_difference_new(evaluation, simulation)
        if self.objective == "national_mean":
            return calculate_national_mean_rmse(evaluation, simulation, self.obs_years, self.obs_weights)
        return calculate_weighted_rmse(evaluation, simulation, self.obs_weights)



def calculate_weighted_rmse(evaluation, simulation, weights=None):
    """
    Weighted RMSE over all pairs where both observation and simulation are available.

    .. math::

        RMSE_weighted = \\sqrt{\\frac{\\sum_{i} w_i (sim_i - obs_i)^2}{\\sum_{i} w_i}}

    With weights=None all pairs get the same weight, i.e. the plain RMSE (same unit as the yields).
    Pairs where obs or sim is NaN (missing or left-out observations) are excluded from the
    numerator and the denominator.

    :param evaluation: observed values
    :param simulation: simulated values, same length and order as evaluation
    :param weights: optional weight per pair (e.g. region crop-area fraction), same length
    :return: weighted RMSE
    :raises RuntimeError: if there is no valid obs/sim pair (nothing to calibrate against)
    """
    obs = np.asarray(evaluation, dtype=float)
    sim = np.asarray(simulation, dtype=float)
    w = np.ones_like(obs) if weights is None else np.asarray(weights, dtype=float)
    if not (len(obs) == len(sim) == len(w)):
        raise ValueError(f"evaluation ({len(obs)}), simulation ({len(sim)}) and weights ({len(w)}) "
                         f"must have the same length")

    valid = ~np.isnan(obs) & ~np.isnan(sim) & ~np.isnan(w)
    if not valid.any() or w[valid].sum() <= 0:
        raise RuntimeError("No valid observation/simulation pair: simulated region-years do not match "
                           "any observed region-year.")
    return np.sqrt(np.sum(w[valid] * (sim[valid] - obs[valid]) ** 2) / np.sum(w[valid]))



def national_means(evaluation, simulation, years, weights):
    """
    Weighted means of observations and simulations per year, using only the regions (pairs) where both
    values are available in that year, so that both means are built from the same regions.

    :return: list of (year, obs_mean, sim_mean) for the years with at least one valid pair
    """
    obs = np.asarray(evaluation, dtype=float)
    sim = np.asarray(simulation, dtype=float)
    w = np.asarray(weights, dtype=float)
    years = np.asarray(years)
    if not (len(obs) == len(sim) == len(w) == len(years)):
        raise ValueError("evaluation, simulation, years and weights must have the same length")
    valid = ~np.isnan(obs) & ~np.isnan(sim) & ~np.isnan(w) & (w > 0)
    means = []
    for year in np.unique(years):
        sel = valid & (years == year)
        if sel.any():
            means.append((int(year), np.sum(w[sel] * obs[sel]) / np.sum(w[sel]),
                          np.sum(w[sel] * sim[sel]) / np.sum(w[sel])))
    return means


def calculate_national_mean_rmse(evaluation, simulation, years, weights):
    """RMSE over the years between the yearly weighted national means of observations and simulations."""
    means = national_means(evaluation, simulation, years, weights)
    if not means:
        raise RuntimeError("No year with a valid observation/simulation pair.")
    return float(np.sqrt(np.mean([(s - o) ** 2 for _, o, s in means])))


def calculate_percentage_difference_new(evaluation: list, simulation: list) -> float:
    """
    Calculate the mean absolute percentage difference between observed (evaluation) and simulated values.

    :param evaluation: Observed data to compare with simulation data.
    :type evaluation: list of numeric values

    :param simulation: Simulation data to compare with evaluation data.
    :type simulation: list of numeric values

    :return: Mean absolute percentage difference, or np.nan if division by zero occurs or
             if the lengths do not match.
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        percentage_differences = []
        for eval_val, sim_val in zip(evaluation, simulation):
            if eval_val != 0:
                percentage_difference = abs((eval_val - sim_val) / eval_val * 100)
                percentage_differences.append(percentage_difference)
            else:
                percentage_differences.append(np.nan)  # Handle division by zero
        
        # Convert list to numpy array to handle np.nan and compute the mean
        percentage_differences = np.array(percentage_differences)
        
        # Return the mean, ignoring nan values
        return np.nanmean(percentage_differences)
    else:
        logging.warning("evaluation and simulation lists do not have the same length or are empty.")
        return np.nan



def rmse(evaluation, simulation):
    """
    Root Mean Squared Error

        .. math::

         RMSE=\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Root Mean Squared Error
    :rtype: float
    """
    if len(evaluation) == len(simulation) > 0:
        return np.sqrt(mse(evaluation, simulation))
    else:
        logging.warning("evaluation and simulation lists do not have the same length.")
        return np.nan



def mse(evaluation, simulation):
    """
    Mean Squared Error

        .. math::

         MSE=\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2

    :evaluation: Observed data to compared with simulation data.
    :type: list

    :simulation: simulation data to compared with evaluation data
    :type: list

    :return: Mean Squared Error
    :rtype: float
    """

    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mse = np.nanmean((obs - sim) ** 2)
        return mse
    else:
        #logging.warning(
        #    "evaluation and simulation lists does not have the same length."
        #)
        return np.nan


def calculate_mbe(evaluation, simulation):
    """
    Calculate the Mean Bias Error.

    Args:
    y_true (array-like): True values.
    y_pred (array-like): Predicted values.

    Returns:
    float: Mean Bias Error.
    """
    if len(evaluation) == len(simulation):
        obs, sim = np.array(evaluation), np.array(simulation)
        mbe = np.nanmean(sim-obs)#swapped obs-sim to sim-obs
    return mbe


#RB addition

def unbiased_rmse_RB(evaluation, simulation):
    """
    Calculate RMSE with prior bias correction.

    Args:
    y_true (array-like): True values.
    y_pred (array-like): Predicted values.

    Returns:
    float: RMSE after bias correction.
    """

    mbe = calculate_mbe(evaluation, simulation)
    print("      mean bias estimate        =",mbe)
    #subtract mean bias estimate from simulation. 
    #This will give bias corrected simulation.
    y_pred_adjusted = simulation - mbe
    return rmse(evaluation, y_pred_adjusted)
