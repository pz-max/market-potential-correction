# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging

import numpy as np
import pandas as pd
import re

import xarray as xr
import geopandas as gpd
import powerplantmatching as pm
from powerplantmatching.export import map_country_bus

from vresutils.costdata import annuity
from vresutils.load import timeseries_opsd
from vresutils import transfer as vtransfer

import pypsa
from pypsa.linopf import (get_var, define_constraints, linexpr, join_exprs,
                          network_lopf, ilopf)

idx = pd.IndexSlice

from pathlib import Path
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    import logging

    kwargs = snakemake.config.get('logging', dict())
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath('..', 'logs', f"{snakemake.rule}.log")
        logfile = snakemake.log.get('python', snakemake.log[0] if snakemake.log
                                    else fallback_path)
        kwargs.update(
            {'handlers': [
                # Prefer the 'python' log, otherwise take the first log for each
                # Snakemake rule
                logging.FileHandler(logfile),
                logging.StreamHandler()
                ]
            })
    logging.basicConfig(**kwargs)


def load_costs(Nyears, tech_costs, config, elec_config):
    if tech_costs is None:
        tech_costs = snakemake.input.tech_costs

    if config is None:
        config = snakemake.config['costs']

    # set all asset costs and other parameters
    costs = pd.read_csv(tech_costs, index_col=list(range(3))).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"),"value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"),"value"] *= config['USD2013_to_EUR2013']

    costs = (costs.loc[idx[:,config['year'],:], "value"]
             .unstack(level=2).groupby("technology").sum(min_count=1))

    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                          "discount rate" : config['discountrate'],
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : 25})

    costs["capital_cost"] = ((annuity(costs["lifetime"], costs["discount rate"]) +
                             costs["FOM"]/100.) *
                             costs["investment"] * Nyears)

    costs.at["fuel cell","capital_cost"] = ((annuity(costs.at["fuel cell","lifetime"], costs.at["fuel cell","discount rate"]) +
                             costs.at["fuel cell","FOM"]/100.) *
                             costs.at["fuel cell","investment"] * curves.at[int(snakemake.wildcards.iterations), 'parameter'] * Nyears)

    costs.at["electrolysis","capital_cost"] = ((annuity(costs.at["electrolysis","lifetime"], costs.at["electrolysis","discount rate"]) +
                             costs.at["electrolysis","FOM"]/100.) *
                             costs.at["electrolysis","investment"] * curves.at[int(snakemake.wildcards.iterations), 'parameter'] * Nyears)

    costs.at['OCGT', 'fuel'] = costs.at['gas', 'fuel']
    costs.at['CCGT', 'fuel'] = costs.at['gas', 'fuel']

    costs['marginal_cost'] = costs['VOM'] + costs['fuel'] / costs['efficiency']

    costs = costs.rename(columns={"CO2 intensity": "co2_emissions"})

    costs.at['OCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']
    costs.at['CCGT', 'co2_emissions'] = costs.at['gas', 'co2_emissions']

    costs.at['solar', 'capital_cost'] = 0.5*(costs.at['solar-rooftop', 'capital_cost'] +
                                             costs.at['solar-utility', 'capital_cost'])


    def costs_for_storage(store, link1, link2=None, max_hours=1.):
        capital_cost = link1['capital_cost'] + max_hours * store['capital_cost']
        if link2 is not None:
            capital_cost += link2['capital_cost']
        return pd.Series(dict(capital_cost=capital_cost,
                              marginal_cost=0.,
                              co2_emissions=0.))

    if elec_config is None:
        elec_config = snakemake.config['electricity']
    max_hours = elec_config['max_hours']
    costs.loc["battery"] = \
        costs_for_storage(costs.loc["battery storage"], costs.loc["battery inverter"],
                          max_hours=max_hours['battery'])
    costs.loc["H2"] = \
        costs_for_storage(costs.loc["hydrogen storage"], costs.loc["fuel cell"],
                          costs.loc["electrolysis"], max_hours=max_hours['H2'])

    for attr in ('marginal_cost', 'capital_cost'):
        overwrites = config.get(attr)
        if overwrites is not None:
            overwrites = pd.Series(overwrites)
            costs.loc[overwrites.index, attr] = overwrites

    return costs


def add_nice_carrier_names(n, config=None):
    if config is None: config = snakemake.config
    carrier_i = n.carriers.index
    nice_names = (pd.Series(config['plotting']['nice_names'])
                  .reindex(carrier_i).fillna(carrier_i.to_series().str.title()))
    n.carriers['nice_name'] = nice_names
    colors = pd.Series(config['plotting']['tech_colors']).reindex(carrier_i)
    if colors.isna().any():
        missing_i = list(colors.index[colors.isna()])
        logger.warning(f'tech_colors for carriers {missing_i} not defined '
                       'in config.')
    n.carriers['color'] = colors


def delete_stores(n, costs):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['Store']

    filtering = n.stores['carrier'] == 'H2'
    stores_i = n.stores[filtering].index

    if 'H2' in carriers:

        n.mremove("Bus", stores_i )

        n.mremove("Store", stores_i)

        n.mremove("Link", stores_i + " Electrolysis")

        n.mremove("Link", stores_i + " Fuel Cell")


def add_stores(n, costs, curves):
    elec_opts = snakemake.config['electricity']
    carriers = elec_opts['extendable_carriers']['Store']

    filtering2 = n.buses['carrier'] == 'AC'
    buses_i = n.buses[filtering2].index
    buses2 = n.buses[filtering2]
    bus_sub_dict = {k: buses2[k].values for k in ['x', 'y', 'country']}

    if 'H2' in carriers:
        h2_buses_i = n.madd("Bus", buses_i + " H2", carrier="H2", **bus_sub_dict)

        n.madd("Store", h2_buses_i,
               bus=h2_buses_i,
               carrier='H2',
               e_nom_extendable=True,
               e_cyclic=True,
               capital_cost=costs.at["hydrogen storage", "capital_cost"])
#               capital_cost=curves.at[int(snakemake.wildcards.iterations), 'parameter'])


        n.madd("Link", h2_buses_i + " Electrolysis",
               bus0=buses_i,
               bus1=h2_buses_i,
               carrier='H2 electrolysis',
               p_nom_extendable=True,
               efficiency=costs.at["electrolysis", "efficiency"],
               capital_cost=costs.at["electrolysis", "capital_cost"],
#               capital_cost=costs.at["electrolysis", "capital_cost"])
               marginal_cost=costs.at["electrolysis", "marginal_cost"])

        n.madd("Link", h2_buses_i + " Fuel Cell",
               bus0=h2_buses_i,
               bus1=buses_i,
               carrier='H2 fuel cell',
               p_nom_extendable=True,
               efficiency=costs.at["fuel cell", "efficiency"],
               #NB: fixed cost is per MWel
               capital_cost=costs.at["fuel cell", "capital_cost"],
#               capital_cost=costs.at["fuel cell", "capital_cost"] * costs.at["fuel cell", "efficiency"])
               marginal_cost=costs.at["fuel cell", "marginal_cost"])


def prepare_network(n, solve_opts):

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        n.madd("Generator", n.buses.index, " load",
               bus=n.buses.index,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
               )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components(n.one_port_components):
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                t.df['marginal_cost'] += (1e-2 + 2e-3 *
                                          (np.random.random(len(t.df)) - 0.5))

        for t in n.iterate_components(['Line', 'Link']):
            t.df['capital_cost'] += (1e-1 +
                2e-2*(np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760./nhours

    return n


def solve_network(n, config, solver_log=None, opts='', **kwargs):
    solver_options = config['solving']['solver'].copy()
    solver_name = solver_options.pop('name')
    cf_solving = config['solving']['options']
    track_iterations = cf_solving.get('track_iterations', False)
    min_iterations = cf_solving.get('min_iterations', 4)
    max_iterations = cf_solving.get('max_iterations', 6)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get('skip_iterations', False):
        network_lopf(n, solver_name=solver_name, solver_options=solver_options,
                     extra_functionality=None, **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=None, **kwargs)
    return n


if __name__ == "__main__":


    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        Path(tmpdir).mkdir(parents=True, exist_ok=True)

    opts = 'Co2L-24H'
    solve_opts = snakemake.config['solving']['options']
    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:
        n = pypsa.Network(snakemake.input[0])
        Nyears = n.snapshot_weightings.sum() / 8760.
        tech_costs = snakemake.input.tech_costs
        config = snakemake.config['costs']
        elec_config = snakemake.config['electricity']
        curves= pd.read_csv(snakemake.input.curves, index_col=list(range(1))).sort_index()
        costs = load_costs(Nyears, tech_costs=snakemake.input.tech_costs,
                           config=snakemake.config['costs'],
                           elec_config=snakemake.config['electricity'])
        delete_stores(n,costs)
        add_stores(n,costs,curves)
        add_nice_carrier_names(n, config=snakemake.config)
        n = prepare_network(n, solve_opts)
        n = solve_network(n, config=snakemake.config, solver_dir=tmpdir,
                        solver_log=snakemake.log.solver, opts=opts)
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
