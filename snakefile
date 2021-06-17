# SPDX-FileCopyrightText: : 2017-2020 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: GPL-3.0-or-later

from os.path import normpath, exists
from shutil import copyfile
import pandas as pd
import numpy as np

configfile: "config.yaml"

COSTS="data/costs.csv"
CURVE="data/curve.csv"


rule solve_all_iterations:
    input: expand("iterations/elec_s_5_ec_lcopt_Co2L-24H_{iterations}.nc", iterations=range(3))

# Write the CORRECT corresponding name of the resulting network.

def input_solve_iterations(w):
    # It's mildly hacky to include the separate costs input as first entry
    n = int(w.iterations)
    if n == 0:
        return "pypsa-eur/networks/elec_s_5_ec_lcopt_Co2L-24H.nc"
    else:
        Z = n - 1
        return ("pypsa-eur/networks/elec_s_5_ec_lcopt_Co2L-24H.nc", "iterations/elec_s_5_ec_lcopt_Co2L-24H_{}.nc".format(Z))


rule solve_iterations:
    input:
        input_solve_iterations,
        curves=CURVE,
        tech_costs=COSTS
    output: "iterations/elec_s_5_ec_lcopt_Co2L-24H_{iterations}.nc"
    log:
        solver=normpath("logs/solve_iterations/elec_s_5_ec_lcopt_Co2L-24H_{iterations}_solver.log"),
        python="logs/solve_iterations/elec_s_5_ec_lcopt_Co2L-24H_{iterations}_python.log",
        memory="logs/solve_iterations/elec_s_5_ec_lcopt_Co2L-24H_{iterations}_memory.log"
    threads: 4
    script: "scripts/solve_iterations.py"
    #notebook: "scripts/solve_iterations.py.ipynb"
