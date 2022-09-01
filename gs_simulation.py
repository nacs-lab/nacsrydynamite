from dynamite.operators import sigmax, sigmay, sigmaz, zero, identity
from dynamite.states import State
from dynamite.tools import mpi_print
from dynamite import config
config.initialize()
from petsc4py import PETSc

from libryd import iotools, science, utils, types
import os

#################################################################
# CONFIGURATIONS AND HELPERS
# READ IN PARAMETERS
path_to_config_file = os.path.dirname(os.path.realpath(__file__)) + "/gs_simulation_params.yml"

config_dict = iotools.get_fdict(path_to_config_file)
coords = iotools.get_arr_params(config_dict)
fname = iotools.get_config_params(config_dict)
C6, omega, delta = iotools.get_ryd_params(config_dict)
# get scan params from file
#unique_scan_params, act_param_list = iotools.get_ryd_scan_params(config_dict)

# make scan params on your own
scan_dict = utils.make_scan_dict()
scan_dict["scan"]["C6"] = dict()
scan_dict["scan"]["C6"]["dim"] = 1
scan_dict["scan"]["C6"]["vals"] = [-30, 0, 30]
unique_scan_params, act_param_list = iotools.get_ryd_scan_params(scan_dict)

config.L = len(coords)

#################################################################
# SCIENCE

scan_param_idxs = utils.find_param_locs([types.ParamType.C6, types.ParamType.Omega, types.ParamType.Delta], unique_scan_params)

gs_energies = []
gs_measures = []

for i in range(len(act_param_list)):
    mpi_print("Working on " + str(act_param_list[i]))
    real_res = []
    gs_energy = -1
    if scan_param_idxs[0] == -1:
        C6ToUse = C6
    else:
        C6ToUse = C6 + act_param_list[i][scan_param_idxs[0]]
    if scan_param_idxs[1] == -1:
        omegaToUse = omega
    else:
        omegaToUse = omega + act_param_list[i][scan_param_idxs[1]]
    if scan_param_idxs[2] == -1:
        deltaToUse = delta
    else:
        deltaToUse = delta + act_param_list[i][scan_param_idxs[2]]

    H = science.HRyd(coords, C6ToUse, omegaToUse, deltaToUse)
    eigvals, eigvecs = H.eigsolve(getvecs=True, nev=4)
    if len(eigvals) > 0:
        gs = eigvecs[0]
        gs_energy = eigvals[0]
        res, real_res = science.measure_sigmaz(gs)
    gs_energies.append(gs_energy)
    gs_measures.append(real_res)

names = ["scan_params", "param_list","energies", "measurez"]
datas = [unique_scan_params, act_param_list, gs_energies, gs_measures]
if PETSc.COMM_WORLD.rank == 0:
    iotools.dump_to_file(fname, path_to_config_file, coords, names, datas)
