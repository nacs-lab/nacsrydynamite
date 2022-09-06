from dynamite.operators import sigmax, sigmay, sigmaz, zero, identity
from dynamite.states import State
from dynamite.tools import mpi_print
from dynamite import config
config.initialize()
from petsc4py import PETSc

from libryd import iotools, science, utils, types
import os
import math

#################################################################
# CONFIGURATIONS AND HELPERS
# READ IN PARAMETERS
path_to_config_file = os.path.dirname(os.path.realpath(__file__)) + "/td_simulation_params.yml"

config_dict = iotools.get_fdict(path_to_config_file)
coords = iotools.get_arr_params(config_dict)
fname = iotools.get_config_params(config_dict)

start_t = 0
end_t = 0.0512 * 40 - 0.01
step_t = 0.0512
fi = -4
ff = 6
FT = 2
omega = 1.481
C6 = 26

ts = utils.gen_t_array(start_t, end_t, step_t)

det_fn = science.lin_ramp(fi, ff, FT)
omega_fn = science.constant(omega)
C6_fn = science.constant(C6)
dets = det_fn(ts)
omegas = omega_fn(ts)
C6s = C6_fn(ts)

config.L = len(coords)

#################################################################
# SCIENCE

init_state = State(state='D'*config.L)

measures = []
result = init_state.copy()
for i in range(len(ts)):
    #mpi_print("Calculating t: " + str(ts[i]))
    if i == 0:
        delta_t = ts[i]
        H = science.HRyd(coords, C6_fn(0), omega_fn(0), det_fn(0))
    else:
        delta_t = ts[i] - ts[i - 1]
        H = science.HRyd(coords, C6s[i-1], omegas[i-1], dets[i-1])
    H.evolve(init_state, t=delta_t, result=result)
    res, real_res = science.measure_n_op(result)
    measures.append(real_res)
    init_state, result = result, init_state

#################################################################
#OUTPUT
names = ["ts", "dets","omegas", "C6s", "sigmazs"]
datas = [ts, dets, omegas, C6s, measures]
if PETSc.COMM_WORLD.rank == 0:
    iotools.dump_to_file(fname, path_to_config_file, coords, names, datas)
