# import yaml which is not installed in Greg's directory

from libryd import types, utils
from dynamite.tools import mpi_print
import numpy as np
import sys
import os

path_to_pyyaml_dir = "/.local/lib/python3.8/site-packages/PyYAML-6.0-py3.8-linux-x86_64.egg/"

sys.path.append(os.path.expanduser("~") + path_to_pyyaml_dir)

import yaml

def get_fdict(fname):
    fhdl = open(fname, 'r')
    res = yaml.load(fhdl, Loader=yaml.FullLoader)
    fhdl.close()
    return res

def get_arr_params(param_dict):
    # array parameters
    arr_info = param_dict["arr"]
    arr_type = arr_info["type"]
    if arr_type == "rect":
        arr_size = arr_info["size"]
        arr_spacing = arr_info["spacing"]
        if len(arr_size) != 2:
            raise Exception("Please input two values for rect array size")
        if len(arr_spacing) != 2:
            raise Exception("Please input two values for rect array spacing")
        # build coordinates for rect array
        coords = []
        for i in range(arr_size[0]):
            for j in range(arr_size[1]):
                coords.append([arr_spacing[0] * i, arr_spacing[1] * j])
    elif arr_type == "circ":
        arr_rad = arr_info["radius"]
        arr_tot_sites = arr_info["n_sites"]
        arr_phase = arr_info["phase"]
        arr_fill_sites = arr_info["filled_sites"]
        arr_center = arr_info["center"]
        if len(arr_center) != 2:
            raise Exception("Please input two values for circle center")
        coords = []
        dtheta = 2 * np.pi / (arr_tot_sites - 1)
        theta = arr_phase
        for i in range(arr_tot_sites):
            if i in arr_fill_sites:
                coords.append([arr_rad * np.cos(theta) + arr_center[0], arr_rad * np.sin(theta) + arr_center[1]])
            theta = theta + dtheta
    else:
        coords = arr_info["coords"]
    return coords

def get_config_params(param_dict):
    config_params = param_dict["config"]
    fname = config_params["out_file"]
    return fname

# read in params here
def get_ryd_params(param_dict):
    parameters = param_dict["params"]
    # get other parameters
    omega = parameters["omega"]
    delta = parameters["delta"]
    C6 = parameters["C6"]
    return C6, omega, delta

def get_ryd_scan_params_raw(param_dict):
    params = [] # nested list with param types
    vals = [] # nested list of values for each param
    if "scan" in param_dict.keys():
        scan_dict = param_dict["scan"]
        for key in scan_dict.keys():
            info = scan_dict[key]
            if "dim" not in info.keys():
                mpi_print("WARNING: dimension not specified for scan param " + str(idx) + ". Ignoring this param.")
                continue
            if "vals" not in info.keys():
                mpi_print("WARNING: parameter vals not specified for scan param " + str(idx) + ". Ignoring this param.")
                continue
            dim = info["dim"]
            param = key
            these_vals = info["vals"]
            utils.fill_entry(params, dim - 1, types.ParamType.fromStr(param))
            utils.fill_entry(vals, dim - 1, these_vals)
    return params, vals

def get_ryd_scan_params(config_dict):
    scan_param_names,scan_vals = get_ryd_scan_params_raw(config_dict)
    unique_scan_params, param_list = utils.get_scan_params(scan_param_names, scan_vals)
    act_param_list = []
    if len(unique_scan_params) > 0:
        utils.make_oned_list(act_param_list, param_list)
    return unique_scan_params, act_param_list


def dump_to_file(fname, config_fname, coords, names, datas):
    # open both files
    config_fhdl = open(config_fname, "r")
    fhdl = open(fname, "w")
    # read content from first file
    for line in config_fhdl:
        # append content to second file
        fhdl.write(line)
    #yaml.dump(input_dict,fhdl, default_flow_style=None)
    fhdl.write("\n" + "out_data: " + "\n" + "  coords: " + str(coords) + "\n")
    for i in range(len(names)):
        fhdl.write("  " + names[i] + ": " + str(datas[i]) + "\n")
    fhdl.close()
    config_fhdl.close()
