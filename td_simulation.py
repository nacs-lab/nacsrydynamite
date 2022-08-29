from dynamite.operators import sigmax, sigmay, sigmaz
from dynamite.states import State
from dynamite import config

# import yaml which is not installed in Greg's directory

import sys
import os

path_to_pyyaml_dir = "/.local/lib/python3.8/site-packages/PyYAML-6.0-py3.8-linux-x86_64.egg/"

sys.path.append(os.path.expanduser("~") + path_to_pyyaml_dir)

import yaml

#################################################################
path_to_config_file = os.path.dirname(os.path.realpath(__file__)) + "/td_simulation_params.yml"

def get_fdict(fname):
    fhdl = open(fname, 'r')
    res = yaml.load(fhdl, Loader=yaml.FullLoader)
    fhdl.close()
    return res

params = get_fdict(path_to_config_file)

print(params["length"][0])
