from dynamite.tools import mpi_print

def fill_entry(this_list, idx, val):
    while (idx + 1) > len(this_list):
        this_list.append([])
    this_list[idx].append(val)

def vec_fill_entry(this_list, idx, vec):
    create = False
    while (idx + 1) > len(this_list):
        this_list.append([])
        create = True
    entry = this_list[idx]
    if create:
        for i in range(len(vec)):
            entry.append([vec[i]])
    else:
        if len(entry) != len(vec):
            raise Exception("Cannot add vector to list")
        else:
            for i in range(len(entry)):
                entry[i].append(vec[i])

def get_scan_params(params, vals):
    scanned_params = []
    dim_n = []
    param_vals = []
    for dimidx in range(len(params)):
        for paramidx in range(len(params[dimidx])):
            param_name = params[dimidx][paramidx]
            these_vals = vals[dimidx][paramidx]
            if param_name in scanned_params:
                raise Exception("Parameter scanned multiple times")
            else:
                if len(dim_n) > dimidx:
                    if len(these_vals) != dim_n[dimidx]:
                        raise Exception("Each scan dimension must have same number of scan parameters")
                else:
                    dim_n.append(len(these_vals))
                scanned_params.append(param_name)
                vec_fill_entry(param_vals, dimidx, these_vals)
    return scanned_params, param_vals

def make_oned_list(res, inlist):
    if len(inlist) == 1:
        # only one dimension
        new_entries = []
        if len(res) == 0:
            res.extend(inlist[0])
        else:
            for i in range(len(res)):
                for j in range(len(inlist[0])):
                    new_entries.append(res[i] + inlist[0][j])
            res.clear()
            res.extend(new_entries)
    else:
        make_oned_list(res, [inlist[0]])
        make_oned_list(res, inlist[1:])

def find_param_locs(desired_params, param_names):
    res = []
    for i in range(len(desired_params)):
        try:
            val = param_names.index(desired_params[i])
            res.append(val)
        except:
            res.append(-1)
    return res

def make_scan_dict():
    scan_dict = dict()
    scan_dict["scan"] = dict()
    return scan_dict

def gen_t_array(start_t, end_t, step_t):
    res = []
    t = start_t
    while t < end_t:
        res.append(t)
        t = t + step_t
    return res
