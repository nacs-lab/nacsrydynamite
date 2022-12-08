from dynamite.operators import sigmax, sigmay, sigmaz, zero, identity
from dynamite import config
import math

def n_op(i):
    return (identity() + sigmaz(i)) / 2

def sigma_bond_op(i):
    return (-1)**i * (n_op(i) - n_op(i + 1))

def HRyd(coords, C6, omega, delta):
    # delta, omega and C6 should all be in units of real frequency
    # coords and C6 need to agree on their units
    n_atoms = len(coords)
    if type(omega) != list:
        omega = [omega for i in range(n_atoms)]
    if type(delta) != list:
        delta = [delta for i in range(n_atoms)]
    res = zero()
    for i in range(n_atoms):
        res = res + 2 * math.pi * omega[i] / 2 * sigmax(i) - 2 * math.pi * delta[i] * n_op(i)
        for j in range(i+1, n_atoms):
            res = res + 2 * math.pi * C6 / (math.dist(coords[i], coords[j]))**6 * n_op(i) * n_op(j)
    return res

def HRydWithPhase(coords, C6, omega, phase, delta):
    # delta, omega and C6 should all be in units of real frequency
    # coords and C6 need to agree on their units
    n_atoms = len(coords)
    if type(omega) != list:
        omega = [omega for i in range(n_atoms)]
    if type(delta) != list:
        delta = [delta for i in range(n_atoms)]
    if type(phase) != list:
        phase = [phase for i in range(n_atoms)]
    res = zero()
    for i in range(n_atoms):
        res = res + 2 * math.pi * omega[i] * math.cos(phase[i]) / 2 * sigmax(i) + 2 * math.pi * omega[i] * math.sin(phase[i]) / 2 * sigmay(i) - 2 * math.pi * delta[i] * n_op(i)
        for j in range(i+1, n_atoms):
            res = res + 2 * math.pi * C6 / (math.dist(coords[i], coords[j]))**6 * n_op(i) * n_op(j)
    return res

def measure_sigmaz(state):
    res = []
    real_res = []
    for i in range(config.L):
        val = state.dot(sigmaz(i) * state)
        res.append(val)
        real_res.append(val.real)
    return res, real_res

def measure_n_op(state):
    res = []
    real_res = []
    for i in range(config.L):
        val = state.dot(n_op(i) * state)
        res.append(val)
        real_res.append(val.real)
    return res, real_res

def measure_sigma_field_open(state):
    res = []
    real_res = []
    for i in range(config.L - 1):
        op = (-1)**i * (n_op(i) - n_op(i + 1))
        val = state.dot(op * state)
        res.append(val)
        real_res.append(val.real)
    return res, real_res

def measure_sigma_field_corr_open(state):
    res = []
    real_res = []
    for i in range(config.L - 1):
        entry = []
        real_entry = []
        for j in range(i, config.L - 1):
            op = sigma_bond_op(i) * sigma_bond_op(j)
            val = state.dot(op * state)
            entry.append(val)
            real_entry.append(val.real)
        res.append(entry)
        real_res.append(real_entry)
    return res, real_res

def measure_sum_bond(state):
    op = zero()
    for i in range(config.L - 1):
        op = op + (-1)**i * (n_op(i) - n_op(i + 1))
    res = state.dot(op*state)
    real_res = res.real
    return res, real_res

def vectorize(func):
    def f(ts):
        if type(ts) == list:
            res = []
            for t in ts:
                res.append(func(t))
            return res
        else:
            return func(ts)
    return f

def lin_ramp(start, stop, T):
    # ramp from start to stop in time T
    @vectorize
    def fn(t):
        if t < 0:
            return start
        elif t < T:
            return start + t / T * (stop - start)
        else:
            return stop
    return fn

def constant(val):
    @vectorize
    def fn(t):
        return val
    return fn

def exp_in_exp_out(start, stop1, t1, tau1, stop2, t2, tau2):
    @vectorize
    def fn(t):
        if t < 0:
            return start
        elif t < t1:
            exp1 = math.exp(t1 / tau1)
            exp1m = exp1 - 1
            return -exp1 * (stop1 - start) / exp1m * math.exp(-t / tau1) + (exp1 * stop1 - start) / exp1m
        elif t < t1 + t2:
            exp2 = math.exp(t2 /tau2)
            exp2m = exp2 - 1
            exp2a = math.exp((t1 + t2) / tau2)
            exp2b = math.exp(t1 / tau2)
            return (stop2 - stop1) / exp2m * math.exp((t - t1) / tau2) + (-exp2a * stop1 + exp2b * stop2) / (exp2b-exp2a)
        else:
            return stop2
    return fn

def lila_in(E0, Ec, delta_i, delta_c, tf):
    @vectorize
    def fn(t):
        if t < 0:
            return delta_i
        elif t < tf:
            return (E0 * delta_c * t + Ec * delta_i * (tf - t)) / (E0 * t + Ec * (tf - t))
        else:
            return delta_c
    return fn

def lila_out(Ef, Ec, delta_f, delta_c, tf):
    @vectorize
    def fn(t):
        if t < 0:
            return delta_c
        elif t < tf:
            return (Ef * delta_c * (tf - t) + Ec * delta_f * t) / (Ef * (tf - t) + Ec * t)
        else:
            return delta_f
    return fn

def lila_in_and_out(Ei, Ec, Ef, delta_i, delta_c, delta_f, t1, t2):
    @vectorize
    def fn(t):
        if t < 0:
            return delta_i
        elif t < t1:
            fn_hdl = lila_in(Ei, Ec, delta_i, delta_c, t1)
            return fn_hdl(t)
        elif t < t1 + t2:
            fn_hdl = lila_out(Ef, Ec, delta_f, delta_c, t2)
            return fn_hdl(t - t1)
        else:
            return delta_f
    return fn

def sine(amp, freq, phase):
    @vectorize
    def fn(t):
        return amp * math.sin(2 * math.pi * freq * t + phase)
    return fn

def piecewise_lin(fs, ts):
    #fs and ts are both lists. first ramp will be from fs[0] to fs[1] in time ts[0] to ts[1], etc.
    # assume ts are sorted.
    @vectorize
    def fn(t):
        idx = 0
        while True:
            if idx == 0:
                if t < ts[idx]:
                    return fs[idx]
                else:
                    idx = idx + 1
            elif idx >= len(ts):
                return fs[-1]
            else:
                if t < ts[idx]:
                    # t is already  greater than ts[idx - 1]
                    return fs[idx - 1] + (fs[idx] - fs[idx - 1]) / (ts[idx] - ts[idx - 1]) * (t - ts[idx - 1])
                else:
                    idx = idx + 1
    return fn

def piecewise_jump(vals, ts):
    # should have one more vals element than ts elements
    @vectorize
    def fn(t):
        idx = 0
        while True:
            if idx == 0:
                if t < ts[idx]:
                    return vals[idx]
                else:
                    idx = idx + 1
            elif idx >= len(ts):
                return vals[-1]
            else:
                if t < ts[idx]:
                    # t is already  greater than ts[idx - 1]
                    return vals[idx]
                else:
                    idx = idx + 1
    return fn

def piecewise_fn(fns, ts):
    # should have one more fns element than ts elements
    # fn use the time from the beginning of the particular segment
    @vectorize
    def fn(t):
        idx = 0
        while True:
            if idx == 0:
                if t < ts[idx]:
                    return fns[0](t)
                else:
                    idx = idx + 1
            elif idx >= len(ts):
                return fns[-1](t - ts[-1])
            else:
                if t < ts[idx]:
                    # t is already  greater than ts[idx - 1]
                    return fns[idx](t - ts[idx - 1])
                else:
                    idx = idx + 1
    return fn

