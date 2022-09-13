from dynamite.operators import sigmax, sigmay, sigmaz, zero, identity
from dynamite import config
import math

def n_op(i):
    return (identity() + sigmaz(i)) / 2

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

def measure_sigmaz(state):
    res = []
    real_res = []
    for i in range(config.L):
        val = state.dot(sigmaz(i) * state)
        res.append(state.dot(sigmaz(i)*state))
        real_res.append(val.real)
    return res, real_res

def measure_n_op(state):
    res = []
    real_res = []
    for i in range(config.L):
        val = state.dot(n_op(i) * state)
        res.append(state.dot(n_op(i)*state))
        real_res.append(val.real)
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
