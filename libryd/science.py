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
