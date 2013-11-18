"""
Simulate dark-counts timestamps.
"""

from numpy.random import random_integers

def gen_dcr(cps, time_s, clk_p):
    return random_integers(low=0, high=time_s/clk_p, size=int(cps*time_s))
