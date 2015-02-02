#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014 Antonino Ingargiola <tritemio@gmail.com>
#
"""
From this module import the find_optimal_T() function to calculate the
T value for burst search that gives a probability to detect noise as burst
below a threshold (P_user), for a given value of m (burst search parameter).

Implementation details:
----------------------
Fixed a rate for a Poisson process, I compare the probability to exceed a
minimum rate using different window sizes (T).

Poisson process rate: lam or rate or bg_rate
Time window duration: T
Poisson distribution parameter: n = lam*T
"""
from __future__ import print_function
from builtins import range, zip

import numpy as np
from scipy.stats import poisson, chi2, erlang

def find_optimal_T_bga(bg_array, m, P):
    """Return T so that m-ph delay from pure BG will be < T with prob. < P.

    This calls find_optimal_T() for each value in `bg_array`.
    """
    TT = np.array([find_optimal_T(bg_i, m, P) for bg_i in bg_array])
    return TT

def find_optimal_T(bg_rate, m, P):
    """Return T so that m-ph delay from pure BG will be < T with prob. < P.

    This is equivalent to find_optimal_T_chi2 but conceptually simpler.
    """
    T = erlang.ppf(P, m, scale=1./bg_rate)
    return T

def find_optimal_T_chi2(bg_rate, m, P):
    """Returns the min. T so that bg_rate is <= than lower_conf_rate(m,T,P).

    This is equivalent but much faster than find_optimal_T_iter().
    Note: This is based on the confidence intervall of multiple exponential
    """
    T = 0.5*chi2(2*m).ppf(P)/bg_rate
    return T

def find_optimal_threshold(m, P):
    """Returns the min. threshold to have prob. < P to be BG (averaging m ph).

    Same formula as find_optimal_T() (must be multiplied by bg to have the rate.
    """
    return m/(0.5*chi2(2*m).ppf(P))

def prob_noise_above_th(rate, T, m):
    """Returns the probability that noise is above the burst search threshold.

    Basically is the probability that a poisson process with rate "rate" had
    "m" or more events in a time window "T".
    """
    return poisson(rate*T).sf(m-1)
_p = prob_noise_above_th

def prob_noise_above_th_test_version(rate, T, m):
    """Returns the probability that noise is above the burst search threshold.

    Same as prob_noise_above_th() but compute ARITHMETICALLY the probability
    adding up the PDF from m to a very high integer (mean + 9std_dev).

    WARNING: Test version, significantly slower!!
    """
    rate = float(rate)
    T = float(T)
    n_poiss = rate*T
    if m < (n_poiss) - 9*np.sqrt(n_poiss): return 1. #

    s = 0.
    for i in np.arange(m,int(rate*T + 9*np.sqrt(rate*T))+1):
        s += reduce(lambda x,y: x*y, [(n_poiss/j) for j in range(1, i+1)])
        assert np.isfinite(s)
    p = np.exp(-(rate*T)) * s
    assert np.isfinite(p)
    return p
_p2 = prob_noise_above_th_test_version

def find_optimal_T_iter(bg_rate, m, P_user, max_iter=int(1e6), debug=False):
    """Returns the T (sec.) that assure prob. P_user or less to have noise.

    Given a background rate (bg_rate) find a T such as Poiss{bg_rate*T} has
    probability to have m or more events (k >= m) lesser than P_user.

    In other words find the max T that satisfies: poisson(lam*T).sf(m) < P_user
    """

    # Initial T in seconds
    T = 1.  # start from a huge T window that gives an enormous rate such as
            # lam*T is so high that poisson(lam*T).sf(m) would be for sure >
            # P_user, then decrease T until poisson(lam*T).sf(m) < T

    converged = False
    for i in range(max_iter):
        prob = prob_noise_above_th(bg_rate,T,m)
        if prob > P_user:
            T /= 2.
        else:
            converged = True
            break
    if not converged: raise StopIteration
    if debug: print("T_min = %.3f ms, T_max = %.3f ms" % (T*1e3, 2*T*1e3))

    step = T/1000.
    T = 2*T
    converged = False
    for i in range(max_iter):
        prob = prob_noise_above_th(bg_rate,T,m)
        if prob > P_user:
            T -= step
        else:
            converged = True
            break
    if not converged: raise StopIteration
    if debug: print(" >T_min = %.3f ms, T_max = %.3f ms" % (T*1e3, (T+step)*1e3))
    return T

def test_find_optimal_T_iter(P_user, debug=False):
    Lam = np.array([100, 1000, 2000, 5000])
    M = np.array([3,5,10,20,50])

    TT = np.zeros((Lam.size,M.size))
    for y,lam in enumerate(Lam):
        for x,m in enumerate(M):
            TT[y,x] = find_optimal_T_iter(lam, m, P_user, debug=debug)
            if debug: print("  '-> rate = %d, m = %d" % (lam, m))
    print("\n P_user = %.2f %%, T values in ms:" % (P_user*1e2))
    print(((" "*8) + " %7d " * M.size) % tuple(M))
    for y,lam in enumerate(Lam):
        print(" %7d " % lam,)
        print((" %7.3f "*M.size) % tuple(TT[y,:]*1e3))
    return TT

def test_find_optimal_T(P_user, debug=False):
    Lam = np.array([100, 1000, 2000, 5000])
    M = np.array([3,5,10,20,50])

    TT = np.zeros((Lam.size,M.size))
    for y,lam in enumerate(Lam):
        for x,m in enumerate(M):
            TT[y,x] = find_optimal_T(lam, m, P_user)
            if debug: print("  '-> rate = %d, m = %d" % (lam, m))
    print("\n P_user = %.2f %%, T values in ms:" % (P_user*1e2))
    print(((" "*8) + " %7d " * M.size) % tuple(M))
    for y,lam in enumerate(Lam):
        print(" %7d " % lam,)
        print((" %7.3f "*M.size) % tuple(TT[y,:]*1e3))
    return TT


def old_hard_threshold_T(factor=2.5):
    Lam = np.array([100, 1000, 2000, 5000])
    M = np.array([3,5,10,20,50])

    TT = np.zeros((Lam.size,M.size))
    for y,lam in enumerate(Lam):
        for x,m in enumerate(M):
            TT[y,x] = float(m)/(factor*lam)
    print("\n T values in ms, OLD HARD-THRESHOLD METHOD:")
    print(((" "*8) + " %7d " * M.size) % tuple(M))
    for y,lam in enumerate(Lam):
        print(" %7d " % lam,)
        print((" %7.3f "*M.size) % tuple(TT[y,:]*1e3))
    return TT


if __name__ == '__main__':
    TT_10 = test_find_optimal_T(0.1)    # 10%
    TT_01 = test_find_optimal_T(0.01)   # 1%
    TT_001 = test_find_optimal_T(0.001)  # 0.1%

    TT_old = old_hard_threshold_T()

    plot(TT_10.ravel(), lw=1.5, label='P = 10%')
    plot(TT_01.ravel(), lw=1.5, label='P = 1%')
    plot(TT_001.ravel(), lw=1.5, label='P = 0.1%')
    plot(TT_old.ravel(), lw=2, color='k', label='T = m/(2.5*bg_rate)')
    legend(loc='best')
    grid(True)
    xlabel("m = "+str([3,5,10,20,50])+"; rate = "+str([100, 1000, 2000, 5000]))
    ylabel('T (s)')
    title("T selection comparison")

    # Reference results for P = 10%, 1%, 0.1%
    Output = """
 P_user = 10.00 %, T values in ms:
               3        5       10       20       50
     100    11.016   24.313   62.188  145.250  411.750
    1000     1.102    2.432    6.219   14.523   41.156
    2000     0.551    1.216    3.109    7.262   20.578
    5000     0.220    0.486    1.244    2.904    8.234

 P_user = 1.00 %
               3        5       10       20       50
     100     4.359   12.789   41.281  110.813  350.250
    1000     0.436    1.278    4.129   11.078   35.031
    2000     0.218    0.639    2.064    5.539   17.516
    5000     0.087    0.256    0.826    2.215    7.004

 P_user = 0.10 %
               3        5       10       20       50
     100     1.905    7.391   29.594   89.563  309.500
    1000     0.190    0.739    2.959    8.953   30.953
    2000     0.095    0.370    1.479    4.477   15.477
    5000     0.038    0.148    0.592    1.791    6.191

"""
