# encoding: utf-8
#
# FRETBursts - A single-molecule FRET burst analysis toolkit.
#
# Copyright (C) 2014-2016 The Regents of the University of California,
#               Antonino Ingargiola <tritemio@gmail.com>
#
"""
The `fretmath` module contains functions to compute corrected FRET efficiency
from the proximity ratio and vice-versa.

For derivation see notebook: "Derivation of FRET and S correction formulas.ipynb"
(`link <http://nbviewer.jupyter.org/github/tritemio/notebooks/blob/master/Derivation%20of%20FRET%20and%20S%20correction%20formulas.ipynb>`__).

"""

from __future__ import division, print_function
from builtins import range, zip

import numpy as np

##
# Comple corrections
#
def correct_E_gamma_leak_dir(Eraw, gamma=1, leakage=0, dir_ex_t=0):
    """Compute corrected FRET efficency from proximity ratio `Eraw`.

    For the inverse function see :func:`uncorrect_E_gamma_leak_dir`.

    Arguments:
        Eraw (float or array): proximity ratio (only background correction,
            no gamma, leakage or direct excitation)
        gamma (float): gamma factor
        leakage (float): leakage coefficient
        dir_ex_t (float): coefficient expressing the direct excitation as
            n_dir = dir_ex_t * (na + gamma*nd). In terms of physical
            parameters it is the ratio of acceptor over donor absorption
            cross-sections at the donor-excitation wavelength.

    Returns
        Corrected FRET effciency
    """
    if isinstance(Eraw, (list, tuple)):
        Eraw = np.asarray(Eraw)
    return (Eraw*(leakage + dir_ex_t*gamma + 1) - leakage - dir_ex_t*gamma) \
           / (Eraw*(leakage -  gamma + 1) - leakage + gamma)

def uncorrect_E_gamma_leak_dir(E, gamma=1, leakage=0, dir_ex_t=0):
    """Compute proximity ratio from corrected FRET efficiency `E`.

    This function is the inverse of :func:`correct_E_gamma_leak_dir`.

    Arguments:
        E (float or array): corrected FRET efficiency
        gamma (float): gamma factor
        leakage (float): leakage coefficient
        dir_ex_t (float): direct excitation coefficient expressed as
            n_dir = dir_ex_t * (na + gamma*nd). In terms of physical
            parameters it is the ratio of absorption cross-section at
            donor-excitation wavelengths of acceptor over donor.

    Returns
        Proximity ratio (reverses gamma, leakage and direct excitation)
    """
    if isinstance(E, (list, tuple)):
        E = np.asarray(E)
    return (E*(gamma - leakage) + leakage + dir_ex_t*gamma) \
           / (E*(gamma - leakage - 1) + leakage + dir_ex_t*gamma + 1)


##
# Single parameter correction
#
# NOTE: Applying the single parameter correction one after the other
#       IS NOT equivalent to applying the comple correction in one step.
#       In fact chaining the single correction gives mathematically the wrong
#       result, even if it may be in some cases close enough to the
#       mathematically accurate value.
#
def gamma_correct_E(Eraw, gamma):
    """Apply gamma correction to the uncorrected FRET `Eraw`.

    For the inverse see :func:`gamma_uncorrect_E`.
    """
    if isinstance(Eraw, (list, tuple)):
        Eraw = np.asarray(Eraw)
    return Eraw / (gamma - gamma*Eraw + Eraw)

def gamma_uncorrect_E(E, gamma):
    """Reverse gamma correction and return uncorrected FRET.

    For the inverse see :func:`gamma_correct_E`.
    """
    if isinstance(E, (list, tuple)):
        E = np.asarray(E)
    return gamma*E/(1 - E + gamma*E)


def leakage_correct_E(Eraw, leakage):
    """Apply leakage correction to the uncorrected FRET `Eraw`.

    For the inverse see :func:`leakage_uncorrect_E`.
    """
    if isinstance(Eraw, (list, tuple)):
        Eraw = np.asarray(Eraw)
    return (Eraw*(leakage + 1) - leakage) / (Eraw*leakage - leakage + 1)

def leakage_uncorrect_E(E, leakage):
    """Reverse leakage correction and return uncorrected FRET.

    For the inverse see :func:`leakage_correct_E`.
    """
    if isinstance(E, (list, tuple)):
        E = np.asarray(E)
    return (E + leakage - E*leakage) / (1 + leakage - E*leakage)


def dir_ex_correct_E(Eraw, dir_ex_t):
    """Apply direct excitation correction to the uncorrected FRET `Eraw`.

    The coefficient `dir_ex_t` expresses the direct excitation as
    n_dir = dir_ex_t * (na + gamma*nd). In terms of physical
    parameters it is the ratio of acceptor over donor absorption
    cross-sections at the donor-excitation wavelength.

    For the inverse see :func:`dir_ex_uncorrect_E`.
    """
    if isinstance(Eraw, (list, tuple)):
        Eraw = np.asarray(Eraw)
    return Eraw*(dir_ex_t + 1) - dir_ex_t

def dir_ex_uncorrect_E(E, dir_ex_t):
    """Reverse direct excitation correction and return uncorrected FRET.

    For the inverse see :func:`dir_ex_correct_E`.
    """
    if isinstance(E, (list, tuple)):
        E = np.asarray(E)
    return (E + dir_ex_t) / (dir_ex_t + 1)

def correct_S(Eraw, Sraw, gamma, leakage, dir_ex_t):
    """Correct S values for gamma, leakage and direct excitation.

    Arguments:
        Eraw (scalar or array): uncorrected ("raw") E after only background
            correction (no gamma, leakage or direct excitation).
        Sraw (scalar or array): uncorrected ("raw") S after only background
            correction (no gamma, leakage or direct excitation).
        gamma (float): gamma factor.
        leakage (float): donor emission leakage into the acceptor channel.
        dir_ex_t (float): direct acceptor excitation by donor laser.
            Defined as ``n_dir = dir_ex_t * (na + g nd)``. The dir_ex_t
            coefficient is the ratio between D and A absorbtion cross-sections
            at the donor-excitation wavelength.

    Returns
        Corrected S (stoichiometry), same size as `Sraw`.
    """
    if isinstance(Eraw, (list, tuple)):
        Eraw = np.asarray(Eraw)
    if isinstance(Sraw, (list, tuple)):
        Sraw = np.asarray(Sraw)
    return (Sraw*(Eraw*leakage - Eraw*gamma + Eraw - leakage + gamma) /
            (Eraw*leakage*Sraw - Eraw*Sraw*gamma + Eraw*Sraw - leakage*Sraw -
             Sraw*dir_ex_t + Sraw*gamma - Sraw + dir_ex_t + 1))

def uncorrect_S(E_R, S, gamma, L_k, d_dirT):
    """Function used to test :func:`correct_S`."""
    return (S*(d_dirT + 1) /
            (-E_R*L_k*S + E_R*L_k + E_R*S*gamma - E_R*S - E_R*gamma +
             E_R + L_k*S - L_k + S*d_dirT - S*gamma + S + gamma))

def test_fretmath():
    """Run a few consistency checks for the correction functions.
    """
    Ex = np.arange(-0.2, 1.2, 0.01)

    # Tests tautology for gamma = 1
    assert np.allclose(Ex, correct_E_gamma_leak_dir(Ex, gamma=1))
    assert np.allclose(Ex, uncorrect_E_gamma_leak_dir(Ex, gamma=1))
    assert np.allclose(Ex, gamma_correct_E(Ex, gamma=1))
    assert np.allclose(Ex, gamma_uncorrect_E(Ex, gamma=1))

    # Tests tautology for leakage = 0
    assert np.allclose(Ex, leakage_correct_E(Ex, leakage=0))
    assert np.allclose(Ex, leakage_uncorrect_E(Ex, leakage=0))

    # Tests tautology for dir_ex_t = 0
    assert np.allclose(Ex, dir_ex_correct_E(Ex, dir_ex_t=0))
    assert np.allclose(Ex, dir_ex_uncorrect_E(Ex, dir_ex_t=0))

    # Test round-trip consistency
    for leakage in [0.01, 0.04, 0.2]:
        for dir_ex_t in [0.01, 0.03, 0.09]:
            for gamma in [0.2, 0.5, 0.8, 1.2, 1.8]:
                Ec = correct_E_gamma_leak_dir(Ex, gamma, leakage, dir_ex_t)
                Eu = uncorrect_E_gamma_leak_dir(Ec, gamma, leakage, dir_ex_t)
                assert np.allclose(Ex, Eu)

    for gamma in [0.2, 0.5, 0.8, 1.2, 1.8]:
        # Test gamma correction
        E1 = gamma_correct_E(Ex, gamma)
        E2 = correct_E_gamma_leak_dir(Ex, gamma)
        assert np.allclose(E1, E2)
        E1 = gamma_uncorrect_E(Ex, gamma)
        E2 = uncorrect_E_gamma_leak_dir(Ex, gamma)
        assert np.allclose(E1, E2)

    for leakage in [0.01, 0.04, 0.2]:
        # Test leakage correction
        E1 = leakage_correct_E(Ex, leakage)
        E2 = correct_E_gamma_leak_dir(Ex, leakage=leakage)
        assert np.allclose(E1, E2)
        E1 = leakage_uncorrect_E(Ex, leakage)
        E2 = uncorrect_E_gamma_leak_dir(Ex, leakage=leakage)
        assert np.allclose(E1, E2)

    for dir_ex_t in [0.01, 0.03, 0.09]:
        # Test direct excitation correction
        E1 = dir_ex_correct_E(Ex, dir_ex_t)
        E2 = correct_E_gamma_leak_dir(Ex, dir_ex_t=dir_ex_t)
        assert np.allclose(E1, E2)
        E1 = dir_ex_uncorrect_E(Ex, dir_ex_t)
        E2 = uncorrect_E_gamma_leak_dir(Ex, dir_ex_t=dir_ex_t)
        assert np.allclose(E1, E2)

    # Test S correction
    np.random.seed(1)
    Ex = np.arange(-0.2, 1.2, 0.01)
    Sx = np.arange(-0.2, 1.2, 0.01)
    np.random.shuffle(Ex)
    np.random.shuffle(Sx)

    gamma_ = 1
    leakage_ = 0
    dir_ex_t_ = 0
    S_corr = correct_S(Ex, Sx, gamma_, leakage_, dir_ex_t_)
    S_uncorr = uncorrect_S(Ex, S_corr, gamma_, leakage_, dir_ex_t_)
    assert np.allclose(S_corr, Sx)
    assert np.allclose(S_uncorr, Sx)

    gamma_ = 0.7
    leakage_ = 0.05
    dir_ex_t_ = 0.1
    S_corr = correct_S(Ex, Sx, gamma_, leakage_, dir_ex_t_)
    S_uncorr = uncorrect_S(Ex, S_corr, gamma_, leakage_, dir_ex_t_)
    assert np.allclose(S_uncorr, Sx)
    assert (S_corr.min() > -0.5) and (S_corr.max() < 1.5)

if __name__ == '__main__':
    test_fretmath()
    print('All tests passed.')
