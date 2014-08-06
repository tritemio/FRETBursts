"""
This module contains functions to compute corrected FRET efficiency from
the proximity ratio and vice-versa.

For derivation see notebook: "FRET corrections algebra.ipynb"
"""

from __future__ import division
import numpy as np


def gamma_correct_E(Eraw, gamma):
    """Apply gamma correction to the uncorrected FRET `Eraw`."""
    Eraw = np.asarray(Eraw)
    return Eraw / (gamma - gamma*Eraw + Eraw)

def gamma_uncorrect_E(E, gamma):
    """Reverse gamma correction and return uncorrected FRET."""
    E = np.asarray(E)
    return gamma*E/(1 - E + gamma*E)

def correct_E_gamma_leak_dir(Eraw, gamma, leakage=0, dir_ex_t=0):
    """Compute corrected FRET efficency from proximity ratio `Eraw`.

    Paramaters:
        Eraw (float or array): proximity ratio (only background correction,
            no gamma, leakage or direct excitation)
        gamma (float): gamma factor
        leakage (float): leakage coefficient
        dir_ex_t (float): direct excitation coefficient expressed as
            n_dir = dir_ex_t * (na + gamma*nd). In terms of physical
            parameters it is the ratio of absorption cross-section at
            donor-excitation wavelengths of acceptor over donor.

    Returns
        Corrected FRET effciency
    """
    Eraw = np.asarray(Eraw)
    return (Eraw*(leakage + dir_ex_t*gamma + 1) - leakage - dir_ex_t*gamma) \
           / ( Eraw*(leakage -  gamma + 1) - leakage + gamma )

def uncorrect_E_gamma_leak_dir(E, gamma, leakage=0, dir_ex_t=0):
    """Compute proximity ratio from corrected FRET efficiency `E`.

    This function is the inverse of :func:`correct_E_gamma_leak_dir`.

    Paramaters:
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
    E = np.asarray(E)
    return (E*(gamma - leakage) + leakage + dir_ex_t*gamma) \
           / ( E*(gamma - leakage - 1) + leakage + dir_ex_t*gamma + 1 )


def test_fretmath():
    Ex = np.arange(-0.2, 1.2, 0.01)

    # Tests tautology for gamma = 1
    assert np.allclose(Ex, gamma_correct_E(Ex, gamma=1))
    assert np.allclose(Ex, correct_E_gamma_leak_dir(Ex, gamma=1))
    assert np.allclose(Ex, gamma_uncorrect_E(Ex, gamma=1))
    assert np.allclose(Ex, uncorrect_E_gamma_leak_dir(Ex, gamma=1))

    # Test gamma correction only
    leakage = 0.1
    dir_ex_t = 0.2
    for gamma in [0.2, 0.5, 0.8, 1.2, 1.8]:
        E1 = gamma_correct_E(Ex, gamma),
        E2 = correct_E_gamma_leak_dir(Ex, gamma)
        assert np.allclose(E1, E2)

        E1 = gamma_uncorrect_E(Ex, gamma),
        E2 = uncorrect_E_gamma_leak_dir(Ex, gamma)
        assert np.allclose(E1, E2)

        Ec = correct_E_gamma_leak_dir(Ex, gamma, leakage, dir_ex_t)
        Eu = uncorrect_E_gamma_leak_dir(Ec, gamma, leakage, dir_ex_t)
        assert np.allclose(Ex, Eu)


if __name__ == '__main__':
    test_fretmath()
    print 'All tests passed.'