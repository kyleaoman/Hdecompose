"""Provide implementation of Springel & Hernquist 2003 neutral gas fraction algorithm."""

from astropy import units as U
from astropy.constants import m_p, k_B
import numpy as np


def sf_neutral_frac(
    fNeutral: np.ndarray,
    SFR: U.Quantity,
    u: U.Quantity,
    rho: U.Quantity,
    fH: float = 0.76,
    gamma: float = 5 / 3,
    Tc: U.Quantity = 1.0e3 * U.K,
    Th: U.Quantity = 5.73e7 * U.K,
    factorEVP: float = 573.0,
    rho_thresh: U.Quantity = 1.37e-1 * U.cm**-3,
) -> np.ndarray:
    """
    Compute particle neutral hydrogen fractions.

    Based on the multiphase ISM model of:
    Springel, V. and Hernquist, L. 2013, MNRAS, 339, 289.
    Precise calculation guided by notes provided by F. Marinacci.

    To compute neutral (HI + H_2) mass of particle, multiply NeutralFraction by
    Hydrogen mass fraction and particle mass.

    All arguments should be passed with units as applicable, use :mod:`~astropy.units`.

    Parameters
    ----------
    fNeutral : ~numpy.ndarray
        Gas neutral fractions from ionization model.
    SFR : ~astropy.units.Quantity
        Gas star formation rate.
    u : ~astropy.units.Quantity
        Gas specific internal energy.
    rho : ~astropy.units.Quantity
        Gas density.
    fH : float
        (Primordial) hydrogen abundance.
    gamma : float
        Adiabatic index, default 5/3.
    Tc : ~astropy.units.Quantity
        Temperature of cold clouds (Auriga: 1E3K).
    Th : ~astropy.units.Quantity
        Supernova temperature (Auriga: 5.73E7K).
    factorEVP : float
        Supernova evaporation parameter (Auriga: 573).
    rho_thresh : ~astropy.units.Quantity
        Density threshold for star formation (Auriga: 1.37E-1cm^-3).

    Returns
    -------
    ~numpy.ndarray
        An array of the same shape as particle property inputs containing the neutral mass
        fractions.
    """
    mu_neutral = 4 / (1 + 3 * fH)  # assumes fNeutral = 1.0
    uc = (k_B * Tc / mu_neutral / (gamma - 1) / m_p).to((U.km / U.s) ** 2)

    mu_ionized = 4 / (8 - 5 * (1 - fH))  # assumes fNeutral = 0.0
    uh = (k_B * Th / mu_ionized / (gamma - 1) / m_p).to((U.km / U.s) ** 2)

    uSN = (
        uh
        / (
            1
            + factorEVP
            * np.power((rho / m_p / rho_thresh).to(U.dimensionless_unscaled), -0.8)
        )
        + uc
    )

    retval = U.Quantity.copy(fNeutral)
    mask = SFR > 0
    retval[mask] = ((uSN - u) / (uSN - uc))[mask].to(U.dimensionless_unscaled)
    return retval
