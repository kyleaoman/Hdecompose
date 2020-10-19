from astropy import units as U
from astropy.constants import m_p, k_B
import numpy as np


def sf_neutral_frac(
        fNeutral,
        SFR,
        u,
        rho,
        fH=.76,
        gamma=5 / 3,
        Tc=1.E3 * U.K,
        Th=5.73E7 * U.K,
        factorEVP=573,
        rho_thresh=1.37E-1 * U.cm ** -3
):

    """
    Computes particle neutral hydrogen fractions based on the multiphase ISM
    model of:
    Springel, V. and Hernquist, L. 2013, MNRAS, 339, 289.
    Precise calculation guided by notes provided by F. Marinacci.

    To compute neutral (HI + H_2) mass of particle, multiply NeutralFraction by
    Hydrogen mass fraction and particle mass.

    All arguments should be passed with units as applicable, use astropy.units.

    fNeutral:          Gas neutral fractions from ionization model.
    SFR:               Gas star formation rate.
    u:                 Gas specific internal energy.
    fH:                (Primordial) hydrogen abundance.
    gamma:             Adiabatic index, default 5/3.
    Tc:                Temperature of cold clouds (Auriga: 1E3K).
    Th:                Supernova temperature (Auriga: 5.73E7K).
    factorEVP:         Supernova evaporation parameter (Auriga: 573).
    rho_thresh:        Density threshold for star formation (Auriga:
                       1.37E-1cm^-3).

    Returns an array of the same shape as particle property inputs containing
    the neutral mass fractions.
    """

    mu_neutral = 4 / (1 + 3 * fH)  # assumes fNeutral = 1.0
    uc = (k_B * Tc / mu_neutral / (gamma - 1) / m_p).to((U.km / U.s) ** 2)

    mu_ionized = 4 / (8 - 5 * (1 - fH))  # assumes fNeutral = 0.0
    uh = (k_B * Th / mu_ionized / (gamma - 1) / m_p).to((U.km / U.s) ** 2)

    uSN = uh / (1 + factorEVP * np.power(
        (rho / m_p / rho_thresh).to(U.dimensionless_unscaled), -0.8)) + uc

    retval = U.Quantity.copy(fNeutral)
    mask = SFR > 0
    retval[mask] = ((uSN - u) / (uSN - uc))[mask].to(U.dimensionless_unscaled)
    return retval
