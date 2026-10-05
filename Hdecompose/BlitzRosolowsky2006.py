"""Provide an implementation of the Blitz & Rosolowsky molecular fraction algorithm."""

import numpy as np
from astropy import units as U
from astropy.constants import m_p as proton_mass


def molecular_frac(
    T: U.Quantity,
    rho: U.Quantity,
    mu: float = 1.22,
    EAGLE_corrections: bool = False,
    TNG_corrections: bool = False,
    Auriga_corrections: bool = False,
    SFR: U.Quantity | None = None,
    fNeutral: np.ndarray | None = None,
    gamma: float = 4.0 / 3.0,
    fH: float = 0.752,
    T0: U.Quantity = 8.0e3 * U.K,
) -> np.ndarray:
    """
    Compute particle molecular hydrogen fractions.

    To compute molecular mass of particle, multiply particle mass by hydrogen mass
    fraction and molecular fraction.

    All arguments should be passed with units as applicable, use :mod:`~astropy.units`.

    Parameters
    ----------
    T : ~astropy.units.Quantity
        Particle temperatures.
    rho : ~astropy.units.Quantity
        Particle densities.
    mu : float
        Mean molecular weight, default 1.22.
    EAGLE_corrections : bool
        Determine which particles are on the EoS and adjust values accordingly.
    TNG_corrections : bool
        Determine which particles have density > .1cm^-3 and give them a neutral fraction
        of 1.
    Auriga_corrections : bool
        Adjust pressures of star-forming gas particles.
    SFR : ~astropy.units.Quantity
        Particle star formation rates (required with EAGLE_corrections &
        Auriga_corrections).
    fNeutral : ~numpy.ndarray
        Particle neutral fraction (required with Auriga_corrections).
    gamma : float
        Polytropic index, default 4/3.
    fH : float
         Primordial hydrogen abundance, default 0.752.
    T0 : ~astropy.units.Quantity
         EoS critical temperature, default 8000 K.

    Returns
    -------
    ~numpy.ndarray
        An array of the same shape as particle property inputs containing the molecular
        mass fractions.

    Notes
    -----
    Based on the partitioning scheme of:
    Blitz, L., & Rosolowski, E. 2006, ApJ, 650, 933.

    Kyle Oman c. December 2015, updated October 2017.
    """
    # cast to float64 to avoid underflow
    P = U.Quantity(rho * T / mu, dtype=np.float64) / proton_mass

    if EAGLE_corrections:
        SFR = U.quantity.Quantity(SFR, copy=True)
        rho0 = 0.1 * U.cm**-3 * proton_mass / fH
        rho0 = rho0.to(U.Msun * U.kpc**-3)  # avoid overflow
        # cast to float64 to avoid underflow
        P0 = U.Quantity(rho0 * T0 / mu, dtype=np.float64) / proton_mass
        P_jeans = P0 * np.power(rho / rho0, gamma)
        P_margin = np.log10(P / P_jeans)
        SFR[P_margin > 0.5] = 0
        return np.where(
            SFR > 0, 1.0 / (1.0 + np.power(P / (4.3e4 * U.K * U.cm**-3), -0.92)), 0.0
        )
    elif Auriga_corrections:
        assert SFR is not None
        P[SFR > 0] = (P * fNeutral)[SFR > 0]
        return 1.0 / (1.0 + np.power(P / (1.7e4 * U.K * U.cm**-3), -0.8))
    else:
        return 1.0 / (1.0 + np.power(P / (4.3e4 * U.K * U.cm**-3), -0.92))
