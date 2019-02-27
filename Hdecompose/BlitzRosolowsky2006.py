import numpy as np
from astropy import units as U
from astropy.constants import m_p as proton_mass


def molecular_frac(
        T,
        rho,
        mu=1.22,
        EAGLE_corrections=False,
        TNG_corrections=False,
        SFR=None,
        gamma=4./3.,
        fH=0.752,
        T0=8.0E3*U.K
):
    """
    Computes particle molecular hydrogen fractions.

    To compute molecular mass of particle, multiply particle mass by hydrogen
    mass fraction and molecular fraction.

    All arguments should be passed with units as applicable, use astropy.units.

    T:                 Particle temperatures.
    rho:               Particle densities.
    mu:                Mean molecular weight, default 1.22.
    EAGLE_corrections: Determine which particles are on the EoS and adjust
                       values accordingly.
    SFR:               Particle star formation rates (required with
                       EAGLE_corrections).
    gamma:             Polytropic index, default 4/3.
    fH:                Primordial hydrogen abundance, default 0.752.
    T0:                EoS critical temperature, default 8000 K.

    Returns an array of the same shape as particle property inputs containing
    the molecular mass fractions.

    Based on the partitioning scheme of:
    Blitz, L., & Rosolowski, E. 2006, ApJ, 650, 933.

    Kyle Oman c. December 2015, updated October 2017.
    """

    # cast to float64 to avoid underflow
    P = U.Quantity(rho * T / mu, dtype=np.float64) / proton_mass

    if EAGLE_corrections and TNG_corrections:
        raise ValueError

    if EAGLE_corrections or TNG_corrections:
        SFR = U.quantity.Quantity(SFR, copy=True)
        if EAGLE_corrections:
            rho0 = 0.1 * U.cm ** -3 * proton_mass / fH
        elif TNG_corrections:
            rho0 = 0.13 * U.cm ** -3 * proton_mass / fH
        rho0 = rho0.to(U.Msun * U.kpc ** -3)  # avoid overflow
        # cast to float64 to avoid underflow
        P0 = U.Quantity(rho0 * T0 / mu, dtype=np.float64) / proton_mass
        P_jeans = P0 * np.power(rho / rho0, gamma)
        P_margin = np.log10(P / P_jeans)
        SFR[P_margin > .5] = 0
        return np.where(
            SFR > 0,
            1. / (1. + np.power(P / (4.3E4 * U.K * U.cm ** -3), -.92)),
            0.
        )
    else:
        return 1. / (1. + np.power(P / (4.3E4 * U.K * U.cm ** -3), -.92))
