import numpy as np
from astropy import units as U
from astropy.constants import m_p as proton_mass


def neutral_frac(
        redshift,
        nH,
        T,
        onlyA1=False,
        noCol=False,
        onlyCol=False,
        SSH_Thresh=False,
        local=False,
        EAGLE_corrections=False,
        TNG_corrections=False,
        SFR=None,
        mu=1.22,
        gamma=4./3.,
        fH=0.752,
        Habundance=None,
        T0=8.0E3 * U.K,
        rho=None
):

    """
    Computes particle neutral hydrogen fractions based on the fitting functions
    of:
    Rahmati, A., Pawlik, A. H., Raicevic, M., & Schaye, J. 2013, MNRAS, 430,
    2427.
    By default, it uses the parameters of Table A1 (based on small cosmological
    volumes) for z > 1, and of Table A2 (based on a 50Mpc volume) for z < 1, to
    better account for the effects of collisional ionisation on the
    self-shielding density.

    To compute neutral (HI + H_2) mass of particle, multiply NeutralFraction by
    Hydrogen mass fraction and particle mass.

    All arguments should be passed with units as applicable, use astropy.units.

    redshift:          Snapshot redshift.
    nH:                Hydrogen number density of the gas.
    T:                 Temperature of the gas.

    onlyA1:            Routine will use Table A1 parameters for z < 0.5.
    noCol:             The contribution of collisional ionisation to the
                       overall ionisation rate is neglected.
    onlyCol:           The contribution of photoionisation to the overall
                       ionisation rate is neglected.
    SSH_Thresh:        All particles above this density are assumed to be fully
                       shielded, i.e. f_neutral=1.
    local:             Compute the local polytropic index.
    EAGLE_corrections: Determine which particles are on the EoS and adjust
                       values accordingly.
    TNG_corrections:   Determine which particles have density > .1cm^-3 and
                       give them a neutral fraction of 1.
    SFR:               Particle star formation rates (required with
                       EAGLE_corrections).
    mu:                Mean molecular weight, default 1.22 (required with
                       EAGLE_corrections).
    gamma:             Polytropic index, default 4/3 (required with
                       EAGLE_corrections).
    fH:                Primordial hydrogen abundance, default 0.752 (required
                       with EAGLE_corrections).
    Habundance:        Particle Hydrogen mass fractions (required with
                       EAGLE_corrections).
    T0:                EoS critical temperature, default 8000 K (required with
                       EAGLE_corrections).
    rho:               Gas particle density (required with EAGLE_corrections).

    Returns an array of the same shape as particle property inputs containing
    the neutral mass fractions.

    Original IDL function written by Rob Crain, Leiden, March 2014, with input
    from Ali Rahmati. Based on Rahmati et al. (2013). Translated from IDL to
    Python 2.7 by Kyle Oman, Victoria, December 2015, updated October 2017 for
    Python 3.
    """

    if EAGLE_corrections and TNG_corrections:
        raise ValueError

    # EAGLE pre-treatment for gas temperature
    if EAGLE_corrections:
        T = U.quantity.Quantity(T, copy=True)
        SFR = U.quantity.Quantity(SFR, copy=True)
        # cast to float64 to avoid underflow
        P = U.Quantity(rho * T / mu, dtype=np.float64) / proton_mass
        rho0 = 0.1 * U.cm ** -3 * proton_mass / fH
        rho0 = rho0.to(U.Msun * U.kpc ** -3)  # avoid underflow
        P0 = U.Quantity(rho0 * T0 / mu, dtype=np.float64) / proton_mass
        P_jeans = P0 * np.power(rho / rho0, gamma)
        P_margin = np.log10(P / P_jeans)
        SFR[P_margin > .5] = 0
        T_jeans = mu * Habundance * P_jeans / nH
        T[
            np.logical_or(
                SFR > 0,
                np.logical_and(
                    P_margin < .5,
                    T_jeans > 1.E4 * U.K
                )
            )
        ] = 1.E4 * U.K

    if ((redshift >= 0.0) and (redshift < 1.0)) or \
       np.isclose(redshift, 0.0, atol=1.E-3):
        dz = redshift
        if onlyA1:
            lg_n0_lo = -2.94
            gamma_uvb_lo = 8.34E-14
            alpha1_lo = -3.98
            alpha2_lo = -1.09
            beta_lo = 1.29
            f_lo = 0.01
        else:
            lg_n0_lo = -2.56
            gamma_uvb_lo = 8.34E-14
            alpha1_lo = -1.86
            alpha2_lo = -0.51
            beta_lo = 2.83
            f_lo = 0.01

        lg_n0_hi = -2.29
        gamma_uvb_hi = 7.3E-14
        alpha1_hi = -2.94
        alpha2_hi = -0.90
        beta_hi = 1.21
        f_hi = 0.03

    elif (redshift >= 1.0) and (redshift < 2.0):
        dz = redshift - 1.0

        lg_n0_lo = -2.29
        gamma_uvb_lo = 7.3E-14
        alpha1_lo = -2.94
        alpha2_lo = -0.90
        beta_lo = 1.21
        f_lo = 0.03

        lg_n0_hi = -2.06
        gamma_uvb_hi = 1.50E-12
        alpha1_hi = -2.22
        alpha2_hi = -1.09
        beta_hi = 1.75
        f_hi = 0.03

    elif (redshift >= 2.0) and (redshift < 3.0):
        dz = redshift - 2.0

        lg_n0_lo = -2.06
        gamma_uvb_lo = 1.50E-12
        alpha1_lo = -2.22
        alpha2_lo = -1.09
        beta_lo = 1.75
        f_lo = 0.03

        lg_n0_hi = -2.13
        gamma_uvb_hi = 1.16E-12
        alpha1_hi = -1.99
        alpha2_hi = -0.88
        beta_hi = 1.72
        f_hi = 0.04

    elif (redshift >= 3.0) and (redshift < 4.0):
        dz = redshift - 3.0

        lg_n0_lo = -2.13
        gamma_uvb_lo = 1.16E-12
        alpha1_lo = -1.99
        alpha2_lo = -0.88
        beta_lo = 1.72
        f_lo = 0.04

        lg_n0_hi = -2.23
        gamma_uvb_hi = 7.91E-13
        alpha1_hi = -2.05
        alpha2_hi = -0.75
        beta_hi = 1.93
        f_hi = 0.02

    elif (redshift >= 4.0) and (redshift < 5.0):
        dz = redshift - 4.0

        lg_n0_lo = -2.23
        gamma_uvb_lo = 7.91E-13
        alpha1_lo = -2.05
        alpha2_lo = -0.75
        beta_lo = 1.93
        f_lo = 0.02

        lg_n0_hi = -2.35
        gamma_uvb_hi = 5.43E-13
        alpha1_hi = -2.63
        alpha2_hi = -0.57
        beta_hi = 1.77
        f_hi = 0.01

    else:
        print("Invalid redshift > 5.0 or < 0.0")
        print("Redshift is: {:.5e}".format(redshift))
        raise ValueError

    lg_n0 = lg_n0_lo + dz * (lg_n0_hi - lg_n0_lo)
    n0 = np.power(10, lg_n0) * U.cm ** -3
    gamma_uvb = gamma_uvb_lo + dz * (gamma_uvb_hi - gamma_uvb_lo)
    alpha1 = alpha1_lo + dz * (alpha1_hi - alpha1_lo)
    alpha2 = alpha2_lo + dz * (alpha2_hi - alpha2_lo)
    beta = beta_lo + dz * (beta_hi - beta_lo)
    f = f_lo + dz * (f_hi - f_lo)

    gamma_ratio = (1.0 - f) * np.power(1.0 + np.power(nH / n0, beta), alpha1)
    gamma_ratio = gamma_ratio + f * np.power(1.0 + nH / n0, alpha2)
    gamma_phot = gamma_uvb * gamma_ratio

    if local:
        gamma_local = 1.3E-13 * np.power(nH, 0.2) * np.power(T / 1.0E4, 0.2)
        gamma_phot = gamma_phot + gamma_local

    Lambda = 315614.0 * U.K / T
    AlphaA = 1.269E-13 * np.power(Lambda, 1.503)
    AlphaA = AlphaA / np.power(1.0 + np.power(Lambda / 0.522, 0.470), 1.923)
    LambdaT = 1.17E-10 * np.sqrt(T / U.K) * np.exp(-157809.0 * U.K / T) / \
        (1.0 + np.sqrt(T / (1.0E5 * U.K)))

    if noCol:
        LambdaT = 0.0
    if onlyCol:
        gamma_phot = 0.0

    A = AlphaA + LambdaT
    B = 2.0 * AlphaA + (gamma_phot * U.cm ** -3 / nH) + LambdaT
    sqrt_arg = np.power(B, 2) - 4.0 * A * AlphaA
    sqrt_arg[sqrt_arg < 0.0] = 0.0
    sqrt_term = np.sqrt(sqrt_arg)
    f_neutral = (B - sqrt_term) / (2.0 * A)

    if SSH_Thresh:
        f_neutral[nH > SSH_Thresh] = 1.0

    if TNG_corrections:
        f_neutral[nH > 0.1 * U.cm ** -3] = 1.

    return f_neutral
