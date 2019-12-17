from astropy import units as U
from astropy.constants import m_p, k_B


def auriga_correct_neutral_frac(
        fNeutral,
        SFR,
        u,
        mu=1.22,
        gamma=4 / 3,
):

    """
    Computes particle neutral hydrogen fractions based on the multiphase ISM
    model of:
    Springel, V. and Hernquist, L. 2013, MNRAS, 339, 289.

    To compute neutral (HI + H_2) mass of particle, multiply NeutralFraction by
    Hydrogen mass fraction and particle mass.

    All arguments should be passed with units as applicable, use astropy.units.

    fNeutral:          Gas neutral fractions from ionization model.
    SFR:               Gas star formation rate.
    u:                 Gas specific internal energy.
    mu:                Mean molecular weight (default 1.22).
    gamma:             Polytropic index, default 4/3.

    Returns an array of the same shape as particle property inputs containing
    the neutral mass fractions.
    """

    Th = 1E6 * U.K  # guided by Springel & Hernquist 2003
    Tc = 1E3 * U.K  # guided by Springel & Hernquist 2003
    retval = U.Quantity.copy(fNeutral)
    mask = SFR > 0
    uh = (k_B * Th / (gamma - 1) / mu / m_p).to(U.km ** 2 / U.s ** 2)
    uc = (k_B * Tc / (gamma - 1) / mu / m_p).to(U.km ** 2 / U.s ** 2)
    retval[mask] = ((uh - u) / (uh - uc)).to(U.dimensionless_unscaled)
    return retval
