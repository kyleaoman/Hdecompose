"""Provide convenience method to compute the atomic mass fraction of particles."""

import numpy as np
from astropy import units as U
from .BlitzRosolowsky2006 import molecular_frac as calc_molecular_frac
from .RahmatiEtal2013 import neutral_frac as calc_neutral_frac


def atomic_frac(
    redshift: float | None = None,
    nH: U.Quantity | None = None,
    T: U.Quantity | None = None,
    rho: U.Quantity | None = None,
    Habundance: np.ndarray | None = None,
    onlyA1: bool = False,
    noCol: bool = False,
    onlyCol: bool = False,
    SSH_Thresh: bool = False,
    local: bool = False,
    EAGLE_corrections: bool = False,
    TNG_corrections: bool = False,
    SFR: U.Quantity = None,
    mu: float = 1.22,
    gamma: float = 4.0 / 3.0,
    fH: float = 0.752,
    T0: U.Quantity = 8.0e3 * U.K,
    neutral_frac: np.ndarray | None = None,
    molecular_frac: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute particle atomic hydrogen mass fractions.

    All arguments should be passed with units as applicable, use :mod:`~astropy.units`.

    Parameters
    ----------
    redshift : float
        Snapshot redshift.
    nH : ~astropy.units.Quantity
        Hydrogen number density of the gas.
    T : ~astropy.units.Quantity
        Temperature of the gas.
    rho : ~astropy.units.Quantity
        Gas particle density.
    Habundance : ~numpy.ndarray
        Particle Hydrogen mass fractions.
    onlyA1 : bool
        Routine will use Table A1 parameters for z < 0.5.
    noCol : bool
        The contribution of collisional ionisation to the overall ionisation rate is
        neglected.
    onlyCol : bool
        The contribution of photoionisation to the overall ionisation rate is neglected.
    SSH_Thresh : ~astropy.units.Quantity
        All particles above this density are assumed to be fully shielded, i.e.
        f_neutral=1.
    local : bool
        Compute the local polytropic index.
    EAGLE_corrections : bool
        Determine which particles are on the EoS and adjust values accordingly.
    TNG_corrections : bool
        Determine which particles have density > .1cm^-3 and give them a neutral fraction
        of 1.
    SFR : ~astropy.units.Quantity
        Particle star formation rates.
    mu : float
        Mean molecular weight, default 1.22 (required with EAGLE_corrections).
    gamma : float
        Polytropic index, default 4/3 (required with EAGLE_corrections).
    fH : float
        Primordial hydrogen abundance, default 0.752 (required with EAGLE_corrections).
    T0 : ~astropy.units.Quantity
        EoS critical temperature, default 8000 K (required with EAGLE_corrections).
    neutral_frac : ~numpy.ndarray
        Previously computed neutral fractions can be provided. In this case can omit
        redshift, nH, onlyA1, noCol, onlyCol, SSH_Thresh, local, TNG_corrections,
        Habundance - will be ignored.
    molecular_frac : ~numpy.ndarray
        Previously computed molecular fractions can be provided.

    Returns
    -------
    ~numpy.ndarray
        An array of the same shape as particle property inputs containing the atomic (HI)
        mass fractions.

    See Also
    --------
    Hdecompose.BlitzRosolowsky2007.molecular_frac
    Hdecompose.RahmatiEtal2013.neutral_frac
    Hdecompose.SpringelHernquist2003.sf_neutral_frac
    """
    assert redshift is not None
    if neutral_frac is None:
        neutral_frac = calc_neutral_frac(
            redshift,
            nH,
            T,
            onlyA1=onlyA1,
            noCol=noCol,
            onlyCol=onlyCol,
            SSH_Thresh=SSH_Thresh,
            local=local,
            EAGLE_corrections=EAGLE_corrections,
            TNG_corrections=TNG_corrections,
            SFR=SFR,
            mu=mu,
            gamma=gamma,
            fH=fH,
            Habundance=Habundance,
            T0=T0,
            rho=rho,
        )
    if molecular_frac is None:
        molecular_frac = calc_molecular_frac(
            T,
            rho,
            EAGLE_corrections=EAGLE_corrections,
            SFR=SFR,
            mu=mu,
            gamma=gamma,
            fH=fH,
            T0=T0,
        )
    return (1.0 - molecular_frac) * neutral_frac
