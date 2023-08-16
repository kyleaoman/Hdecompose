from astropy import units as U
from .BlitzRosolowsky2006 import molecular_frac as calc_molecular_frac
from .RahmatiEtal2013 import neutral_frac as calc_neutral_frac


def atomic_frac(
    redshift=None,
    nH=None,
    T=None,
    rho=None,
    Habundance=None,
    onlyA1=False,
    noCol=False,
    onlyCol=False,
    SSH_Thresh=False,
    local=False,
    EAGLE_corrections=False,
    TNG_corrections=False,
    SFR=None,
    mu=1.22,
    gamma=4.0 / 3.0,
    fH=0.752,
    T0=8.0e3 * U.K,
    neutral_frac=None,
    molecular_frac=None,
):
    """
    Computes particle atomic hydrogen mass fractions. See also molecular_frac
    and neutral_frac in this module.

    All arguments should be passed with units as applicable, use astropy.units.

    redshift:          Snapshot redshift.
    nH:                Hydrogen number density of the gas.
    T:                 Temperature of the gas.
    SFR:               Particle star formation rates.
    rho:               Gas particle density.
    Habundance:        Particle Hydrogen mass fractions.
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
    mu:                Mean molecular weight, default 1.22 (required with
                       EAGLE_corrections).
    gamma:             Polytropic index, default 4/3 (required with
                       EAGLE_corrections).
    fH:                Primordial hydrogen abundance, default 0.752 (required
                       with EAGLE_corrections).
    T0:                EoS critical temperature, default 8000 K (required with
                       EAGLE_corrections).
    neutral_frac:      Previously computed neutral fractions can be provided.
                       In this case can omit redshift, nH, onlyA1, noCol,
                       onlyCol, SSH_Thresh, local, TNG_corrections, Habundance
                       - will be ignored.
    molecular_frac:    Previously computed molecular fractions can be provided.

    Returns an array of the same shape as particle property inputs containing
    the atomic (HI) mass fractions.
    """

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
