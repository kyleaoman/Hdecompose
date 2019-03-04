from astropy import units as U
from .BlitzRosolowsky2006 import molecular_frac
from .RahmatiEtal2013 import neutral_frac


def atomic_frac(
        redshift,
        nH,
        T,
        rho,
        Habundance,
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
        T0=8.E3 * U.K,
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
    TNG_correctiosn:   Determine which particles are colder than 10^4K and give
                       them a neutral fraction of 1.
    mu:                Mean molecular weight, default 1.22 (required with
                       EAGLE_corrections).
    gamma:             Polytropic index, default 4/3 (required with
                       EAGLE_corrections).
    fH:                Primordial hydrogen abundance, default 0.752 (required
                       with EAGLE_corrections).
    T0:                EoS critical temperature, default 8000 K (required with
                       EAGLE_corrections).

    Returns an array of the same shape as particle property inputs containing
    the atomic (HI) mass fractions.
    """

    return (1. - molecular_frac(
        T,
        rho,
        EAGLE_corrections=EAGLE_corrections,
        SFR=SFR,
        mu=mu,
        gamma=gamma,
        fH=fH,
        T0=T0
    )) * \
        neutral_frac(
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
            rho=rho
        )
