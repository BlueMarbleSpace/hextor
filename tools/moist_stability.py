#!/usr/bin/env python3
"""
moist_stability.py — the explicit-diffusion timestep limit, with moist static
energy diffusion taken into account.

HEXTOR integrates the diffusion term explicitly on 18 belts in x = sin(lat),
so it needs D*dt/(C*dx^2) below about 1/2.  Exceeding it raises no error: the
model integrates quietly to NaN over thousands of years, which the SAMOSA
harness would then label a runaway.

With &ebm::moistdiff the diffused quantity is h/cp = T + (L/cp) q rather than
T, and h/cp responds to a temperature perturbation faster than T does, by

    beta = d(h/cp)/dT = 1 + (L/cp) dq/dT,

so the same D carries beta times more heat per kelvin and the stable timestep
falls by beta.  beta is 2.7 at 290 K and 13 at 343 K for an N2 atmosphere at
1 bar, so it is not a small correction anywhere the surface is warm.

The functions below mirror driver.f (moistprops, esatw, qmoist) so that the
timestep the harness picks is the one the Fortran actually needs.
"""

XLV = 2.501e6            # latent heat of vaporization [J/kg], driver.f
NBELTS = 18
CFL_SAFETY = 0.4         # explicit diffusion needs D*dt/(C*dx^2) below ~1/2


def moistprops(pn2, pco2, ph2=0.0, pch4=0.0):
    """Dry-air specific heat cpd [J/kg/K] and eps = Mw/Md, as driver.f does it.

    Partial pressures in bar; only their ratios matter.
    """
    wn2, wco2, wh2, wch4 = 28.0 * pn2, 44.0 * pco2, 2.0 * ph2, 16.0 * pch4
    wsum = wn2 + wco2 + wh2 + wch4
    cpd = 4184. * (0.2484 * wn2 + 0.2105 * wco2 + 3.420 * wh2
                   + 0.5271 * wch4) / wsum
    eps = 18.016 * (pn2 + pco2 + ph2 + pch4) / wsum
    return cpd, eps


def esat(t):
    """Saturation vapour pressure of water [Pa], matching driver.f's esatw:
    Wagner & Pruss (2002) over liquid, Clausius-Clapeyron over ice."""
    import math
    if t < 273.16:
        tuse = max(t, 50.)
        return 611.2 * math.exp((2.8347e6 / 461.50) * (1. / 273.16 - 1. / tuse))
    tuse = min(t, 647.096)
    th = 1. - tuse / 647.096
    return 22.064e6 * math.exp((647.096 / tuse) * (
        -7.85951783 * th + 1.84408259 * th ** 1.5 - 11.7866497 * th ** 3
        + 22.6807411 * th ** 3.5 - 15.9618719 * th ** 4
        + 1.80122502 * th ** 7.5))


def qmoist(t, pdry, eps, rh):
    """Specific humidity [kg/kg] over a dry surface pressure pdry [bar]."""
    r = eps * rh * esat(t) / (pdry * 1.e5)
    return r / (1. + r)


def beta(t, pdry, eps, cpd, rh, dt=0.05):
    """d(h/cp)/dT at temperature t, by centred difference on driver.f's own q.

    Differencing the same expressions the Fortran evaluates is safer than an
    analytic Clausius-Clapeyron derivative, which would not match the
    Wagner-Pruss branch above the triple point.
    """
    if rh <= 0.0:
        return 1.0
    dq = (qmoist(t + dt, pdry, eps, rh) - qmoist(t - dt, pdry, eps, rh)) / (2 * dt)
    return 1.0 + (XLV / cpd) * dq


def stable_dt(D, heatcap, dt_cap, nbelts=NBELTS, cfl=CFL_SAFETY, beta_max=1.0):
    """Timestep satisfying the explicit-diffusion limit, capped at dt_cap.

    beta_max is the largest d(h/cp)/dT on the grid (1 for dry diffusion).
    """
    if D <= 0:
        return dt_cap
    dx = 2.0 / nbelts
    return min(dt_cap, cfl * heatcap * dx * dx / (D * max(beta_max, 1.0)))


def moist_dt(D, pdry, pn2, pco2, t_hot, rh, heatcap, dt_cap, **kw):
    """stable_dt with beta evaluated at the hottest belt, t_hot."""
    cpd, eps = moistprops(pn2, pco2)
    b = beta(t_hot, pdry, eps, cpd, rh)
    return stable_dt(D, heatcap, dt_cap, beta_max=b, **kw), b


if __name__ == '__main__':
    cpd, eps = moistprops(1.0, 4.0e-4)
    print('N2 + 400 ubar CO2:  cp_dry = %.1f J/kg/K   eps = %.4f' % (cpd, eps))
    print('%8s %10s %10s %12s' % ('T [K]', 'q [g/kg]', 'beta', 'dt_max(D=3.3)'))
    for t in (200, 240, 273, 290, 300, 310, 320, 343, 360, 400, 420):
        b = beta(t, 1.0, eps, cpd, 0.8)
        print('%8.0f %10.3f %10.2f %12.0f'
              % (t, 1e3 * qmoist(t, 1.0, eps, 0.8), b,
                 stable_dt(3.3, 4.0e6, 1350., beta_max=b)))
