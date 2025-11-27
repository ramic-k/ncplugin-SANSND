# The SANSND plugin for  NCrystal

This repository contains code for the NCrystal plugin named `SANSND`, which implement Small Angles Neutron Scattering for NanoDiamonds.

This plugin is mainly based on the work of José Ignacio Marquez Damian, Monte Carlo Simulation Scientist at ESS Spallation Physics, who implemented Nanodiamond support in MCNP and PHITS.

**It is currently under development and not ready for general usage**

The simple model implemented here has its basis in the SANS structure factor measured by [Teshigawara et al.](https://doi.org/10.1016/j.nima.2019.03.038) and combines it with a scattering kernel for bulk diamond computed from DFT as explained by [Granada et al.](https://doi.org/10.1051/epjconf/202023104002). The structure factor in the Teshigawara paper is fitted by a power-exponential law, which is then further simplified into a piecewise power fit.

## theta_min_deg option

You can optionally set a minimum scattering angle (in degrees) inside the `@CUSTOM_SANSND` block via an extra line:

```
theta_min_deg 0.25
```

- Default is `0.0`, matching the previous full angular coverage.
- For `theta_min_deg>0`, the SANS removal cross section integrates only q values above `q_min = 2k sin(theta_min/2)`.
- Scattering samples are drawn from the same truncated distribution, so no angles below `theta_min_deg` are emitted.
