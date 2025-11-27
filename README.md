# The SANSND plugin for  NCrystal

This repository contains code for the NCrystal plugin named `SANSND`, which implement Small Angles Neutron Scattering for NanoDiamonds.

This plugin is mainly based on the work of José Ignacio Marquez Damian, Monte Carlo Simulation Scientist at ESS Spallation Physics, who implemented Nanodiamond support in MCNP and PHITS.

**It is currently under development and not ready for general usage**

The simple model implemented here has its basis in the SANS structure factor measured by [Teshigawara et al.](https://doi.org/10.1016/j.nima.2019.03.038) and combines it with a scattering kernel for bulk diamond computed from DFT as explained by [Granada et al.](https://doi.org/10.1051/epjconf/202023104002). The structure factor in the Teshigawara paper is fitted by a power-exponential law, which is then further simplified into a piecewise power fit.

### Minimum Scattering Angle

The plugin supports an optional minimum scattering angle, specified in the `.ncmat` file as `theta_min_deg`. This parameter, measured in degrees, sets a lower cutoff for the scattering angle. When `theta_min_deg` is provided, any scattering events with an angle smaller than this value are not considered "removed" from the beam. This affects both the SANS removal cross section and the sampling of scattering angles.

If `theta_min_deg` is not specified, it defaults to 0.0, meaning all scattering angles are considered.

