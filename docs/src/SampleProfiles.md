# Sample Vertical Profiles of Wind Speed and Temperature

SurfaceFluxes.jl provides profile recovery functions for the roughness sublayer using the Physick and Garratt (1995) formulation. Here, we use these functions to reproduce vertical profiles of wind speed and temperature for different conditions of atmospheric stability, reproducing Bonan 2019 Figure 6.4 with and without canopy correction.

```@example
include("plot_profiles.jl")
```

# Fig 6.4 (a) 

![](Fig6.4a_profile.svg)
![](Fig6.4a_canopy_profile.svg)

# Fig 6.4 (b) 

![](Fig6.4b_profile.svg)
![](Fig6.4b_canopy_profile.svg)

# Fig 6.4 (c)

![](Fig6.4c_profile.svg)
![](Fig6.4c_canopy_profile.svg)

# Fig 6.4 (d) 

![](Fig6.4d_profile.svg)
![](Fig6.4d_canopy_profile.svg)

