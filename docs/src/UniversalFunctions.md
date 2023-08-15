# Universal Functions

UniversalFunctions.jl provides universal functions for SurfaceFluxes.jl. The functions are defined in `src/UniversalFunctions.jl`, while the plots are generated in `docs/src/plot_universal_functions.jl` and saved in `docs/src/assets`. Here, we reproduce some plots from literature, specifically from Gryanik et al. 2020, Businger 1971, and Bonan 2019. Note that Bonan uses the forms of $\phi$ and $\psi$ from Dyer and Hicks 1970; Dyer 1974; Brutsaert 1982, pp. 68–71; Garratt 1992, pp. 52–54. 

```@example
include("plot_universal_functions.jl")
```

## Figs 1,2 (Gryanik)

![](Gryanik12_phi_h.svg)
![](Gryanik12_phi_m.svg)
![](Gryanik12_psi_h.svg)
![](Gryanik12_psi_m.svg)

## Fig 3 (Gryanik)

![](Gryanik3_phi_h.svg)
![](Gryanik3_phi_m.svg)


## Figs 1,2 (Businger)

![](Businger_phi_h.svg)
![](Businger_phi_m.svg)


## Figs 1,2 (Bonan)

![](Bonan_phi_h.svg)
![](Bonan_phi_m.svg)
![](Bonan_psi_h.svg)
![](Bonan_psi_m.svg)
