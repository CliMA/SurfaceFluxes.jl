[PR212] Refactor of SurfaceFluxes.jl: Consistently use stability parameter in all solvers and as inputs to many functions. Added functionality for wind speed dependent roughness lengths (Charnock, COARE3) and option to use functions to compute surface temperature/humidity. 

[PR 206] Updates UniversalFunctions.jl: Update functions to avoid catastrophic cancellations. Add continuity and linearisation tests + fix bug in the near-neutral limit. 

[PR 186] Removes unused stability function types (Holtslag, Cheng, Beljaars). Currently supports (Businger, Grachev, Gryanik) types.  
