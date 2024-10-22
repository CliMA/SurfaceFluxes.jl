var documenterSearchIndex = {"docs":
[{"location":"UniversalFunctions/#Universal-Functions","page":"Universal Functions","title":"Universal Functions","text":"","category":"section"},{"location":"UniversalFunctions/","page":"Universal Functions","title":"Universal Functions","text":"UniversalFunctions.jl provides universal functions for SurfaceFluxes.jl. Here, we reproduce some plots from literature, specifically from Gryanik et al. 2020, and Businger.","category":"page"},{"location":"UniversalFunctions/","page":"Universal Functions","title":"Universal Functions","text":"include(\"plot_universal_functions.jl\")","category":"page"},{"location":"UniversalFunctions/#Figs-1,2-(Gryanik)","page":"Universal Functions","title":"Figs 1,2 (Gryanik)","text":"","category":"section"},{"location":"UniversalFunctions/","page":"Universal Functions","title":"Universal Functions","text":"(Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"UniversalFunctions/#Fig-3-(Gryanik)","page":"Universal Functions","title":"Fig 3 (Gryanik)","text":"","category":"section"},{"location":"UniversalFunctions/","page":"Universal Functions","title":"Universal Functions","text":"(Image: ) (Image: )","category":"page"},{"location":"UniversalFunctions/#Figs-1,2-(Businger)","page":"Universal Functions","title":"Figs 1,2 (Businger)","text":"","category":"section"},{"location":"UniversalFunctions/","page":"Universal Functions","title":"Universal Functions","text":"(Image: ) (Image: )","category":"page"},{"location":"References/#References","page":"References","title":"References","text":"","category":"section"},{"location":"References/","page":"References","title":"References","text":"","category":"page"},{"location":"#SurfaceFluxes.jl","page":"Home","title":"SurfaceFluxes.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes","category":"page"},{"location":"#SurfaceFluxes","page":"Home","title":"SurfaceFluxes","text":"SurfaceFluxes\n\nInterface\n\nsurface_conditions computes\nMonin-Obukhov length\nPotential temperature flux (if not given) using Monin-Obukhov theory\ntransport fluxes using Monin-Obukhov theory\nfriction velocity/temperature scale/tracer scales\nexchange coefficients\n\nReferences\n\nS Nishizawa, Y Kitamura (2018)\nDaewon W. Byun (1990)\n\n\n\n\n\n","category":"module"},{"location":"#Core-input-types","page":"Home","title":"Core input types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.SurfaceValues\nSurfaceFluxes.InteriorValues","category":"page"},{"location":"#SurfaceFluxes.SurfaceValues","page":"Home","title":"SurfaceFluxes.SurfaceValues","text":"SurfaceValues\n\nInput container for state variables at the ground level.\n\nFields\n\nz\nu\nts\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.InteriorValues","page":"Home","title":"SurfaceFluxes.InteriorValues","text":"InteriorValues\n\nInput container for state variables at the first interior node.\n\nFields\n\nz\nu\nts\n\n\n\n\n\n","category":"type"},{"location":"#Dispatch-types","page":"Home","title":"Dispatch types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.Fluxes\nSurfaceFluxes.FluxesAndFrictionVelocity\nSurfaceFluxes.Coefficients\nSurfaceFluxes.ValuesOnly","category":"page"},{"location":"#SurfaceFluxes.Fluxes","page":"Home","title":"SurfaceFluxes.Fluxes","text":"Fluxes\n\nInput container for state variables, latent and sensible heat fluxes roughness lengths, initial obukhov length and gustiness.\n\nFields\n\nstate_in\nstate_sfc\nshf\nlhf\nz0m\nz0b\nL_MO_init\ngustiness\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.FluxesAndFrictionVelocity","page":"Home","title":"SurfaceFluxes.FluxesAndFrictionVelocity","text":"FluxesAndFrictionVelocity\n\nInput container, given surface state variables, latent and sensible heat fluxes, and the friction velocity, roughness lengths, initial obukhov length and gustiness.\n\nFields\n\nstate_in\nstate_sfc\nshf\nlhf\nustar\nz0m\nz0b\ngustiness\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.Coefficients","page":"Home","title":"SurfaceFluxes.Coefficients","text":"Coefficients\n\nInput container, given surface state variables, and exchange coefficients,roughness lengths, initial obukhov length and gustiness.\n\nFields\n\nstate_in\nstate_sfc\nCd\nCh\nz0m\nz0b\ngustiness\nbeta\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.ValuesOnly","page":"Home","title":"SurfaceFluxes.ValuesOnly","text":"ValuesOnly\n\nInput container, given only surface state variables, roughness lengths, initial obukhov length and gustiness.\n\nFields\n\nstate_in\nstate_sfc\nz0m\nz0b\nL_MO_init\ngustiness\nbeta\n\n\n\n\n\n","category":"type"},{"location":"#User-facing-methods","page":"Home","title":"User-facing methods","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.surface_conditions\nSurfaceFluxes.recover_profile","category":"page"},{"location":"#SurfaceFluxes.surface_conditions","page":"Home","title":"SurfaceFluxes.surface_conditions","text":"surface_conditions(\n    param_set::AbstractSurfaceFluxesParameters,\n    sc::SurfaceFluxes.AbstractSurfaceConditions{FT},\n    scheme::SurfaceFluxes.SolverScheme = FVScheme();\n    tol::RS.AbstractTolerance = RS.SolutionTolerance(FT(Δz(sc) / 50)),\n    tol_neutral::FT = SFP.cp_d(param_set) / 100,\n    maxiter::Int = 10,\n    soltype::RS.SolutionType = RS.CompactSolution(),\n) where {FT}\n\nThe main user facing function of the module. It computes the surface conditions based on the Monin-Obukhov similarity functions. Requires information about thermodynamic parameters (param_set) the surface state sc, the universal function type and the discretisation scheme. Default tolerance for  Monin-Obukhov length is absolute (i.e. has units [m]). Returns the RootSolvers CompactSolution by default.\n\nResult struct of type SurfaceFluxConditions{FT} contains:\n\nL_MO:   Monin-Obukhov lengthscale\nshf:    Sensible Heat Flux\nlhf:    Latent Heat Flux\nρτxz:   Momentum Flux (Eastward component)\nρτyz:   Momentum Flux (Northward component)\nustar:  Friction velocity\nCd:     Momentum Exchange Coefficient\nCh:     Thermal Exchange Coefficient\n\n\n\n\n\n","category":"function"},{"location":"#SurfaceFluxes.recover_profile","page":"Home","title":"SurfaceFluxes.recover_profile","text":"recover_profile(param_set, sc, L_MO, Z, X_in, X_sfc, transport, uft, scheme)\n\nRecover profiles of variable X given values of Z coordinates. Follows Nishizawa equation (21,22)\n\nArguments\n\nparam_set: Abstract Parameter Set containing physical, thermodynamic parameters.\nsc: Container for surface conditions based on known combination     of the state vector, and {fluxes, friction velocity, exchange coefficients} for a given experiment\nL_MO: Monin-Obukhov length\nZ: Z coordinate(s) (within surface layer) for which variable values are required\nXin,Xsfc: For variable X, values at interior and surface nodes\ntransport: Transport type, (e.g. Momentum or Heat, used to determine physical scale coefficients)\nuft: A Universal Function type, (returned by, e.g., Businger())\nscheme: Discretization scheme (currently supports FD and FV)\n\nTODO: add tests\n\n\n\n\n\n","category":"function"},{"location":"#Universal-Functions","page":"Home","title":"Universal Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.UniversalFunctions","category":"page"},{"location":"#SurfaceFluxes.UniversalFunctions","page":"Home","title":"SurfaceFluxes.UniversalFunctions","text":"UniversalFunctions\n\nUniversal stability and stability correction functions for SurfaceFluxes module. Supports universal functions:\n\nBusinger\nGryanik\nGrachev\n\n\n\n\n\n","category":"module"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.UniversalFunctions.Gryanik\nSurfaceFluxes.UniversalFunctions.Grachev\nSurfaceFluxes.UniversalFunctions.Businger","category":"page"},{"location":"#SurfaceFluxes.UniversalFunctions.Gryanik","page":"Home","title":"SurfaceFluxes.UniversalFunctions.Gryanik","text":"Gryanik <: AbstractUniversalFunction{FT}\n\nReferences\n\nVladimir M Gryanik, Christof L{\\\"u}pkes, Andrey Grachev, Dmitry Sidorenko (2020)\n\nEquations in reference:\n\n`ϕ_m`: Eq. 13\n`ϕ_h`: Eq. 13\n`ψ_m`: Eq. 14\n`ψ_h`: Eq. 14\n\nGryanik et al. (2020) functions are used in stable conditions\n\nIn unstable conditions the functions of Businger (1971) are\n\nassigned by default.\n\nFields\n\nL\n\n: Monin-Obhukov Length\n\nparams\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.UniversalFunctions.Grachev","page":"Home","title":"SurfaceFluxes.UniversalFunctions.Grachev","text":"Grachev <: AbstractUniversalFunction{FT}\n\nReferences\n\nAndrey A Grachev, Edgar L Andreas, Christopher W Fairall, Peter S Guest, P Ola G Persson (2007)\n\nEquations in reference:\n\n`ϕ_m`: Eq. 13\n`ϕ_h`: Eq. 13\n`ψ_m`: Eq. 14\n`ψ_h`: Eq. 14\n\nGrachev (2007) functions are applicable in the\n\nstable b.l. regime (ζ >= 0). Businger (1971) functions\n\nare applied in the unstable b.l. (ζ<0) regime by\n\ndefault.\n\nFields\n\nL\n\n: Monin-Obhukov Length\n\nparams\n\n\n\n\n\n","category":"type"},{"location":"#SurfaceFluxes.UniversalFunctions.Businger","page":"Home","title":"SurfaceFluxes.UniversalFunctions.Businger","text":"Businger\n\nReference\n\nS Nishizawa, Y Kitamura (2018)\n\nOriginal research\n\nJoost A Businger, John C Wyngaard, Yꎬ Izumi, Edward F Bradley (1971)\n\nEquations in reference:\n\n`ϕ_m`: Eq. A1\n`ϕ_h`: Eq. A2\n`ψ_m`: Eq. A3\n`ψ_h`: Eq. A4\n\nFields\n\nL\n\n: Monin-Obhukov Length\n\nparams\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"SurfaceFluxes.UniversalFunctions.phi\nSurfaceFluxes.UniversalFunctions.psi\nSurfaceFluxes.UniversalFunctions.Psi","category":"page"},{"location":"#SurfaceFluxes.UniversalFunctions.phi","page":"Home","title":"SurfaceFluxes.UniversalFunctions.phi","text":"phi\n\nUniversal stability function for wind shear (ϕ_m) and temperature gradient (ϕ_h)\n\n\n\n\n\n","category":"function"},{"location":"#SurfaceFluxes.UniversalFunctions.psi","page":"Home","title":"SurfaceFluxes.UniversalFunctions.psi","text":"psi\n\nUniversal stability correction function for momentum (ψ_m) and heat (ψ_h)\n\n\n\n\n\n","category":"function"},{"location":"#SurfaceFluxes.UniversalFunctions.Psi","page":"Home","title":"SurfaceFluxes.UniversalFunctions.Psi","text":"Psi\n\nIntegral of universal stability correction function for momentum (ψ_m) and heat (ψ_h)\n\n\n\n\n\n","category":"function"}]
}
