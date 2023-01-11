## Surface Fluxes

### Monin-Obukhov Similarity Theory (MOST)
Surface fluxes of momentum, energy and moisture are computed using the Monin-Obukhov Similarity Theory (MOST), which follows from dimensional analysis and similarity laws in stratified boundary layers. The theory gives the fluxes in an unresolved surface layer based on the self-similarity of stratified boundary layers, and can be used to compute exchange coefficients for bulk-formulae for the fluxes at the bottom boundary in simulations of the atmosphere. The MOST similarity profiles depend on one length scale, given by the ratio of velocity and buoyancy scales and known as the Obukhov length:
```math
\begin{equation}
\label{eq:monin_obukhov_length}
L_O =  \frac{u_{\star}^2 }{\kappa b_{\star}}.
\end{equation}
```
Here, ${\star}$ subscripts indicate a characteristic physical scale of the variable in question in the surface layer. The buoyancy ($b$) is defined as
```math
\begin{equation}
b = g \frac{\mathrm{DSE}_v'}{\mathrm{DSE}_v},
\end{equation}
```
with $\mathrm{DSE}_v'$ as the perturbation from a reference virtual dry static energy value, $\mathrm{DSE}$. These mean and perturbation quantities are implicit in the model and can be used to compute the relation between $b_{\star}$ and $\phi_{\star}$, where $\phi$ represents a thermoynamic scale variable. Following equations (8) and (9) in \citep{Nishizawa18a} we write the physical scale for such a variable $\phi$ as:

```math
\begin{equation}
\label{eq:thermodynamic_scale}
\phi_{\star} = \frac{\kappa / Pr}{F_h(\Delta z,z_{0\phi}, L_O)} \Delta \phi,
\end{equation}
```
and a corresponding momentum scale is given by:
```math
\begin{equation}
\label{eq:u_star}
u_{\star} = \frac{\kappa}{F_m(\Delta z,z_{0m}, L_O)} \Delta u.
\end{equation}
```

Here $\kappa = 0.4$ is the von-Karman constant, $Pr$ is the Prandtl number and $\Delta$ corresponds to the difference in values between the two input heights (i.e. $\Delta u = u(z_{in}) - u(z_{srf})$). The formulation of $F_m$ and $F_h$ is subjective to the choice of discretization as shown in \citep{Nishizawa18a}. In the typical finite-difference form: 
```math
\begin{equation}
F_h^{(FD)} = \log \left(\frac{\Delta z}{z_{0b}} \right) - \psi_h \left(\frac{\Delta z}{L_O} \right) + \psi_h \left(\frac{z_{0b}}{L_O} \right),
\end{equation}
```
and 
```math
\begin{equation}
F_m^{(FD)} = \log \left(\frac{\Delta z}{z_{0m}} \right) - \psi_m \left(\frac{\Delta z}{L_O} \right) + \psi_m \left(\frac{z_{0m}}{L_O} \right).
\end{equation}
```
Subscripts $h, m$ are used to represent the equations corresponding to the exchange of heat and momentum respectively. The characteristics of tracer diffusion are assumed to be identical to those of the thermal diffusion in this system. The expressions for $\psi_m$ and $\psi_h$ are defined in appendix A in \citep{Nishizawa18a} for various universal functions. Our current approach uses the Businger universal functions by default. 

The correction for the finite volume form $F_h, F_m$ are given respectively as:
```math
\begin{equation}
F_h^{(FV)} = \log \left(\frac{\Delta z}{z_{0b}} \right) - \Psi_h \left(\frac{\Delta z}{L_O} \right) + \frac{z_{0b}}{\Delta z} \Psi_h \left(\frac{z_{0b}}{L_O} \right) + R_{z0h}\left[\Psi_h \left(\frac{z_{0b}}{L_O} \right) - 1\right]
\end{equation}
```
and 

```math
\begin{equation}
F_m^{(FV)} = \log \left(\frac{\Delta z}{z_{0m}} \right) - \Psi_m \left(\frac{\Delta z}{L_O} \right) + \frac{z_{0m}}{\Delta z} \Psi_m \left(\frac{z_{0m}}{L_O} \right) + R_{z0m}\left[\Psi_m \left(\frac{z_{0m}}{L_O} \right) - 1\right],
\end{equation}
```
with $\Psi_m$ and $\Psi_h$ as the FV corrected form. Their expressions are found in appendix A in \citep{Nishizawa18a} and here as well we use the Businger universal functions formulation as a default. 

Note that both $F_h$ and $F_m$ (for FV and FD) are functions of known parameters ($z_{0b}, z_{0m}, \Delta z$) and $L_O$ only. Furthermore we can compute $b_{\star}$ from the known thermodynamic states. This fact makes (1) a transcendental equation with $L_O$ as the only unknown. A numerical solver that iterate over $L_O$ can thus find the value for which the error in (1) is smaller than some desired tolerance. 

The most direct computation of $b_{\star}$ follows by writing (3) for virtual dry static energy and multiplying by a $(g/\mathrm{DSE}_v)$ factor:

```math
\begin{equation}
b_{\star} = \frac{g}{ \mathrm{DSE}_v }\frac{\kappa / Pr}{F_h(z,z_{0b}, L_O)} \Delta \mathrm{DSE}_v,
\end{equation}
```
with $\Delta \mathrm{DSE}_v$ computed from the thermodynamic states at the corresponding levels. This approach forces us to assume that the various scalars that take part in the virtual dry static energy has the same roughness length $z_{0b}$ but provides a straightforward implementation where any tracer (humidity or other) which plays a role in the virtual dry static energy of air is accounted for in the computation of $L_O$. 

An alternative method common in atmospheric models, is to write $b_{\star}$ using the definition of virtual potential temperature as:

```math
\begin{equation}
b_{\star} = (1+(\eps_{dv}-1) q_t ) \theta_{\star} +  (\eps_{dv}-1) \theta q_{t,\star},
\end{equation}
```
with $\theta_{\star}$ and $q_{t,\star}$ given by (3). Here separate roughness lengths can be assumed for heat and humidity but the implementation would have to change if any additional scalar is added to the air density. For these reasons we implement (10). 

The flux from the top of the boundary layer is given by:

```math
\begin{equation}
\overline{w'u'} = - u_{\star}^2,
\end{equation}
```

```math
\begin{equation}
\overline{w'\phi'} = - u_{\star} \phi_{\star},
\end{equation}
```
which could be viewed as alternatives to bulk formula expressions for the sensible and latent heat fluxes if the physical scales are known. Correspondingly, the ratios between the physical scale and its change in the surface layer provides and expression for the exchange coefficient:
```math
\begin{equation}
    C_{m} = \frac{u_{\star}^2}{\Delta u^2},
\end{equation}
```
```math
\begin{equation}
\label{eq:exchange_coeff_scalar}
    C_{\phi} = \frac{u_{\star} \phi_{\star}}{\Delta u \Delta\phi}.
\end{equation}
```
Equations (13) and (14) are singular if $\Delta u = 0$, and typically a gustiness parameter $w_{\star}$ is added to represent the SGS convective velocity and avoid this singularity. A typical formulation for gustiness is 
```math
\begin{equation}
    w_{\star} = \left(\overline{w'b'} h_{fc} \right)^{\frac{1}{3}},
\end{equation}
```
here $h_{fc}$ is the height of free convection, a non local estimation of the maximum height (above the surface layer) possible of convective motions based on the atmospheric profile. This formulation of gustiness requires a priori knowledge of the surface fluxes and thus can only be used in cases with prescribed fluxes. More generally gustiness is often taken to be a some assumed constant based on scaling arguments.

At the same time equation (14) is also singular if $\Delta\phi = 0$. If $\phi$ corresponds to buoyancy ($\mathrm{DSE}_v$) then the conditions are neutral, and the buoyancy flux is zero. $L_O$ is set to infinity, and the momentum exchange coefficient and heat exchange coefficient are given by the law of the wall:

```math
\begin{equation}
C_{m,N} = \left( \frac{\kappa}{ln [\Delta z/z_{0m}]} \right)^2.
\end{equation}
```

```math
\begin{equation}
C_{\phi,N} = \frac{\kappa^2}{ln [\Delta z/z_{0b}] ln [\Delta z/z_{0m}]}.
\end{equation}
```
In neutral conditions 
Once the surface conditions are known, the computation of the profile of any variable within the surface layer (for any $z_{sfc} < z < z_{in}$) is done by rewriting (3) and (4) as:

```math
\begin{equation}
\phi(z) = \phi_{\star}\frac{F_h(z,z_{0\phi}, L_O)}{\kappa / Pr}  + \phi_{sfc},
\end{equation}
```

```math
\begin{equation}
u(z) = u_{\star} \frac{F_m(z,z_{0m}, L_O)}{\kappa} + u_{sfc},
\end{equation}
```

with $\phi_{\star}$ and $u_{\star}$ given by (3) and (4) respectively. Here, the same choice of discretization should be used for $F_h$ and $F_m$ as used to obtain the surface conditions. 
