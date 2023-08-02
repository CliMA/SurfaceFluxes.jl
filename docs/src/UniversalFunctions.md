# Universal Functions

UniversalFunctions.jl provides universal functions for SurfaceFluxes.jl. The functions are defined in `src/UniversalFunctions.jl` and the plots are generated in `docs/src/plot_universal_functions.jl`. 

## Businger Universal Functions

Nishizawa (2018) Equation (A7)
```math
\begin{equation}
    f_m = (1 - 15\zeta)^{\frac{1}{4}}
\end{equation}
```

Nishizawa (2018) Equation (A8)
```math
\begin{equation}
    f_h = (1 - 9\zeta)^{\frac{1}{2}}
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
\phi_{m}(\zeta)=
    \begin{cases}
        a\zeta + 1 & \text{for } L\geq0\\
        (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A13)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        -\frac{15\zeta}{8} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A5)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        
         \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{1 - f_m^3}{12\zeta} + \frac{\pi}{2} - 1 & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A14)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
         -\frac{9\zeta}{4} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A6)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right) + \frac{2(1-f_h)}{9\zeta} - 1& \text{for } L < 0
        
    \end{cases}
\end{equation}
```

## Gryanik Universal Functions

Gryanik (2020) Equation (32)
```math
\begin{equation}
    \phi_{m}(\zeta) = 1 + \frac{a_m\zeta}{(1 + b_m\zeta)^\frac{2}{3}}
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
\phi_{m}(\zeta)=
    \begin{cases}
        a\zeta + 1 & \text{for } L\geq0\\
        (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Gryanik (2020) Equation (33)
```math
\begin{equation}
    \phi_{h}(\zeta) = Pr_{0}\left(1 + \frac{a_h\zeta}{1 + b_h\zeta}\right)
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Gryanik (2020) Equation (34)
```math
\begin{equation}
    \psi_{m}(\zeta) = -3\frac{a_m}{b_m}\left[(1+b_m)^\frac{1}{3} - 1 \right]
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Gryanik (2020) Equation (35)
```math
\begin{equation}
    \psi_{h}(\zeta) = 
    -Pr_{0}\frac{a_h}{b_h}\ln\left(1 + b_h\zeta\right)
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A13)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        -\frac{15\zeta}{8} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A5)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        
         \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{1 - f_m^3}{12\zeta} + \frac{\pi}{2} - 1 & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A14)
```math
\begin{equation}
    \Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
         -\frac{9\zeta}{4} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A6)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right) + \frac{2(1-f_h)}{9\zeta} - 1& \text{for } L < 0
        
    \end{cases}
\end{equation}
```

## Grachev Universal Functions

Grachev (2007) Equation (9a)
```math
\begin{equation}
\phi_{m}(\zeta) = 1 + \frac{a_m\zeta(1 + \zeta)^{\frac{1}{3}}}{1 + b_m\zeta}
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
\phi_{m}(\zeta)=
    \begin{cases}
        a\zeta + 1 & \text{for } L\geq0\\
        (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Grachev (2007) Equation (9b)
```math
\begin{equation}
\phi_{h}(\zeta) = 1 + \frac{a_h\zeta + b_h\zeta^2}{1 + c_h\zeta + \zeta^2}
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Grachev (2007) Equation 12
```math
\begin{equation}
\begin{split}
    \psi_{m}(\zeta) = -\frac{3a_m}{b_m}(\zeta - 1) + \frac{a_mb_m}{2b_m}\left[2\ln\left(\frac{\zeta + b_m}{1 + b_m}\right) - \ln\left(\frac{\zeta^2 - \zeta b_m + b_m^2}{1 - b_m + b_m^2}\right) + \\
    2\sqrt{3}\left(\arctan\left(\frac{2\zeta - b_m}{\sqrt{3}b_m}\right) - \arctan\left(\frac{2-b_m}{\sqrt{3}b_m} \right)\right) \right]
\end{split}
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Grachev (2007) Eq. 13
```math
\begin{equation}
    \psi_{h}(\zeta)= -\frac{b_h}{2}\ln\left(1 + c_h\zeta + \zeta^2\right) + \left(-\frac{a_h}{b_h} + \frac{b_hc_h}{2b_h}\right)\left[\ln\left(\frac{2\zeta + c_h - b_h}{2\zeta + c_h +b_h}\right) - \ln\left(\frac{c_h - b_h}{c_h + b_h}\right)\right] 
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

## Cheng Universal Functions

Cheng (2005) Equation 22
```math
\begin{equation}
    \phi_{\mathrm{m}}(\zeta)=1+a_m\left(\frac{\zeta+\zeta^{b_m}\left(1+\zeta^{b_m}\right)^{\frac{1-b_m}{b_m}}}{\zeta+\left(1+\zeta^{b_m}\right)^{\frac{1}{b_m}}}\right) .
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
    \phi_{m}(\zeta)=
        \begin{cases}
            a\zeta + 1 & \text{for } L\geq0\\
            (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
        \end{cases}
\end{equation}
```

Cheng (2005) Equation 24
```math
\begin{equation}
    \phi_{\mathrm{h}}(\zeta)=1+a_h\left(\frac{\zeta+\zeta^{b_h}\left(1+\zeta^{b_h}\right)^{\frac{1-b_h}{b_h}}}{\zeta+\left(1+\zeta^{b_h}\right)^{\frac{1}{b_h}}}\right) .
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Cheng (2005) Equation 21
```math
\begin{equation}
    \phi_{\mathrm{m}}(\zeta) = -a_m\ln\left[\zeta + \left(1 + \zeta^{b_m}\right)^{\frac{1}{b_m}}\right]
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Cheng (2005) Equation 23
```math
\begin{equation}
    \phi_{\mathrm{h}}(\zeta) = -a_h\ln\left[\zeta + \left(1 + \zeta^{b_h}\right)^{\frac{1}{b_h}}\right]
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

Volume-averaged form of Cheng (2005) Eq. 23 approximation (ζ  ≈ 0) with Nishizawa (2018) Eq. B7
```math
\begin{equation}
    \Psi_{h}(\zeta) \approx -\frac{a_h\zeta}{2}
\end{equation}
```

Nishizawa (2018) Equation (A14)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
         -\frac{9\zeta}{4} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A6)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right) + \frac{2(1-f_h)}{9\zeta} - 1& \text{for } L < 0
        
    \end{cases}
\end{equation}
```

Volume-averaged form of Cheng (2005) Eq. 21 approximation (ζ  ≈ 0) with Nishizawa (2018) Eq. B7
```math
\begin{equation}
    \Psi_{m}(\zeta) \approx -\frac{a_m\zeta}{2}
\end{equation}
```

Nishizawa (2018) Equation (A13)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        -\frac{15\zeta}{8} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A5)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        
         \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{1 - f_m^3}{12\zeta} + \frac{\pi}{2} - 1 & \text{for } L < 0
    \end{cases}
\end{equation}
```

## Holtslag Universal Functions

Derived from Holtslag (1988) Eq. 12 using Gyranik (2020) Eq. 11
```math
\begin{equation}
    \phi_{\mathrm{m}}(\zeta)= \zeta\left(a_m - b_md_me^{-d_m\zeta}\left(\zeta - \frac{c_m}{d_m}\right) + b_me^{-d_m\zeta}\right) + 1
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
\phi_{m}(\zeta)=
    \begin{cases}
        a\zeta + 1 & \text{for } L\geq0\\
        (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Holtslag (1988) Eq. 12 using Gyranik (2020) Eq. 11 
```math
\begin{equation}
    \phi_{\mathrm{h}}(\zeta)= \zeta\left(a_h - b_hd_he^{-d_h\zeta}\left(\zeta - \frac{c_h}{d_h}\right) + b_he^{-d_h\zeta}\right) + Pr_0
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Holtslag (1988) Eq. 12 
```math
\begin{equation}
    \psi_{m}(\zeta) = -a_m\zeta - b_m\left(\zeta - \frac{c_m}{d_m}\right)e^{-d_m\zeta} - \frac{b_mc_m}{d_m}
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Holtslag (1988) Eq. 12 (& 10)
```math
\begin{equation}
    \psi_{h}(\zeta) = -a_h\zeta - b_h\left(\zeta - \frac{c_h}{d_h}\right)e^{-d_h\zeta} - \frac{b_hc_h}{d_h}
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Holtslag (1988) Eq. 12 using Nishizawa (2018) Equation (14)
```math
\begin{equation}
    \Psi_{m}(\zeta) = b_m\left(c_m\left(\frac{1 - e^{-d_m\zeta}}{d_m^2\zeta} - \frac{1}{d_m}\right) + \frac{e^{-d_m\zeta} - 1}{d_m^2\zeta} + \frac{e^{-d_m\zeta}}{d_m}\right) - \frac{a_m\zeta}{2}
\end{equation}
```

Nishizawa (2018) Equation (A13)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        -\frac{15\zeta}{8} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A5)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        
         \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{1 - f_m^3}{12\zeta} + \frac{\pi}{2} - 1 & \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Holtslag (1988) Eq. 12 using Nishizawa (2018) Eq. 15
```math
\begin{equation}
    \Psi_{h}(\zeta) = b_h\left(c_h\left(\frac{1 - e^{-d_h\zeta}}{d_h^2\zeta} - \frac{1}{d_h}\right) + \frac{e^{-d_h\zeta} - 1}{d_h^2\zeta} + \frac{e^{-d_h\zeta}}{d_h}\right) - \frac{a_h\zeta}{2}
\end{equation}
```

Nishizawa (2018) Equation (A14)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
         -\frac{9\zeta}{4} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A6)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right) + \frac{2(1-f_h)}{9\zeta} - 1& \text{for } L < 0
        
    \end{cases}
\end{equation}
```

## Beljaars Universal Functions

Derived from Beljaars (1991) Eq. 28 using Gyranik (2020) Eq. 11
```math
\begin{equation}
    \phi_{\mathrm{m}}(\zeta)= \zeta\left(a_h - b_hd_he^{-d_h\zeta}\left(\zeta - \frac{c_h}{d_h}\right) + b_he^{-d_h\zeta}\right) + 1
\end{equation}
```

Nishizawa (2018) Equation (A1)
```math
\begin{equation}
\phi_{m}(\zeta)=
    \begin{cases}
        a\zeta + 1 & \text{for } L\geq0\\
        (1 - 15\zeta)^{-\frac{1}{4}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Beljaars (1991) Eq. 32 using Gyranik (2020) Eq. 11
```math
\begin{equation}
    \phi_{\mathrm{h}}(\zeta) = \zeta\left[b_he^{-d_h\zeta} + a_h\sqrt{\frac{2a_h\zeta}{3} + 1} - b_hd_he^{-d_h\zeta}\left(\zeta - \frac{c_h}{d_h}\right)\right] + Pr_0
\end{equation}
```

Nishizawa (2018) Equation (A2)
```math
\begin{equation}
\phi_{h}(\zeta)=
    \begin{cases}
        \frac{a\zeta}{\Pr} + 1 & \text{for } L\geq0\\
        (1 - 9\zeta)^{-\frac{1}{2}} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Beljaars (1991) Eq. 28
```math
\begin{equation}
    \psi_{m}(\zeta) = -a_m\zeta - b_m\left(\zeta - \frac{c_m}{d_m}\right)e^{-d_m\zeta} - \frac{b_mc_m}{d_m}
\end{equation}
```

Nishizawa (2018) Equation (A3)
```math
\begin{equation}
\psi_{m}(\zeta)=
    \begin{cases}
        -a\zeta & \text{for } L\geq0\\
        \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{\pi}{2} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Beljaars (1991) Eq. 32
```math
\begin{equation}
    \psi_{h}(\zeta) = -\left(1 + \frac{2a_h}{3}\zeta\right)^\frac{3}{2} + 1 - b_h\left(\zeta - \frac{c_h}{d_h}\right)e^{-d_h\zeta} - \frac{b_hc_h}{d_h}
\end{equation}
```

Nishizawa (2018) Equation (A4)
```math
\begin{equation}
\psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right)& \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Beljaars (1991) Eq. 28 using Nishizawa (2018) Equation (14)
```math
\begin{equation}
    \Psi_{m}(\zeta) = b_m\left(c_m\left(\frac{1 - e^{-d_m\zeta}}{d_m^2\zeta} - \frac{1}{d_m}\right) + \frac{e^{-d_m\zeta} - 1}{d_m^2\zeta} + \frac{e^{-d_m\zeta}}{d_m}\right) - \frac{a_m\zeta}{2}
\end{equation}
```

Nishizawa (2018) Equation (A13)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        -\frac{15\zeta}{8} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A5)
```math
\begin{equation}
\Psi_{m}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2} & \text{for } L\geq0\\
        
         \log\left[ \frac{(1 + f_{m})^2(1 + f_{m}^2)}{8} \right] - 2\tan^{-1}f_m + \frac{1 - f_m^3}{12\zeta} + \frac{\pi}{2} - 1 & \text{for } L < 0
    \end{cases}
\end{equation}
```

Derived from Beljaars (1991) Eq. 32 using Nishizawa (2018) Equation (15)
```math
\begin{equation}
    \Psi_{h}(\zeta) = -\frac{1}{\zeta}\left[
    \frac{\sqrt{3}\left(2a_h\zeta + 3\right)^\frac{5}{2} - 27}{45a_h} + \zeta\left(\frac{b_hc_h}{d_h} - 1\right) + b_h\left(\frac{1}{d_h^2} - \frac{e^{-d_h\zeta}\left(d_h\zeta + 1\right)}{d_h^2}\right) + \frac{b_hc_h\left(e^{-d_h\zeta} - 1\right)}{d_h^2}\right]
\end{equation}
```

Nishizawa (2018) Equation (A14)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
         -\frac{9\zeta}{4} & \text{for } L < 0
    \end{cases}
\end{equation}
```

Nishizawa (2018) Equation (A6)
```math
\begin{equation}
\Psi_{h}(\zeta)=
    \begin{cases}
        -\frac{a\zeta}{2\Pr} & \text{for } L\geq0\\
        
        2\log\left(\frac{1 + f_h}{2} \right) + \frac{2(1-f_h)}{9\zeta} - 1& \text{for } L < 0
        
    \end{cases}
\end{equation}
```

## Plots

Here, we reproduce some plots from literature, specifically from Gryanik et al. 2020, Businger 1971, and Bonan 2019. Note that Bonan uses the forms of $\phi$ and $\psi$ from Dyer and Hicks 1970; Dyer 1974; Brutsaert 1982, pp. 68–71; Garratt 1992, pp. 52–54. 

```@example
include("plot_universal_functions.jl")
```

### Figs 1,2 (Gryanik)

![](Gryanik12_phi_h.svg)
![](Gryanik12_phi_m.svg)
![](Gryanik12_psi_h.svg)
![](Gryanik12_psi_m.svg)

### Fig 3 (Gryanik)

![](Gryanik3_phi_h.svg)
![](Gryanik3_phi_m.svg)


### Figs 1,2 (Businger)

![](Businger_phi_h.svg)
![](Businger_phi_m.svg)


### Figs 1,2 (Bonan)

![](Bonan_phi_h.svg)
![](Bonan_phi_m.svg)
![](Bonan_psi_h.svg)
![](Bonan_psi_m.svg)
