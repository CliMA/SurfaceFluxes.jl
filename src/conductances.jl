@inline heat_conductance(inputs::SurfaceFluxInputs, Ch, gustiness) =
    Ch * windspeed(inputs, gustiness)

