import cdsapi

filename = "test_global.nc"
filename_sfc = "test_global_sfc.nc"

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'pressure_level': '950',
        'variable': [
            'geopotential', 'specific_humidity', 'temperature',
            'u_component_of_wind', 'v_component_of_wind',
        ],
        'year': '2022',
        'month': '11',
        'day': '01',
        'time': [
            '03:00', '15:00',
        ],
    },
    '%s'%filename )


c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'forecast_logarithm_of_surface_roughness_for_heat', 'geopotential', '2m_dewpoint_temperature', '2m_temperature', 'eastward_turbulent_surface_stress',
            'evaporation', 'forecast_surface_roughness', 'friction_velocity',
            'instantaneous_moisture_flux', 'northward_turbulent_surface_stress', 'potential_evaporation',
            'surface_latent_heat_flux', 'surface_pressure', 'surface_sensible_heat_flux',
        ],
        'time': [
            '03:00', '15:00',
        ],
        'year': '2022',
        'month': '11',
        'day': '01',
    },
    '%s'%filename_sfc )