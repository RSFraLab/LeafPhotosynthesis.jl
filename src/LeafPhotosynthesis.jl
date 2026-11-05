module LeafPhotosynthesis

export AbstractPhotosynthesis, C3Photosynthesis, LeafModel, LeafEnvironment,
    Arrhenius, PeakedArrhenius, temperature_response, C3TemperatureResponse,
    C3Rates, biochemical_rates, ElectronTransportParameters, electron_transport,
    C3Stoichiometry, MinimumLimitation, QuadraticLimitation, colimit,
    SignedAssimilation, NonnegativeAssimilation, carbon_assimilation, photosynthesize

# Follow this list from physical responses to a complete prescribed-CO₂ leaf call.
include("temperature_response.jl")
include("rate_limitation.jl")
include("electron_transport.jl")
include("carbon_assimilation.jl")
include("photosynthesis.jl")
include("leaf_model.jl")

end
