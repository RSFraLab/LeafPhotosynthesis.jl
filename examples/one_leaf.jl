# Start here: evaluate one leaf with prescribed internal CO₂ and temperature.
# Run: julia --project=. examples/one_leaf.jl
using LeafPhotosynthesis
include("setup/read_parameters.jl")
include("setup/bonan_leaf.jl")

FT = Float64                         # Change to Float32 to use single precision.
settings = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "bonan_11_1.toml"))
leaf = bonan_leaf(FT, settings)
forcing = settings["experiment"]

environment = bonan_environment(FT, settings;
    incident_par = forcing["incident_par"],
    internal_co2 = forcing["internal_co2_fraction"] * forcing["atmospheric_co2"],
    temperature_celsius = forcing["temperature_celsius"])

result = photosynthesize(leaf, environment)
println("Net CO₂ uptake: ", result.net_assimilation, " μmol CO₂ m⁻² leaf s⁻¹")
