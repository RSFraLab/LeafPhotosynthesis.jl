# Bonan 11.1: vary one forcing while keeping the others fixed.
# Run: julia --project=. examples/bonan_11_1.jl [output_directory]
using LeafPhotosynthesis
include("setup/read_parameters.jl")
include("setup/bonan_leaf.jl")
include("setup/write_tables.jl")

function bonan_curves(::Type{FT}, settings) where {FT}
    leaf = bonan_leaf(FT, settings)
    e = settings["experiment"]
    light = e["incident_par"]
    co2 = e["internal_co2_fraction"] * e["atmospheric_co2"]
    temperature = e["temperature_celsius"]

    # Evaluate the production model at a single forcing point.
    function response(incident_par, internal_co2, temperature_celsius)
        environment = bonan_environment(FT, settings;
            incident_par, internal_co2, temperature_celsius)
        result = photosynthesize(leaf, environment)
        return (; incident_par, temperature_celsius, internal_co2,
            atmospheric_co2 = internal_co2 / e["internal_co2_fraction"],
            rubisco_net = result.rubisco_limited - result.respiration,
            electron_net = result.electron_limited - result.respiration,
            net_assimilation = result.net_assimilation)
    end

    # Each curve changes just one input. No photosynthesis equations are duplicated.
    return Dict(
        "light" => [response(x, co2, temperature)
            for x in e["light_start"]:e["light_step"]:e["light_stop"]],
        "co2" => [response(light, x, temperature)
            for x in e["co2_start"]:e["co2_step"]:e["co2_stop"]],
        "temperature" => [response(light, co2, x)
            for x in e["temperature_start"]:e["temperature_step"]:e["temperature_stop"]],
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    output = isempty(ARGS) ? joinpath(@__DIR__, "output", "bonan_11_1") : only(ARGS)
    for FT in (Float32, Float64)
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "bonan_11_1.toml"))
        curves = bonan_curves(FT, settings)
        write_tables(joinpath(output, string(FT)), curves)
    end
    println("Response curves written to ", output)
end
