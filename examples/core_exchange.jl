# A leaf at equilibrium: assimilation, stomata, oxygen, and sensitivities.
# Illustrative parameters; not a reproduction of Bonan 12.1.
# Run: julia --project=. examples/core_exchange.jl
include("bonan_11_1.jl")
include("setup/core_leaf.jl")

function equilibrium_example(::Type{FT}, biochemistry, settings) where {FT}
    leaf, diffusion = core_leaf(FT, biochemistry, settings)
    environment = core_environment(FT, biochemistry, settings)
    steady_leaf = LeafModel(; photosynthesis = leaf.photosynthesis,
        stomata = SteadyStateStomata(leaf.stomata.target), oxygen = leaf.oxygen)

    # Solve the coupled assimilation–diffusion–stomatal equilibrium.
    exchange = equilibrium_gas_exchange(steady_leaf, environment, diffusion, core_solver(FT, settings))
    exchange.converged || error("Equilibrium failed: $(exchange.status)")
    state = StomatalState(exchange.stomatal_conductance)

    # Evaluate carbon and water responses to opening the stomata, at fixed T.
    carbon = assimilation_sensitivity(steady_leaf.photosynthesis, environment, state, diffusion, exchange)
    water_settings = settings["water"]
    water = transpiration_sensitivity(WaterVaporDiffusion(water_settings["boundary_water"]),
        state, water_settings["leaf_to_air_vpd"], environment.pressure)

    # Oxygen uses the solved CO₂ concentration inside the leaf.
    site = LeafEnvironment(environment.temperature, environment.absorbed_par,
        exchange.biochemical_co2, environment.oxygen)
    oxygen = oxygen_exchange(steady_leaf, site)
    return (; exchange, carbon, water, oxygen)
end

if abspath(PROGRAM_FILE) == @__FILE__
    for FT in (Float32, Float64)
        biochemistry = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "bonan_11_1.toml"))
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "core_exchange.toml"))
        result = equilibrium_example(FT, biochemistry, settings)
        println(FT, ": ", result)
    end
end
