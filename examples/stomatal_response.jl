# Stomata open in light and close in darkness, with different time constants.
# Run: julia --project=. examples/stomatal_response.jl [output_directory]
include("core_exchange.jl")

function stomatal_response(::Type{FT}, biochemistry, settings) where {FT}
    leaf, diffusion = core_leaf(FT, biochemistry, settings)
    light = core_environment(FT, biochemistry, settings)
    experiment = settings["experiment"]
    validate_time_steps(experiment)
    dark = GasExchangeEnvironment(light.temperature, experiment["dark_absorbed_par"],
        light.co2, light.oxygen, light.pressure, light.relative_humidity, light.vpd)
    steady = LeafModel(; photosynthesis = leaf.photosynthesis,
        stomata = SteadyStateStomata(leaf.stomata.target))
    solver = core_solver(FT, settings)
    state = StomatalState(experiment["initial_conductance"])
    dt = experiment["step"]
    rows = NamedTuple[]
    elapsed = zero(FT)

    for (environment, duration) in ((light, experiment["duration"]), (dark, experiment["dark_duration"]))
        equilibrium = equilibrium_gas_exchange(steady, environment, diffusion, solver)
        equilibrium.converged || error("Equilibrium failed: $(equilibrium.status)")
        # Save before each update, including the endpoints on both sides of the light switch.
        steps = round(Int, duration / dt)
        for step in 0:steps
            result = leaf_tendencies(leaf, environment, state, diffusion)
            push!(rows, (; time = elapsed + step * dt, absorbed_par = environment.absorbed_par,
                conductance = state.conductance, equilibrium_conductance = equilibrium.stomatal_conductance,
                net_assimilation = result.gas_exchange.carbon.net_assimilation,
                equilibrium_assimilation = equilibrium.carbon.net_assimilation,
                net_oxygen = result.oxygen.net_release))
            # Illustrative Euler step: the caller owns time integration.
            if step < steps
                state = StomatalState(state.conductance + dt * result.conductance_tendency)
            end
        end
        elapsed += duration
    end
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    output = isempty(ARGS) ? joinpath(@__DIR__, "output", "stomatal_response") : only(ARGS)
    for FT in (Float32, Float64)
        biochemistry = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "bonan_11_1.toml"))
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "core_exchange.toml"))
        rows = stomatal_response(FT, biochemistry, settings)
        write_tables(joinpath(output, string(FT)), Dict("stomata" => rows))
    end
    println("Stomatal response written to ", output)
end
