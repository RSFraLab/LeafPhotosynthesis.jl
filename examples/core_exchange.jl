# Minimal leaf-only composition example; not a reproduction of Bonan 12.1 or Liu's experiments.
include("bonan_11_1.jl")

function core_leaf(::Type{FT}, biochemical_settings, settings) where {FT}
    base = bonan_leaf(FT, biochemical_settings).photosynthesis
    photosynthesis = C3Photosynthesis(; temperature = base.temperature,
        electron_transport = base.electron_transport, stoichiometry = base.stoichiometry,
        limitation = MinimumLimitation(), compensation = SignedAssimilation())
    s, d, o = settings["stomata"], settings["diffusion"], settings["oxygen"]
    target = BallBerry{FT}(; minimum_conductance = s["minimum_conductance"],
        slope = s["ball_berry_slope"], micromole_to_mole = s["micromole_to_mole"])
    stomata = PrognosticStomata(; target, opening_timescale = s["opening_timescale"],
        closing_timescale = s["closing_timescale"])
    oxygen = CoreOxygenExchange{FT}(; (Symbol(k) => v for (k, v) in o)...)
    diffusion = SeriesDiffusion{FT}(; (Symbol(k) => v for (k, v) in d)...)
    return LeafModel(; photosynthesis, stomata, oxygen), diffusion
end

function core_environment(::Type{FT}, p, settings) where {FT}
    c, e, x = p["constants"], p["experiment"], settings["experiment"]
    return GasExchangeEnvironment{FT}(;
        temperature = c["celsius_offset"] + e["temperature_celsius"],
        absorbed_par = e["incident_par"] * (one(FT) - e["reflectance"] - e["transmittance"]),
        co2 = e["atmospheric_co2"] * c["pressure"] * c["micromole_per_mole"],
        oxygen = e["oxygen"] * c["pressure"] * c["millimole_per_mole"],
        pressure = c["pressure"], relative_humidity = x["relative_humidity"], vpd = x["vpd"])
end

function core_solver(::Type{FT}, settings; method = SafeguardedNewton(;
    bisection_interval = Int(settings["solver"]["bisection_interval"]))) where {FT}
    key = FT === Float32 ? "solver_float32" : "solver_float64"
    s = settings[key]
    return EquilibriumSolver(; method, co2_atol = s["co2_atol"], co2_rtol = s["co2_rtol"],
        flux_atol = s["flux_atol"], flux_rtol = s["flux_rtol"],
        max_iterations = Int(settings["solver"]["max_iterations"]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "bonan_11_1.toml"))
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "parameters", "core_exchange.toml"))
        leaf, diffusion = core_leaf(FT, p, settings)
        environment = core_environment(FT, p, settings)
        state = StomatalState(settings["experiment"]["initial_conductance"])
        # This illustrative Euler integration belongs to the caller, not the library.
        dt, duration = settings["experiment"]["step"], settings["experiment"]["duration"]
        for time in dt:dt:duration
            result = leaf_tendencies(leaf, environment, state, diffusion)
            state = StomatalState(state.conductance + dt * result.conductance_tendency)
        end
        result = leaf_tendencies(leaf, environment, state, diffusion)
        println(FT, ": g_sw=", state.conductance, ", net CO2 uptake=",
            result.gas_exchange.carbon.net_assimilation, ", net O2 release=", result.oxygen.net_release)
        steady = LeafModel(; photosynthesis = leaf.photosynthesis,
            stomata = SteadyStateStomata(leaf.stomata.target))
        equilibrium = equilibrium_gas_exchange(steady, environment, diffusion, core_solver(FT, settings))
        equilibrium.converged || error("Equilibrium failed: $(equilibrium.status)")
        sensitivity = assimilation_sensitivity(steady.photosynthesis, environment,
            StomatalState(equilibrium.stomatal_conductance), diffusion, equilibrium)
        water = settings["water"]
        transpiration = transpiration_sensitivity(WaterVaporDiffusion(water["boundary_water"]),
            StomatalState(equilibrium.stomatal_conductance), water["leaf_to_air_vpd"], environment.pressure)
        println("  Equilibrium: gs=", equilibrium.stomatal_conductance, ", iterations=", equilibrium.iterations,
            "; dA/dgs=", sensitivity.dnet_assimilation_dgs, "; dE/dgs=", transpiration.dtranspiration_dgs)
    end
end
