# Explicit component assembly for the illustrative gas-exchange examples.
function validate_time_steps(experiment)
    dt = experiment["step"]
    isfinite(dt) && dt > zero(dt) || error("The example time step must be finite and positive")
    for key in ("duration", "dark_duration")
        duration = experiment[key]
        isfinite(duration) && duration > zero(duration) || error("$key must be finite and positive")
        isinteger(duration / dt) || error("The time step must divide $key exactly")
    end
end

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
