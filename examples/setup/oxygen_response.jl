isdefined(@__MODULE__, :core_leaf) || include(joinpath(@__DIR__, "..", "core_exchange.jl"))

"""Replace only the Γ* component while retaining every other leaf component."""
function with_compensation_point(leaf, component)
    p = leaf.photosynthesis
    t = p.temperature
    temperature = C3TemperatureResponse(;
        t.vcmax,
        t.jmax,
        t.respiration,
        t.kc,
        t.ko,
        gamma_star = component,
    )
    photosynthesis = C3Photosynthesis(;
        temperature,
        p.electron_transport,
        p.stoichiometry,
        p.limitation,
        p.compensation,
    )
    return LeafModel(; photosynthesis, leaf.stomata, leaf.oxygen)
end

function oxygen_response_environments(::Type{FT}, case) where {FT}
    x, e = case.settings["experiment"], case.environment
    return [
        GasExchangeEnvironment(
            T,
            FT(x["response_absorbed_par_micromol_m2_s"]),
            FT(x["co2_Pa"]),
            O,
            e.pressure,
            e.relative_humidity,
            e.vpd,
        ) for T in FT.(x["temperature_K"]) for O in FT.(x["oxygen_Pa"])
    ]
end

function write_oxygen_response(rows, ::Type{FT}, directory) where {FT}
    write_tables(joinpath(directory, string(FT)), Dict("oxygen_response" => rows))
end

function setup_oxygen_response(::Type{FT}) where {FT}
    directory = joinpath(@__DIR__, "../parameters")
    biochemistry = read_bonan_parameters(FT, joinpath(directory, "bonan_11_1.toml"))
    core = read_bonan_parameters(FT, joinpath(directory, "core_exchange.toml"))
    settings = TOML.parsefile(joinpath(directory, "oxygen_response.toml"))
    leaf, diffusion = core_leaf(FT, biochemistry, core)
    response = leaf.photosynthesis.temperature.gamma_star
    oxygen_reference = FT(settings["experiment"]["reference_oxygen_Pa"])
    prescribed = with_compensation_point(leaf, PrescribedCompensationPoint(response))
    scaled = with_compensation_point(
        leaf,
        OxygenScaledCompensationPoint(response, oxygen_reference),
    )
    environment = core_environment(FT, biochemistry, core)
    return (;
        legacy = leaf,
        prescribed,
        scaled,
        diffusion,
        environment,
        solver = core_solver(FT, core),
        settings,
        oxygen_reference,
    )
end
