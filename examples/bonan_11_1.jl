# Run from the package root:
# julia --project=. examples/bonan_11_1.jl [output_directory]
# This host-side teaching example owns configuration, arrays, and file I/O.
using LeafPhotosynthesis
using TOML

"""Read all numerical settings as FT. The reusable physics never reads files."""
function read_bonan_parameters(::Type{FT}, path) where {FT}
    raw = TOML.parsefile(path)
    return Dict(section => Dict(key => value isa Real ? FT(value) : value
        for (key, value) in entries) for (section, entries) in raw)
end

"""Compose the leaf from named physical components and explicit parameters."""
function bonan_leaf(::Type{FT}, settings) where {FT}
    c, b = settings["constants"], settings["biochemistry"]
    et, ca = settings["electron_transport"], settings["carbon_assimilation"]
    co2_to_pa = c["pressure"] * c["micromole_per_mole"]
    o2_to_pa = c["pressure"] * c["millimole_per_mole"]
    arr(value, energy) = Arrhenius{FT}(; value_ref = value, activation_energy = energy,
        temperature_ref = c["temperature_ref"], gas_constant = c["gas_constant"])
    peak(value, name) = PeakedArrhenius(; arrhenius = arr(value, b[name * "_activation"]),
        deactivation_energy = b[name * "_deactivation"], entropy = b[name * "_entropy"])
    temperature = C3TemperatureResponse(
        vcmax = peak(b["vcmax_ref"], "vcmax"),
        jmax = peak(b["vcmax_ref"] * b["jmax_vcmax_ratio"], "jmax"),
        respiration = peak(b["vcmax_ref"] * b["respiration_vcmax_ratio"], "respiration"),
        kc = arr(b["kc_ref"] * co2_to_pa, b["kc_activation"]),
        ko = arr(b["ko_ref"] * o2_to_pa, b["ko_activation"]),
        gamma_star = arr(b["gamma_star_ref"] * co2_to_pa, b["gamma_star_activation"]),
    )
    transport = ElectronTransportParameters(; psii_fraction = et["psii_fraction"],
        psii_yield = et["psii_yield"], limitation = QuadraticLimitation(et["curvature"]))
    stoichiometry = C3Stoichiometry{FT}(;
        electrons_per_carboxylation = ca["electrons_per_carboxylation"],
        oxygenation_cost = ca["oxygenation_cost"])
    limitation = if ca["limitation"] == "quadratic"
        QuadraticLimitation(ca["curvature"])
    elseif ca["limitation"] == "minimum"
        MinimumLimitation()
    else
        error("unknown carbon limitation model")
    end
    compensation = if ca["compensation"] == "nonnegative"
        NonnegativeAssimilation()
    elseif ca["compensation"] == "signed"
        SignedAssimilation()
    else
        error("unknown compensation treatment")
    end
    return LeafModel(; photosynthesis = C3Photosynthesis(; temperature,
        electron_transport = transport, stoichiometry, limitation, compensation))
end

function bonan_environment(::Type{FT}, p; incident_par, internal_co2, temperature_celsius) where {FT}
    c, e = p["constants"], p["experiment"]
    absorptance = one(FT) - e["reflectance"] - e["transmittance"]
    return LeafEnvironment{FT}(;
        temperature = temperature_celsius + c["celsius_offset"],
        absorbed_par = incident_par * absorptance,
        co2 = internal_co2 * c["micromole_per_mole"] * c["pressure"],
        oxygen = e["oxygen"] * c["millimole_per_mole"] * c["pressure"])
end

"""Calculate all three response curves, calling the same public API used by a coupled model."""
function bonan_curves(::Type{FT}, settings) where {FT}
    leaf = bonan_leaf(FT, settings)
    e = settings["experiment"]
    curves = Dict{String, Vector{NamedTuple}}()
    for variable in ("light", "co2", "temperature")
        rows = NamedTuple[]
        for x in e[variable * "_start"]:e[variable * "_step"]:e[variable * "_stop"]
            incident_par = variable == "light" ? x : e["incident_par"]
            internal_co2 = variable == "co2" ? x :
                e["internal_co2_fraction"] * e["atmospheric_co2"]
            temperature_celsius = variable == "temperature" ? x : e["temperature_celsius"]
            env = bonan_environment(FT, settings; incident_par, internal_co2, temperature_celsius)
            result = photosynthesize(leaf, env)
            push!(rows, (; incident_par, temperature_celsius,
                atmospheric_co2 = internal_co2 / e["internal_co2_fraction"], internal_co2,
                rubisco_net = result.rubisco_limited - result.respiration,
                electron_net = result.electron_limited - result.respiration,
                net_assimilation = result.net_assimilation))
        end
        curves[variable] = rows
    end
    return curves
end

function write_bonan_curves(output_directory, curves)
    mkpath(output_directory)
    for (name, rows) in curves
        open(joinpath(output_directory, name * ".csv"), "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for row in rows
                println(io, join(values(row), ','))
            end
        end
    end
    return output_directory
end

if abspath(PROGRAM_FILE) == @__FILE__
    parameter_file = joinpath(@__DIR__, "parameters", "bonan_11_1.toml")
    output_directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "bonan_11_1") : only(ARGS)
    for FT in (Float32, Float64)
        settings = read_bonan_parameters(FT, parameter_file)
        curves = bonan_curves(FT, settings)
        write_bonan_curves(joinpath(output_directory, string(FT)), curves)
    end
    println("Response curves written to ", output_directory)
end
