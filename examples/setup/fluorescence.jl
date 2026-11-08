isdefined(@__MODULE__, :core_leaf) || include(joinpath(@__DIR__, "..", "core_exchange.jl"))

function setup_fluorescence(::Type{FT}) where {FT}
    directory = joinpath(@__DIR__, "../parameters")
    settings = read_bonan_parameters(FT, joinpath(directory, "fluorescence.toml"))
    bio = read_bonan_parameters(FT, joinpath(directory, "bonan_11_1.toml"))
    core = read_bonan_parameters(FT, joinpath(directory, "core_exchange.toml"))
    base, _ = core_leaf(FT, bio, core)
    p, x = base.photosynthesis, settings["experiment"]
    transport = ElectronTransportParameters(;
        psii_fraction = x["psii_fraction"],
        psii_yield = x["psii_low_light_yield"],
        limitation = p.electron_transport.limitation,
    )
    photosynthesis = C3Photosynthesis(;
        p.temperature,
        electron_transport = transport,
        p.stoichiometry,
        p.limitation,
        p.compensation,
    )
    # No oxygen diagnostic is needed to obtain realized linear electron flow.
    leaf = LeafModel(; photosynthesis)
    rates, kn, validation =
        settings["rates"], settings["saturation_quenching"], settings["validation"]
    key = FT === Float32 ? "float32" : "float64"
    models = map((
        ConstantQuenching(rates["constant_quenching"]),
        SaturationQuenching(kn["maximum_rate"], kn["exponent"], kn["half_saturation"]),
    )) do quenching
        model = RateCompetitionFluorescence(;
            fluorescence_rate = rates["fluorescence"],
            constitutive_heat_rate = rates["constitutive_heat"],
            maximum_photochemical_rate = rates["maximum_photochemistry"],
            quenching,
            yield_tolerance = validation[key*"_yield_tolerance"],
        )
        validate_fluorescence_coupling(model, transport)
        model
    end
    environments = [
        LeafEnvironment(x["temperature_K"], light, x["biochemical_co2_Pa"], x["oxygen_Pa"]) for
        light = x["absorbed_par_start"]:x["absorbed_par_step"]:x["absorbed_par_stop"]
    ]
    return (; leaf, unquenched = models[1], regulated = models[2], environments, settings)
end

function fluorescence_row(environment, carbon, inputs, unquenched, regulated)
    fluorescence_valid(unquenched) && fluorescence_valid(regulated) ||
        error("Incompatible fluorescence rates: $(unquenched.status), $(regulated.status)")
    # Only valid dark yields are absent. Invalid sources never become plot gaps.
    yield_value(result, key) =
        result.status == FluorescenceDark ? nothing : getproperty(result, key)
    return (;
        absorbed_par = environment.absorbed_par,
        excitation = inputs.excitation,
        linear_electron_flow = inputs.linear_electron_flow,
        net_assimilation = carbon.net_assimilation,
        photochemical_yield = yield_value(regulated, :photochemical_yield),
        unquenched_yield = yield_value(unquenched, :fluorescence_yield),
        regulated_yield = yield_value(regulated, :fluorescence_yield),
        unquenched_photons = unquenched.photon_flux,
        regulated_photons = regulated.photon_flux,
        dark = regulated.status==FluorescenceDark,
        boundary_roundoff = unquenched.boundary_roundoff || regulated.boundary_roundoff,
    )
end

function write_fluorescence(rows, ::Type{FT}, directory) where {FT}
    write_tables(joinpath(directory, string(FT)), Dict("fluorescence"=>rows))
end
