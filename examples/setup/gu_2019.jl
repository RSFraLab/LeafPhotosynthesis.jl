using LeafPhotosynthesis: PSIIOpenFractionInputs, ObservedFluorescenceInputs, lake_model_yields,
    psii_maximum_photochemical_yield
isdefined(@__MODULE__, :read_bonan_parameters) || include(joinpath(@__DIR__, "read_parameters.jl"))
isdefined(@__MODULE__, :write_tables) || include(joinpath(@__DIR__, "write_tables.jl"))
isdefined(@__MODULE__, :setup_flexas_2002) || include(joinpath(@__DIR__, "flexas_2002.jl"))

function setup_gu_2019(::Type{FT}) where {FT<:AbstractFloat}
    settings = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "parameters", "gu_2019.toml"))
    key = FT === Float32 ? "float32" : "float64"
    yield_tol = settings["validation"][key * "_yield_tolerance"]
    lake_cfg = settings["lake_model"]
    # The parameters of the paper, for its own figures.
    paper_lake = LakeModelFluorescence(;
        maximum_photochemical_yield = lake_cfg["maximum_photochemical_yield"],
        heat_to_fluorescence_ratio = lake_cfg["heat_to_fluorescence_ratio"],
        yield_tolerance = yield_tol,
    )
    # The same PSII rate constants as the Johnson & Berry (2021) model of the Flexas
    # benchmark, so that only the fluorescence formulation differs in that comparison.
    flexas = setup_flexas_2002(FT)
    psii = flexas.settings["johnson_berry_psii"]
    host_lake = LakeModelFluorescence(;
        maximum_photochemical_yield = psii_maximum_photochemical_yield(flexas.fl_cyto),
        heat_to_fluorescence_ratio =
            (psii["constitutive_heat_rate"] + psii["sustained_quenching_rate"]) / psii["fluorescence_rate"],
        yield_tolerance = yield_tol,
    )
    # The lake relations applied to PAM quantities as measured: with Fo' from Oxborough &
    # Baker they hold identically between the detector levels Fo, Fm, Fs and Fm', whatever
    # PSI adds to them, when the maximum yield is the apparent Fv/Fm of the detector.
    apparent_lake = LakeModelFluorescence(;
        maximum_photochemical_yield = flexas.dark_pam.fv_fm,
        heat_to_fluorescence_ratio = lake_cfg["heat_to_fluorescence_ratio"],
        yield_tolerance = yield_tol,
    )
    return (; paper_lake, host_lake, apparent_lake, flexas, settings)
end

# Piecewise-linear q_L(PAR) through the digitized Fig. 1c points, from the dark value at
# PAR = 0 and held at the last point beyond it.
function measured_open_fraction(settings, par::FT) where {FT<:AbstractFloat}
    dig = settings["digitized_figure1c"]
    xs, ys = FT.(dig["par"]), FT.(dig["open_fraction"])
    par <= zero(FT) && return dig["dark_open_fraction"]
    par >= xs[end] && return ys[end]
    i = something(findlast(<=(par), xs), 0)
    x0, y0 = i == 0 ? (zero(FT), dig["dark_open_fraction"]) : (xs[i], ys[i])
    return FT(y0 + (ys[i + 1] - y0) * (par - x0) / (xs[i + 1] - x0))
end

fs_over_fo(model, result) =
    result.fluorescence_yield / lake_model_yields(model, one(result.fluorescence_yield), zero(result.fluorescence_yield)).fluorescence_yield

# Leaf, environment, stomatal state and diffusion with which the package's `gas_exchange`
# solves Gu et al. (2019), Eq. 24 at r_ch = 0 for a prescribed electron flow: the
# electron-limited rate of Eq. 4 with series diffusion through `boundary_stomatal_resistance`
# (r_sb) and the wall resistance r_wp [Pa m^2 s micromol^-1].
function prescribed_electron_flow_exchange(
    lake::LakeModelFluorescence{FT}, electron_flow::FT, boundary_stomatal_resistance::FT, settings,
) where {FT<:AbstractFloat}
    f2, host = settings["figure2"], settings["figure2_host"]
    photosynthesis = C3Photosynthesis(;
        temperature = C3TemperatureResponse(;
            vcmax = ConstantResponse(host["nonlimiting_vcmax"]),
            jmax = ConstantResponse(electron_flow),
            respiration = ConstantResponse(f2["respiration"]),
            kc = ConstantResponse(host["kc_Pa"]),
            ko = ConstantResponse(host["ko_Pa"]),
            gamma_star = ConstantResponse(f2["gamma_star_Pa"]),
        ),
        electron_transport = ElectronTransportParameters(;
            psii_fraction = host["psii_fraction"],
            psii_yield = lake.maximum_photochemical_yield,
            limitation = MinimumLimitation(),
        ),
        stoichiometry = C3Stoichiometry{FT}(;
            electrons_per_carboxylation = f2["electrons_per_carboxylation"],
            oxygenation_cost = f2["oxygenation_cost"],
        ),
        limitation = MinimumLimitation(),
        assimilation_sign = SignedAssimilation(),
    )
    # Molar conductance [mol m^-2 s^-1] of a resistance [Pa m^2 s micromol^-1].
    conductance(resistance) = host["pressure_Pa"] / (resistance * host["mole_to_micromole"])
    diffusion = SeriesDiffusion{FT}(;
        boundary_co2 = FT(Inf),
        mesophyll_co2 = conductance(f2["wall_resistance"]),
        diffusivity_ratio = host["diffusivity_ratio"],
        mole_to_micromole = host["mole_to_micromole"],
    )
    state = StomatalState(host["diffusivity_ratio"] * conductance(boundary_stomatal_resistance))
    environment = GasExchangeEnvironment(
        host["temperature_K"], host["nonlimiting_absorbed_par"], f2["ambient_co2_Pa"],
        host["oxygen_Pa"], host["pressure_Pa"], host["relative_humidity"], host["vpd_Pa"],
    )
    return (; photosynthesis, environment, state, diffusion)
end
