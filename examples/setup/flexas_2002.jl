isdefined(@__MODULE__, :read_bonan_parameters) || include(joinpath(@__DIR__, "read_parameters.jl"))
isdefined(@__MODULE__, :write_tables) || include(joinpath(@__DIR__, "write_tables.jl"))
isdefined(@__MODULE__, :exact_dual_excitation) || include(joinpath(@__DIR__, "johnson_berry_2021.jl"))

function nonstomatal_drought_factors(::Type{FT}, ns_cfg, g_sw_mmol::FT; enable::Bool = true) where {FT<:AbstractFloat}
    enable || return (one(FT), one(FT))
    g_half = ns_cfg["half_conductance_mmol"]
    n = ns_cfg["hill_exponent"]
    gn = g_sw_mmol^n
    hill = gn / (g_half^n + gn)
    vc_floor = ns_cfg["vcmax_severe_floor"]
    vq_floor = ns_cfg["vqmax_severe_floor"]
    f_vc = vc_floor + (one(FT) - vc_floor) * hill
    f_vq = vq_floor + (one(FT) - vq_floor) * hill
    return (FT(f_vc), FT(f_vq))
end

function setup_flexas_2002(::Type{FT}; g_sw_mmol::Union{Nothing,FT} = nothing, enable_nonstomatal::Bool = true) where {FT<:AbstractFloat}
    settings = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "parameters", "flexas_2002.toml"))
    bio = settings["grapevine_biochemistry"]
    stoich_cfg = settings["stoichiometry"]
    ns_cfg = settings["nonstomatal_drought"]
    diff_cfg = settings["diffusion"]
    psii_cfg = settings["johnson_berry_psii"]
    psi_cfg = settings["johnson_berry_psi"]
    pam_cfg = settings["pam_observer"]
    emp_rates = settings["empirical_rates"]
    emp_cot = settings["empirical_cotton"]
    emp_drt = settings["empirical_drought"]
    exp_cfg = settings["experiment"]
    val_cfg = settings["validation"]
    key = FT === Float32 ? "float32" : "float64"
    yield_tol = val_cfg[key * "_yield_tolerance"]
    g_sw = something(g_sw_mmol, FT(exp_cfg["fig2_conductances_mmol"][1]))

    f_vc, f_vq = nonstomatal_drought_factors(FT, ns_cfg, g_sw; enable = enable_nonstomatal)

    fl_cyto = CytochromeFluorescence(;
        psii = PSIIRateParameters(;
            fluorescence_rate = psii_cfg["fluorescence_rate"],
            constitutive_heat_rate = psii_cfg["constitutive_heat_rate"],
            maximum_photochemical_rate = psii_cfg["maximum_photochemical_rate"],
            exciton_transfer_rate = psii_cfg["exciton_transfer_rate"],
            sustained_quenching_rate = psii_cfg["sustained_quenching_rate"],
        ),
        psi = PSIRateParameters(;
            fluorescence_rate = psi_cfg["fluorescence_rate"],
            constitutive_heat_rate = psi_cfg["constitutive_heat_rate"],
            maximum_photochemical_rate = psi_cfg["maximum_photochemical_rate"],
            closed_quenching_rate = psi_cfg["closed_quenching_rate"],
        ),
        yield_tolerance = yield_tol,
    )

    gamma_comp = OxygenScaledCompensationPoint(
        ConstantResponse(bio["gamma_star_ref_Pa"]),
        bio["reference_oxygen_Pa"],
    )
    cyto_temp = CytochromeTemperatureResponse(;
        vcmax = ConstantResponse(bio["vcmax_ref"] * f_vc),
        vqmax = CytochromeCapacity(;
            site_density = bio["site_density"] * f_vq,
            turnover_rate = ConstantResponse(bio["turnover_ref"]),
        ),
        respiration = ConstantResponse(bio["respiration_ref"] * f_vc),
        kc = ConstantResponse(bio["kc_Pa"]),
        ko = ConstantResponse(bio["ko_Pa"]),
        gamma_star = gamma_comp,
        linear_atp_per_electron = ConstantResponse(bio["linear_atp_per_electron"]),
        cyclic_atp_per_electron = ConstantResponse(bio["cyclic_atp_per_electron"]),
    )
    cyto_stoich = CytochromeStoichiometry{FT}(;
        electrons_per_carboxylation = stoich_cfg["electrons_per_carboxylation"],
        electrons_per_photorespiratory_co2 = stoich_cfg["electrons_per_photorespiratory_co2"],
        atp_per_carboxylation = stoich_cfg["atp_per_carboxylation"],
        atp_per_photorespiratory_co2 = stoich_cfg["atp_per_photorespiratory_co2"],
    )
    cyto_et = CytochromeElectronParameters(;
        psii_fraction = bio["broadband_psii_fraction"],
        psi_maximum_yield = psi_maximum_photochemical_yield(fl_cyto),
        limitation = SerialLimitation(),
    )
    cyto_photo = CytochromePhotosynthesis(;
        temperature = cyto_temp,
        electron_transport = cyto_et,
        stoichiometry = cyto_stoich,
        limitation = MinimumLimitation(),
        compensation = SignedAssimilation(),
    )
    validate_fluorescence_coupling(
        fl_cyto, cyto_photo, exp_cfg["temperature_K"], exp_cfg["ambient_oxygen_Pa"],
    )

    fvcb_temp = C3TemperatureResponse(;
        vcmax = ConstantResponse(bio["vcmax_ref"] * f_vc),
        jmax = ConstantResponse(bio["jmax_ref"] * f_vq),
        respiration = ConstantResponse(bio["respiration_ref"] * f_vc),
        kc = ConstantResponse(bio["kc_Pa"]),
        ko = ConstantResponse(bio["ko_Pa"]),
        gamma_star = gamma_comp,
    )
    fvcb_et = ElectronTransportParameters(;
        psii_fraction = emp_rates["psii_fraction"],
        psii_yield = emp_rates["psii_low_light_yield"],
        limitation = QuadraticLimitation(bio["fvcb_curvature"]),
    )
    fvcb_stoich = C3Stoichiometry{FT}(;
        electrons_per_carboxylation = stoich_cfg["electrons_per_carboxylation"],
        oxygenation_cost = stoich_cfg["fvcb_oxygenation_cost"],
    )
    fvcb_photo = C3Photosynthesis(;
        temperature = fvcb_temp,
        electron_transport = fvcb_et,
        stoichiometry = fvcb_stoich,
        limitation = MinimumLimitation(),
        compensation = SignedAssimilation(),
    )

    unquenched_fl = RateCompetitionFluorescence(;
        fluorescence_rate = emp_rates["fluorescence_rate"],
        constitutive_heat_rate = emp_rates["constitutive_heat_rate"],
        maximum_photochemical_rate = emp_rates["maximum_photochemical_rate"],
        quenching = ConstantQuenching(zero(FT)),
        yield_tolerance = yield_tol,
    )
    cotton_fl, drought_fl = map((emp_cot, emp_drt)) do cfg
        model = RateCompetitionFluorescence(;
            fluorescence_rate = emp_rates["fluorescence_rate"],
            constitutive_heat_rate = emp_rates["constitutive_heat_rate"],
            maximum_photochemical_rate = emp_rates["maximum_photochemical_rate"],
            quenching = SaturationQuenching(
                cfg["maximum_rate"],
                cfg["exponent"],
                cfg["half_saturation"],
            ),
            yield_tolerance = yield_tol,
        )
        validate_fluorescence_coupling(model, fvcb_et)
        model
    end
    validate_fluorescence_coupling(unquenched_fl, fvcb_et)

    empirical_fo_yield = emp_rates["fluorescence_rate"] / (
        emp_rates["fluorescence_rate"] +
        emp_rates["constitutive_heat_rate"] +
        emp_rates["maximum_photochemical_rate"]
    )

    diffusion = SeriesDiffusion{FT}(;
        boundary_co2 = diff_cfg["boundary_co2"],
        mesophyll_co2 = diff_cfg["mesophyll_co2"],
        diffusivity_ratio = diff_cfg["diffusivity_ratio"],
        mole_to_micromole = diff_cfg["mole_to_micromole"],
    )
    pam_ideal = PAMFluorescenceObserver(;
        psii_sensitivity = pam_cfg["psii_sensitivity"],
        psi_sensitivity = pam_cfg["ideal_psi_sensitivity"],
    )
    pam_pam2000 = PAMFluorescenceObserver(;
        psii_sensitivity = pam_cfg["psii_sensitivity"],
        psi_sensitivity = pam_cfg["pam2000_psi_sensitivity"],
    )
    rule = StateTransitionRule{FT}()
    # Dark-balanced PSII excitation share at ambient air; under connected PSII units
    # (KU2 > 0), Phi_P2,max(q2) > q2 * Phi_P2,max(1) for all lit q2 < 1, keeping KN >= 0
    # across both ambient air and the lit PPFD >= 25 high-CO2 / low-O2 reversion grid.
    beta_dark = dark_balanced_psii_partition(
        fl_cyto,
        cyto_photo,
        exp_cfg["temperature_K"],
        exp_cfg["ambient_co2_Pa"],
        exp_cfg["ambient_oxygen_Pa"],
    )
    absorptance = exp_cfg["leaf_absorptance"]
    alpha2_dark = absorptance * beta_dark
    alpha1_dark = absorptance * (one(FT) - beta_dark)
    dark_pam = dark_pam_levels(pam_pam2000, fl_cyto, alpha2_dark, alpha1_dark)
    dark_psii = dark_pam_levels(pam_ideal, fl_cyto, alpha2_dark, alpha1_dark)
    empirical_fv_fm = FT(emp_rates["psii_low_light_yield"])

    return (;
        cyto_photo,
        fvcb_photo,
        fl_cyto,
        unquenched_fl,
        cotton_fl,
        drought_fl,
        empirical_fo_yield,
        empirical_fv_fm,
        f_vc,
        f_vq,
        g_sw,
        diffusion,
        pam_ideal,
        pam_pam2000,
        rule,
        beta_dark,
        alpha2_dark,
        alpha1_dark,
        dark_pam,
        dark_psii,
        stomatal_state = StomatalState(g_sw * exp_cfg["mmol_to_mol"]),
        settings,
    )
end

function evaluate_flexas_point(
    case,
    incident_par::FT;
    co2_Pa::FT = case.settings["experiment"]["ambient_co2_Pa"],
    oxygen_Pa::FT = case.settings["experiment"]["ambient_oxygen_Pa"],
) where {FT<:AbstractFloat}
    x = case.settings["experiment"]
    T = x["temperature_K"]
    absorptance = x["leaf_absorptance"]
    q_abs = absorptance * incident_par
    genv = GasExchangeEnvironment(
        T, q_abs, co2_Pa, oxygen_Pa, x["pressure_Pa"], x["relative_humidity"], x["vpd_Pa"],
    )
    if iszero(q_abs)
        ge_dark = gas_exchange(
            case.cyto_photo,
            genv,
            case.stomatal_state,
            case.diffusion,
            exact_dual_excitation(case.beta_dark, zero(FT)),
        )
        isfinite(ge_dark.biochemical_co2) || error("Nonfinite dark Cytochrome Cc at PPFD = $incident_par")
        return (;
            incident_par,
            absorbed_par = q_abs,
            pam2000_fs_over_fo = one(FT),
            psii_fs_over_fo = one(FT),
            qt_pam2000_fs_over_fo = one(FT),
            cotton_fs_over_fo = nothing,
            drought_fs_over_fo = nothing,
            unquenched_fs_over_fo = nothing,
            net_assimilation = ge_dark.carbon.net_assimilation,
            biochemical_co2 = ge_dark.biochemical_co2,
            psii_open_fraction = one(FT),
            psii_photochemical_yield = nothing,
            pam_apparent_photochemical_yield = nothing,
            psii_quenching_rate = zero(FT),
            pam2000_npq = zero(FT),
            dark = true,
        )
    end

    static_exc = exact_dual_excitation(case.beta_dark, q_abs)
    ge_cyto = gas_exchange(case.cyto_photo, genv, case.stomatal_state, case.diffusion, static_exc)
    isfinite(ge_cyto.biochemical_co2) || error("Nonfinite Cytochrome Cc at PPFD = $incident_par")
    site_cyto = LeafEnvironment(T, q_abs, ge_cyto.biochemical_co2, oxygen_Pa)
    r_cyto = fluorescence(case.fl_cyto, photosystem_flux_inputs(case.cyto_photo, site_cyto, static_exc))
    fluorescence_valid(r_cyto) || error("Invalid Cytochrome fluorescence at PPFD = $incident_par")
    p_pam = pam_fluorescence(case.pam_pam2000, case.fl_cyto, r_cyto, case.alpha2_dark, case.alpha1_dark)
    p_psii = pam_fluorescence(case.pam_ideal, case.fl_cyto, r_cyto, case.alpha2_dark, case.alpha1_dark)
    (pam_fluorescence_valid(p_pam) && pam_fluorescence_valid(p_psii)) ||
        error("Invalid PAM observation at PPFD = $incident_par")

    # Coupled state-transition (qT) partition and prescribed-conductance CO2 diffusion.
    cc_qt = ge_cyto.biochemical_co2
    exc_qt = static_exc
    beta_qt = case.beta_dark
    for _ in 1:4
        site_iter = LeafEnvironment(T, q_abs, cc_qt, oxygen_Pa)
        exc_qt = state_transition_excitation(case.rule, case.fl_cyto, case.cyto_photo, site_iter)
        beta_qt = state_transition_partition(case.rule, case.fl_cyto, case.cyto_photo, site_iter)
        ge_iter = gas_exchange(case.cyto_photo, genv, case.stomatal_state, case.diffusion, exc_qt)
        isfinite(ge_iter.biochemical_co2) || error("Nonfinite qT Cc at PPFD = $incident_par")
        cc_qt = ge_iter.biochemical_co2
    end
    site_qt = LeafEnvironment(T, q_abs, cc_qt, oxygen_Pa)
    r_qt = fluorescence(case.fl_cyto, photosystem_flux_inputs(case.cyto_photo, site_qt, exc_qt))
    fluorescence_valid(r_qt) || error("Invalid qT fluorescence at PPFD = $incident_par")
    p_qt = pam_fluorescence(
        case.pam_pam2000,
        case.fl_cyto,
        r_qt,
        absorptance * beta_qt,
        absorptance * (one(FT) - beta_qt),
        case.alpha2_dark,
        case.alpha1_dark,
    )
    pam_fluorescence_valid(p_qt) || error("Invalid qT PAM observation at PPFD = $incident_par")

    ge_fvcb = gas_exchange(case.fvcb_photo, genv, case.stomatal_state, case.diffusion)
    isfinite(ge_fvcb.biochemical_co2) || error("Nonfinite FvCB Cc at PPFD = $incident_par")
    site_fvcb = LeafEnvironment(T, q_abs, ge_fvcb.biochemical_co2, oxygen_Pa)
    fvcb_inp = psii_flux_inputs(case.fvcb_photo, site_fvcb)
    r_cot = fluorescence(case.cotton_fl, fvcb_inp)
    r_drt = fluorescence(case.drought_fl, fvcb_inp)
    r_unq = fluorescence(case.unquenched_fl, fvcb_inp)
    all(fluorescence_valid, (r_cot, r_drt, r_unq)) ||
        error("Invalid empirical fluorescence at PPFD = $incident_par")

    return (;
        incident_par,
        absorbed_par = q_abs,
        pam2000_fs_over_fo = FT(p_pam.fs / p_pam.fo_dark),
        psii_fs_over_fo = FT(p_psii.fs / p_psii.fo_dark),
        qt_pam2000_fs_over_fo = FT(p_qt.fs / p_qt.fo_dark),
        cotton_fs_over_fo = FT(r_cot.fluorescence_yield / case.empirical_fo_yield),
        drought_fs_over_fo = FT(r_drt.fluorescence_yield / case.empirical_fo_yield),
        unquenched_fs_over_fo = FT(r_unq.fluorescence_yield / case.empirical_fo_yield),
        net_assimilation = ge_cyto.carbon.net_assimilation,
        biochemical_co2 = ge_cyto.biochemical_co2,
        psii_open_fraction = r_cyto.psii_open_fraction,
        psii_photochemical_yield = r_cyto.psii_photochemical_yield,
        pam_apparent_photochemical_yield = p_pam.apparent_photochemical_yield,
        psii_quenching_rate = r_cyto.psii_quenching_rate,
        pam2000_npq = p_pam.apparent_npq,
        dark = false,
    )
end

function flexas_rmse(::Type{FT}, predicted, observed) where {FT<:AbstractFloat}
    n = length(observed)
    sq_sum = sum((FT(p) - FT(o))^2 for (p, o) in zip(predicted, observed))
    return FT(sqrt(sq_sum / FT(n)))
end

function flexas_lit_peak(::Type{FT}, rows, field::Symbol) where {FT<:AbstractFloat}
    lit_rows = [r for r in rows if !r.dark]
    vals = [FT(getproperty(r, field)) for r in lit_rows]
    peak_val, idx = findmax(vals)
    return (; peak = peak_val, ppfd = lit_rows[idx].incident_par)
end
