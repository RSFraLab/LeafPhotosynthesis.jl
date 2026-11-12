isdefined(@__MODULE__, :read_bonan_parameters) || include(joinpath(@__DIR__, "read_parameters.jl"))
isdefined(@__MODULE__, :write_tables) || include(joinpath(@__DIR__, "write_tables.jl"))

function exact_dual_excitation(beta::FT, q_abs::FT) where {FT<:AbstractFloat}
    half = inv(FT(2))
    if beta >= half
        q2 = beta * q_abs
        return PrescribedDualExcitation(q2, q_abs - q2)
    else
        q1 = (one(FT) - beta) * q_abs
        return PrescribedDualExcitation(q_abs - q1, q1)
    end
end

function setup_johnson_berry_2021(::Type{FT}; capacity_scale::FT = one(FT)) where {FT}
    settings = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "parameters", "johnson_berry_2021.toml"))
    bio = settings["biochemistry"]
    stoich_cfg = settings["stoichiometry"]
    psii_cfg = settings["psii_rates"]
    psi_cfg = settings["psi_rates"]
    pam_cfg = settings["pam_observer"]
    exp_cfg = settings["experiment"]
    val_cfg = settings["validation"]
    key = FT === Float32 ? "float32" : "float64"
    yield_tol = val_cfg[key * "_yield_tolerance"]

    psi_params = PSIRateParameters(;
        fluorescence_rate = psi_cfg["fluorescence_rate"],
        constitutive_heat_rate = psi_cfg["constitutive_heat_rate"],
        maximum_photochemical_rate = psi_cfg["maximum_photochemical_rate"],
        closed_quenching_rate = psi_cfg["closed_quenching_rate"],
    )
    make_fl(ku2) = CytochromeFluorescence(;
        psii = PSIIRateParameters(;
            fluorescence_rate = psii_cfg["fluorescence_rate"],
            constitutive_heat_rate = psii_cfg["constitutive_heat_rate"],
            maximum_photochemical_rate = psii_cfg["maximum_photochemical_rate"],
            exciton_transfer_rate = ku2,
            sustained_quenching_rate = psii_cfg["sustained_quenching_rate"],
        ),
        psi = psi_params,
        yield_tolerance = yield_tol,
    )
    fl_puddle = make_fl(psii_cfg["puddle_transfer_rate"])
    fl_ku2 = make_fl(psii_cfg["moderate_transfer_rate"])
    fl_ku10 = make_fl(psii_cfg["high_transfer_rate"])

    # Rubisco specificity S_c/o = (kc / Kc) / (ko / Ko) and O2-scaled Gamma* = 0.5 * O2 / S_c/o.
    ko_turnover = bio["rubisco_kc_turnover"] * bio["rubisco_ko_kc_ratio"]
    specificity = (bio["rubisco_kc_turnover"] / bio["kc_Pa"]) * (bio["ko_Pa"] / ko_turnover)
    gamma_at_ref_o2 = bio["photorespiratory_co2_per_oxygenation"] * bio["reference_pressure_Pa"] / specificity

    vcmax = bio["vcmax_ref"] * capacity_scale
    cyto_temp = CytochromeTemperatureResponse(;
        vcmax = ConstantResponse(vcmax),
        vqmax = CytochromeCapacity(;
            site_density = bio["site_density"] * capacity_scale,
            turnover_rate = ConstantResponse(bio["turnover_ref"]),
        ),
        respiration = ConstantResponse(vcmax * bio["respiration_fraction"]),
        kc = ConstantResponse(bio["kc_Pa"]),
        ko = ConstantResponse(bio["ko_Pa"]),
        gamma_star = OxygenScaledCompensationPoint(
            ConstantResponse(gamma_at_ref_o2),
            bio["reference_pressure_Pa"],
        ),
        linear_atp_per_electron = ConstantResponse(bio["linear_atp_per_electron"]),
        cyclic_atp_per_electron = ConstantResponse(bio["cyclic_atp_per_electron"]),
    )
    stoich = CytochromeStoichiometry{FT}(;
        electrons_per_carboxylation = stoich_cfg["electrons_per_carboxylation"],
        electrons_per_photorespiratory_co2 = stoich_cfg["electrons_per_photorespiratory_co2"],
        atp_per_carboxylation = stoich_cfg["atp_per_carboxylation"],
        atp_per_photorespiratory_co2 = stoich_cfg["atp_per_photorespiratory_co2"],
    )
    et = CytochromeElectronParameters(;
        psii_fraction = bio["broadband_psii_fraction"],
        psi_maximum_yield = psi_maximum_photochemical_yield(fl_ku2),
        limitation = SerialLimitation(),
    )
    photo = CytochromePhotosynthesis(;
        temperature = cyto_temp,
        electron_transport = et,
        stoichiometry = stoich,
        limitation = MinimumLimitation(),
        compensation = SignedAssimilation(),
    )

    for fl in (fl_puddle, fl_ku2, fl_ku10)
        validate_fluorescence_coupling(fl, photo, exp_cfg["temperature_K"], exp_cfg["oxygen_Pa"])
    end

    rule = StateTransitionRule{FT}()
    pam_ideal = PAMFluorescenceObserver(;
        psii_sensitivity = pam_cfg["psii_sensitivity"],
        psi_sensitivity = pam_cfg["ideal_psi_sensitivity"],
    )
    pam_leaked = PAMFluorescenceObserver(;
        psii_sensitivity = pam_cfg["psii_sensitivity"],
        psi_sensitivity = pam_cfg["fluorometer_psi_sensitivity"],
    )
    beta_dark = dark_balanced_psii_partition(
        fl_ku2, photo, exp_cfg["temperature_K"], exp_cfg["co2_Pa"], exp_cfg["oxygen_Pa"])

    return (;
        photo,
        fl_puddle,
        fl_ku2,
        fl_ku10,
        rule,
        pam_ideal,
        pam_leaked,
        beta_dark,
        settings,
    )
end
