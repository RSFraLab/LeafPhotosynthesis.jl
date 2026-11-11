isdefined(@__MODULE__, :core_leaf) || include(joinpath(@__DIR__, "..", "core_exchange.jl"))

function scale_arrhenius_like(law::Arrhenius{FT}, factor::FT) where {FT}
    return Arrhenius{FT}(;
        value_ref = law.value_ref * factor,
        law.activation_energy,
        law.temperature_ref,
        law.gas_constant,
    )
end

function scale_arrhenius_like(law::PeakedArrhenius{FT}, factor::FT) where {FT}
    return PeakedArrhenius(;
        arrhenius = scale_arrhenius_like(law.arrhenius, factor),
        law.deactivation_energy,
        law.entropy,
    )
end

function setup_cytochrome_stress(::Type{FT}; vcmax_factor::FT = one(FT),
    transport_factor::FT = one(FT)) where {FT}
    directory = joinpath(@__DIR__, "../parameters")
    settings = read_bonan_parameters(FT, joinpath(directory, "cytochrome_stress.toml"))
    bio = read_bonan_parameters(FT, joinpath(directory, "bonan_11_1.toml"))
    core = read_bonan_parameters(FT, joinpath(directory, "core_exchange.toml"))
    base_leaf, diffusion = core_leaf(FT, bio, core)
    base_temp = base_leaf.photosynthesis.temperature

    c = bio["constants"]
    b6f = settings["cytochrome_b6f"]
    stoich_cfg = settings["cytochrome_stoichiometry"]
    psii_cfg = settings["johnson_berry_psii"]
    psi_cfg = settings["johnson_berry_psi"]
    emp_rates = settings["empirical_rates"]
    emp_cot = settings["empirical_cotton"]
    emp_drt = settings["empirical_drought"]
    val = settings["validation"]
    key = FT === Float32 ? "float32" : "float64"
    yield_tol = val[key * "_yield_tolerance"]

    # 1. Johnson & Berry (2021) coupled PSII + PSI fluorescence models (connected vs. puddle).
    psi_params = PSIRateParameters(;
        fluorescence_rate = psi_cfg["fluorescence_rate"],
        constitutive_heat_rate = psi_cfg["constitutive_heat_rate"],
        maximum_photochemical_rate = psi_cfg["maximum_photochemical_rate"],
        closed_quenching_rate = psi_cfg["closed_quenching_rate"],
    )
    jb_connected = CytochromeFluorescence(;
        psii = PSIIRateParameters(;
            fluorescence_rate = psii_cfg["fluorescence_rate"],
            constitutive_heat_rate = psii_cfg["constitutive_heat_rate"],
            maximum_photochemical_rate = psii_cfg["maximum_photochemical_rate"],
            exciton_transfer_rate = psii_cfg["exciton_transfer_rate"],
            sustained_quenching_rate = psii_cfg["sustained_quenching_rate"],
        ),
        psi = psi_params,
        yield_tolerance = yield_tol,
    )
    jb_puddle = CytochromeFluorescence(;
        psii = PSIIRateParameters(;
            fluorescence_rate = psii_cfg["fluorescence_rate"],
            constitutive_heat_rate = psii_cfg["constitutive_heat_rate"],
            maximum_photochemical_rate = psii_cfg["maximum_photochemical_rate"],
            exciton_transfer_rate = psii_cfg["puddle_transfer_rate"],
            sustained_quenching_rate = psii_cfg["sustained_quenching_rate"],
        ),
        psi = psi_params,
        yield_tolerance = yield_tol,
    )

    # 2. Cytochrome b6f C3 photosynthesis sharing exact Rubisco kinetics with FvCB.
    kq_law = Arrhenius{FT}(;
        value_ref = b6f["turnover_ref"],
        activation_energy = b6f["turnover_activation"],
        temperature_ref = c["temperature_ref"],
        gas_constant = c["gas_constant"],
    )
    cyto_temp = CytochromeTemperatureResponse(;
        vcmax = scale_arrhenius_like(base_temp.vcmax, vcmax_factor),
        vqmax = CytochromeCapacity(;
            site_density = b6f["site_density"] * transport_factor,
            turnover_rate = kq_law,
        ),
        respiration = base_temp.respiration,
        kc = base_temp.kc,
        ko = base_temp.ko,
        gamma_star = base_temp.gamma_star,
        linear_atp_per_electron = ConstantResponse(b6f["linear_atp_per_electron"]),
        cyclic_atp_per_electron = ConstantResponse(b6f["cyclic_atp_per_electron"]),
    )
    cyto_stoich = CytochromeStoichiometry{FT}(;
        electrons_per_carboxylation = stoich_cfg["electrons_per_carboxylation"],
        electrons_per_photorespiratory_co2 = stoich_cfg["electrons_per_photorespiratory_co2"],
        atp_per_carboxylation = stoich_cfg["atp_per_carboxylation"],
        atp_per_photorespiratory_co2 = stoich_cfg["atp_per_photorespiratory_co2"],
    )
    cyto_et = CytochromeElectronParameters(;
        psii_fraction = b6f["psii_fraction"],
        psi_maximum_yield = psi_maximum_photochemical_yield(jb_connected),
        limitation = SerialLimitation(),
    )
    cyto_photo = CytochromePhotosynthesis(;
        temperature = cyto_temp,
        electron_transport = cyto_et,
        stoichiometry = cyto_stoich,
        limitation = MinimumLimitation(),
        compensation = SignedAssimilation(),
    )

    x = settings["experiment"]
    for model in (jb_connected, jb_puddle)
        validate_fluorescence_coupling(model, cyto_photo, x["temperature_K"], x["oxygen_Pa"])
    end

    # 3. Empirical van der Tol et al. (2014) Cotton & Drought KN(x) models on FvCB.
    fvcb_temp = C3TemperatureResponse(;
        vcmax = scale_arrhenius_like(base_temp.vcmax, vcmax_factor),
        jmax = scale_arrhenius_like(base_temp.jmax, transport_factor),
        respiration = base_temp.respiration,
        kc = base_temp.kc,
        ko = base_temp.ko,
        gamma_star = base_temp.gamma_star,
    )
    fvcb_et = ElectronTransportParameters(;
        psii_fraction = emp_rates["psii_fraction"],
        psii_yield = emp_rates["psii_low_light_yield"],
        limitation = base_leaf.photosynthesis.electron_transport.limitation,
    )
    fvcb_photo = C3Photosynthesis(;
        temperature = fvcb_temp,
        electron_transport = fvcb_et,
        stoichiometry = base_leaf.photosynthesis.stoichiometry,
        limitation = MinimumLimitation(),
        compensation = SignedAssimilation(),
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

    s = core["stomata"]
    medlyn = SteadyStateStomata(Medlyn{FT}(;
        minimum_conductance = s["minimum_conductance"],
        slope = s["medlyn_slope"],
        diffusivity_ratio = core["diffusion"]["diffusivity_ratio"],
        micromole_to_mole = s["micromole_to_mole"],
        vpd_floor = s["vpd_floor"],
    ))
    solver = core_solver(FT, core)

    return (;
        cyto_leaf = LeafModel(; photosynthesis = cyto_photo, stomata = medlyn, oxygen = base_leaf.oxygen),
        fvcb_leaf = LeafModel(; photosynthesis = fvcb_photo, stomata = medlyn, oxygen = base_leaf.oxygen),
        jb_connected,
        jb_puddle,
        cotton_fl,
        drought_fl,
        diffusion,
        solver,
        settings,
    )
end

function evaluate_fluorescence_comparison(case, site::LeafEnvironment)
    cyto_carbon = photosynthesize(case.cyto_leaf, site)
    fvcb_carbon = photosynthesize(case.fvcb_leaf, site)
    cyto_inputs = photosystem_flux_inputs(case.cyto_leaf.photosynthesis, site)
    fvcb_inputs = psii_flux_inputs(case.fvcb_leaf.photosynthesis, site)
    jb = fluorescence(case.jb_connected, cyto_inputs)
    jb_pud = fluorescence(case.jb_puddle, cyto_inputs)
    cot = fluorescence(case.cotton_fl, fvcb_inputs)
    drt = fluorescence(case.drought_fl, fvcb_inputs)
    all(fluorescence_valid, (jb, jb_pud, cot, drt)) ||
        error("Invalid fluorescence comparison state: $(jb.status), $(cot.status), $(drt.status)")
    val_or_gap(result, field) =
        result.status == FluorescenceDark ? nothing : getproperty(result, field)
    return (;
        absorbed_par = site.absorbed_par,
        biochemical_co2 = site.co2,
        jb_net_assimilation = cyto_carbon.net_assimilation,
        fvcb_net_assimilation = fvcb_carbon.net_assimilation,
        jb_psi_to_psii_ratio = cyto_carbon.psi_to_psii_ratio,
        jb_psii_yield = val_or_gap(jb, :psii_fluorescence_yield),
        jb_puddle_psii_yield = val_or_gap(jb_pud, :psii_fluorescence_yield),
        jb_psi_yield = val_or_gap(jb, :psi_fluorescence_yield),
        jb_total_yield = val_or_gap(jb, :total_fluorescence_yield),
        jb_psii_open_fraction = val_or_gap(jb, :psii_open_fraction),
        jb_npq = val_or_gap(jb, :npq),
        jb_total_photons = jb.total_photon_flux,
        vdt_cotton_yield = val_or_gap(cot, :fluorescence_yield),
        vdt_drought_yield = val_or_gap(drt, :fluorescence_yield),
        vdt_cotton_photons = cot.photon_flux,
        vdt_drought_photons = drt.photon_flux,
        dark = jb.status == FluorescenceDark,
    )
end
