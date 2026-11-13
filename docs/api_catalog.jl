# Canonical biophysical-process organization of the LeafPhotosynthesis public API.
# Consumed by both `docs/build_api_docs.jl` (GitHub/Pages Markdown under `docs/api/`)
# and `docs/make.jl` (local Documenter.jl HTML build with `checkdocs = :exports`).

const API_PAGES = [
    (
        slug = "leaf_composition_and_environments",
        title = "Leaf Composition & Environments",
        walkthroughs = [
            ("Design decisions", "../design.md"),
            ("Core gas exchange", "../core_exchange.md"),
        ],
        summary = """
        `LeafModel` assembles concrete photosynthesis, stomatal, and oxygen components
        behind abstract process interfaces, using explicit no-op types
        (`NoStomatalDynamics`, `NoOxygenExchange`) when optional budgets are omitted.
        Scalar leaf-site (`LeafEnvironment`) and boundary-air (`GasExchangeEnvironment`)
        containers carry forcing in explicit SI/molar units (`K`, `Pa`,
        `μmol m⁻² leaf s⁻¹`, `mol m⁻² leaf s⁻¹`) across `Float32`, `Float64`, GPU
        broadcasts, and `ForwardDiff.Dual` tangents.
        """,
        sections = [
            "Leaf Assembly & Prescribed-CO₂ Evaluation" => [
                :LeafModel,
                :photosynthesize,
            ],
            "Scalar Forcing Environments" => [
                :LeafEnvironment,
                :GasExchangeEnvironment,
            ],
        ],
    ),
    (
        slug = "temperature_and_compensation",
        title = "Temperature Responses & Photorespiratory Compensation",
        walkthroughs = [
            ("Bonan 11.1 curves", "../bonan_exercises.md"),
            ("Oxygen response", "../oxygen_response.md"),
        ],
        summary = """
        Thermal kinetics evaluate scalar rate laws (`ConstantResponse`, `Arrhenius`,
        `PeakedArrhenius`) at leaf temperature `T` [K] without embedded reference
        temperatures or gas constants. Rate bundles assemble all temperature-dependent
        capacities and Michaelis–Menten constants for Farquhar–von Caemmerer–Berry C3
        biochemistry (`C3TemperatureResponse`, `C3Rates`) and Johnson & Berry (2021)
        Cytochrome \$b_6f\$ coupling (`CytochromeCapacity`,
        `CytochromeTemperatureResponse`, `CytochromeRates`), paired with prescribed or
        explicitly oxygen-scaled photorespiratory compensation points \$\\Gamma^*\$ [Pa].
        """,
        sections = [
            "Scalar Thermal Response Laws" => [
                :ConstantResponse,
                :Arrhenius,
                :PeakedArrhenius,
                :temperature_response,
            ],
            "Biochemical Rate Bundles" => [
                :C3TemperatureResponse,
                :C3Rates,
                :CytochromeCapacity,
                :CytochromeTemperatureResponse,
                :CytochromeRates,
                :biochemical_rates,
            ],
            "Photorespiratory Compensation Point" => [
                :AbstractCompensationPoint,
                :PrescribedCompensationPoint,
                :OxygenScaledCompensationPoint,
                :compensation_point,
            ],
        ],
    ),
    (
        slug = "excitation_and_electron_transport",
        title = "Excitation, Electron Transport & Stoichiometry",
        walkthroughs = [
            ("Excitation contract", "../excitation.md"),
            ("Cytochrome b6f mechanics", "../cytochrome_fluorescence.md"),
        ],
        summary = """
        Excitation types (`BroadbandExcitation`, `PrescribedPSIIExcitation`,
        `PrescribedDualExcitation`) pass absorbed photon fluxes from external leaf/canopy
        optics into photosystems II and I without embedding spectral pigment tables inside
        `LeafModel`. Potential electron supply is evaluated either from empirical
        non-rectangular hyperbola parameters (`ElectronTransportParameters`,
        `electron_transport`) or from Cytochrome \$b_6f\$ series kinetics with linear and
        cyclic ATP/NADPH stoichiometry (`CytochromeElectronParameters`,
        `CytochromeStoichiometry`).
        """,
        sections = [
            "Photosystem Excitation Inputs" => [
                :AbstractExcitation,
                :BroadbandExcitation,
                :PrescribedPSIIExcitation,
                :PrescribedDualExcitation,
                :valid_excitation,
                :psii_excitation,
                :psi_excitation,
            ],
            "Potential Electron Supply" => [
                :ElectronTransportParameters,
                :electron_transport,
                :CytochromeElectronParameters,
                :cytochrome_electron_supply,
            ],
            "ATP/NADPH & Cyclic Electron Stoichiometry" => [
                :C3Stoichiometry,
                :CytochromeStoichiometry,
                :cytochrome_cyclic_coefficients,
                :cytochrome_carbon_coefficients,
                :psi_to_psii_electron_ratio,
            ],
        ],
    ),
    (
        slug = "carbon_assimilation_and_turnover",
        title = "Carbon Assimilation, Rate Limitation & Turnover",
        walkthroughs = [
            ("Biochemical turnover", "../biochemical_turnover.md"),
            ("Bonan 11.1 curves", "../bonan_exercises.md"),
        ],
        summary = """
        Photosynthetic demand components (`C3Photosynthesis`, `CytochromePhotosynthesis`)
        combine Rubisco-limited and electron-transport-limited assimilation via explicit
        colimitation rules (`MinimumLimitation`, `SerialLimitation`,
        `QuadraticLimitation`) and respiration sign policies (`SignedAssimilation`,
        `NonnegativeAssimilation`). `biochemical_turnover` diagnoses gross carboxylation,
        photorespiratory oxygenation, and realized linear electron flow separately from
        light-limited electron capacity, remaining regular across the CO₂ compensation
        point.
        """,
        sections = [
            "Photosynthesis Demand Components" => [
                :AbstractPhotosynthesis,
                :C3Photosynthesis,
                :CytochromePhotosynthesis,
            ],
            "Rate Colimitation & Sign Policies" => [
                :MinimumLimitation,
                :SerialLimitation,
                :QuadraticLimitation,
                :colimit,
                :SignedAssimilation,
                :NonnegativeAssimilation,
            ],
            "Carbon Demand & Realized Biochemical Turnover" => [
                :carbon_assimilation,
                :biochemical_turnover,
            ],
        ],
    ),
    (
        slug = "fluorescence_and_pam",
        title = "Chlorophyll Fluorescence, State Transitions & PAM",
        walkthroughs = [
            ("PSII fluorescence", "../fluorescence.md"),
            ("Cytochrome b6f stress fluorescence", "../cytochrome_fluorescence.md"),
        ],
        summary = """
        Internal photosystem fluorescence yields and emitted photon fluxes are computed
        from the same photosystem excitations and realized electron flows used by carbon
        assimilation, either through empirical PSII rate competition with replaceable
        non-photochemical quenching (`RateCompetitionFluorescence`, `ConstantQuenching`,
        `SaturationQuenching`) or through coupled PSII/PSI + Cytochrome \$b_6f\$ redox
        feedback (`CytochromeFluorescence`). Closed-form LHCII state transitions (`qT`,
        `StateTransitionRule`) and apparent Pulse-Amplitude-Modulated fluorometry with
        PSI detector leakage (`PAMFluorescenceObserver`) compose directly on top of the
        steady-state photosystem states.
        """,
        sections = [
            "Fluorescence Models & Quenching Laws" => [
                :AbstractFluorescence,
                :NoFluorescence,
                :RateCompetitionFluorescence,
                :CytochromeFluorescence,
                :AbstractQuenching,
                :ConstantQuenching,
                :SaturationQuenching,
                :quenching_rate,
            ],
            "Photosystem Rate Parameters & Flux Inputs" => [
                :PSIIRateParameters,
                :PSIRateParameters,
                :PSIIFluxInputs,
                :CoupledPhotosystemInputs,
                :psii_flux_inputs,
                :photosystem_flux_inputs,
                :maximum_photochemical_yield,
                :psii_maximum_photochemical_yield,
                :psi_maximum_photochemical_yield,
                :minimum_psii_partition,
                :fluorescence_low_light_yield,
                :validate_fluorescence_coupling,
            ],
            "Steady-State Fluorescence Evaluation & Status Codes" => [
                :RateCompetitionResult,
                :CytochromeFluorescenceResult,
                :FluorescenceStatus,
                :FluorescenceLit,
                :FluorescenceDark,
                :FluorescenceInvalidInput,
                :FluorescenceIncompatibleRates,
                :fluorescence,
                :fluorescence_valid,
            ],
            "Closed-Form LHCII State Transitions (qT)" => [
                :StateTransitionRule,
                :dark_balanced_psii_partition,
                :state_transition_partition,
                :state_transition_excitation,
            ],
            "Pulse-Amplitude-Modulated (PAM) Fluorometry" => [
                :PAMFluorescenceObserver,
                :PAMFluorescenceResult,
                :pam_fluorescence,
                :pam_fluorescence_valid,
                :dark_pam_levels,
            ],
        ],
    ),
    (
        slug = "stomata_diffusion_and_equilibrium",
        title = "Stomatal Conductance, Gas Diffusion, Equilibrium & Oxygen",
        walkthroughs = [
            ("Core exchange", "../core_exchange.md"),
            ("Equilibrium & sensitivities", "../equilibrium_and_sensitivities.md"),
            ("Oxygen response", "../oxygen_response.md"),
        ],
        summary = """
        Gas-exchange components separate stomatal target laws (`BallBerry`, `Medlyn`) and
        time-evolution semantics (`SteadyStateStomata`, `PrognosticStomata`,
        `NoStomatalDynamics`) from boundary+stomatal series diffusion (`SeriesDiffusion`,
        `WaterVaporDiffusion`), safeguarded scalar \$A\\text{–}C_i\\text{–}g_{sw}\$
        equilibrium root-finding (`Bisection`, `SafeguardedNewton`, `EquilibriumSolver`),
        analytic fixed-temperature gas-exchange sensitivities, and explicit gross/net O₂
        budgets (`CoreOxygenExchange`, `NoOxygenExchange`).
        """,
        sections = [
            "Stomatal Conductance Targets & Dynamics" => [
                :AbstractStomatalTarget,
                :BallBerry,
                :Medlyn,
                :AbstractStomatalModel,
                :SteadyStateStomata,
                :PrognosticStomata,
                :NoStomatalDynamics,
                :StomatalState,
                :StomatalDrivers,
                :target_conductance,
                :stomatal_tendency,
            ],
            "Gas Diffusion & Coupled Leaf Tendencies" => [
                :AbstractLeafDiffusion,
                :SeriesDiffusion,
                :WaterVaporDiffusion,
                :gas_exchange,
                :leaf_tendencies,
            ],
            "Coupled A–Ci–gsw Equilibrium Solver" => [
                :AbstractEquilibriumMethod,
                :Bisection,
                :SafeguardedNewton,
                :EquilibriumSolver,
                :EquilibriumSolution,
                :EquilibriumStatus,
                :EquilibriumConverged,
                :EquilibriumIterationLimit,
                :EquilibriumPrecisionLimit,
                :EquilibriumInvalidInput,
                :EquilibriumInvalidBracket,
                :equilibrium_gas_exchange,
            ],
            "Fixed-Temperature Gas-Exchange Sensitivities" => [
                :carbon_assimilation_derivative,
                :assimilation_sensitivity,
                :transpiration_sensitivity,
            ],
            "Leaf Oxygen Budget" => [
                :AbstractOxygenExchange,
                :CoreOxygenExchange,
                :NoOxygenExchange,
                :OxygenFluxInputs,
                :oxygen_exchange,
            ],
        ],
    ),
]
