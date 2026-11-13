"""
    LeafPhotosynthesis

Provide composable, zero-dependency scalar models for C3 leaf photosynthesis
(Farquhar–von Caemmerer–Berry and Johnson & Berry 2021 Cytochrome b₆f + PSI/PSII
coupling), analytical and equilibrium stomatal CO₂/H₂O gas exchange, oxygen
budgets, LHCII state transitions (`qT`), and steady-state chlorophyll and PAM
fluorescence across `Float32`, `Float64`, GPU kernels, and automatic differentiation.
"""
module LeafPhotosynthesis

export AbstractPhotosynthesis, C3Photosynthesis, CytochromePhotosynthesis,
    LeafModel, LeafEnvironment,
    ConstantResponse, Arrhenius, PeakedArrhenius, temperature_response,
    C3TemperatureResponse, C3Rates,
    CytochromeCapacity, CytochromeTemperatureResponse, CytochromeRates,
    biochemical_rates, ElectronTransportParameters, electron_transport,
    CytochromeElectronParameters, cytochrome_electron_supply,
    C3Stoichiometry, CytochromeStoichiometry,
    cytochrome_cyclic_coefficients, cytochrome_carbon_coefficients, psi_to_psii_electron_ratio,
    MinimumLimitation, SerialLimitation, QuadraticLimitation, colimit,
    SignedAssimilation, NonnegativeAssimilation, carbon_assimilation, photosynthesize,
    biochemical_turnover

export AbstractCompensationPoint, PrescribedCompensationPoint,
    OxygenScaledCompensationPoint, compensation_point

export AbstractExcitation, BroadbandExcitation, PrescribedPSIIExcitation,
    PrescribedDualExcitation, valid_excitation, psii_excitation, psi_excitation

export AbstractFluorescence, NoFluorescence, AbstractQuenching, ConstantQuenching,
    SaturationQuenching, RateCompetitionFluorescence, PSIIFluxInputs,
    RateCompetitionResult, PSIIRateParameters, PSIRateParameters,
    CytochromeFluorescence, CoupledPhotosystemInputs, CytochromeFluorescenceResult,
    FluorescenceStatus, FluorescenceLit, FluorescenceDark,
    FluorescenceInvalidInput, FluorescenceIncompatibleRates,
    fluorescence, fluorescence_valid, quenching_rate, maximum_photochemical_yield,
    psi_maximum_photochemical_yield, psii_maximum_photochemical_yield,
    minimum_psii_partition, psii_flux_inputs, photosystem_flux_inputs,
    fluorescence_low_light_yield, validate_fluorescence_coupling,
    StateTransitionRule, dark_balanced_psii_partition,
    state_transition_partition, state_transition_excitation,
    PAMFluorescenceObserver, PAMFluorescenceResult,
    pam_fluorescence_valid, dark_pam_levels, pam_fluorescence

export AbstractStomatalTarget, AbstractStomatalModel, BallBerry, Medlyn,
    SteadyStateStomata, PrognosticStomata, NoStomatalDynamics, StomatalState, StomatalDrivers,
    target_conductance, stomatal_tendency,
    AbstractOxygenExchange, CoreOxygenExchange, NoOxygenExchange, OxygenFluxInputs, oxygen_exchange

export AbstractLeafDiffusion, SeriesDiffusion, GasExchangeEnvironment,
    gas_exchange, leaf_tendencies

export carbon_assimilation_derivative, assimilation_sensitivity,
    WaterVaporDiffusion, transpiration_sensitivity,
    AbstractEquilibriumMethod, Bisection, SafeguardedNewton, EquilibriumSolver,
    EquilibriumSolution, EquilibriumStatus, EquilibriumConverged,
    EquilibriumIterationLimit, EquilibriumPrecisionLimit, EquilibriumInvalidInput,
    EquilibriumInvalidBracket, equilibrium_gas_exchange

# Follow this list from physical responses to a complete prescribed-CO₂ leaf call.
include("temperature_response.jl")
include("photorespiratory_compensation.jl")
include("rate_limitation.jl")
include("excitation.jl")
include("electron_transport.jl")
include("carbon_assimilation.jl")
include("assimilation_derivatives.jl")
include("photosynthesis.jl")
include("fluorescence.jl")
include("stomatal_conductance.jl")
include("gas_diffusion.jl")
include("equilibrium_solver.jl")
include("gas_exchange_sensitivities.jl")
include("oxygen_exchange.jl")
include("cytochrome_photosynthesis.jl")
include("cytochrome_fluorescence.jl")
include("leaf_model.jl")
include("cytochrome_state_transitions.jl")
include("pam_fluorescence.jl")

const _ScalarBroadcastable = Union{
    LeafModel,
    AbstractPhotosynthesis,
    AbstractStomatalModel,
    AbstractStomatalTarget,
    AbstractOxygenExchange,
    AbstractLeafDiffusion,
    AbstractEquilibriumMethod,
    EquilibriumSolver,
    AbstractFluorescence,
    AbstractQuenching,
    AbstractExcitation,
    AbstractCompensationPoint,
    Arrhenius,
    PeakedArrhenius,
    ConstantResponse,
    C3TemperatureResponse,
    C3Rates,
    CytochromeCapacity,
    CytochromeTemperatureResponse,
    CytochromeRates,
    ElectronTransportParameters,
    CytochromeElectronParameters,
    C3Stoichiometry,
    CytochromeStoichiometry,
    MinimumLimitation,
    SerialLimitation,
    QuadraticLimitation,
    SignedAssimilation,
    NonnegativeAssimilation,
    LeafEnvironment,
    GasExchangeEnvironment,
    StomatalState,
    StomatalDrivers,
    OxygenFluxInputs,
    PSIIFluxInputs,
    CoupledPhotosystemInputs,
    PSIIRateParameters,
    PSIRateParameters,
    StateTransitionRule,
    PAMFluorescenceObserver,
    RateCompetitionResult,
    CytochromeFluorescenceResult,
    PAMFluorescenceResult,
    EquilibriumSolution,
}

Base.broadcastable(x::_ScalarBroadcastable) = Ref(x)

end
