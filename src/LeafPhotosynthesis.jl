"""
    LeafPhotosynthesis

Provide composable, zero-dependency scalar models for C3 leaf photosynthesis
(Farquhar–von Caemmerer–Berry and Johnson & Berry 2021 Cytochrome b₆f + PSI/PSII
coupling), analytical and equilibrium stomatal CO₂/H₂O gas exchange, oxygen
budgets, LHCII state transitions (`qT`), and steady-state chlorophyll and PAM
fluorescence across `Float32`, `Float64`, GPU kernels, and automatic differentiation.
"""
module LeafPhotosynthesis

# What a user writes to build and run a leaf model: components, environments, light
# inputs, and the entry points.
export LeafModel, C3Photosynthesis, CytochromePhotosynthesis,
    ConstantResponse, Arrhenius, PeakedArrhenius, C3TemperatureResponse,
    CytochromeTemperatureResponse, CytochromeCapacity,
    PrescribedCompensationPoint, OxygenScaledCompensationPoint,
    ElectronTransportParameters, CytochromeElectronParameters,
    C3Stoichiometry, CytochromeStoichiometry,
    MinimumLimitation, SerialLimitation, QuadraticLimitation,
    SignedAssimilation, NonnegativeAssimilation

export LeafEnvironment, GasExchangeEnvironment, StomatalState,
    BroadbandExcitation, PrescribedPSIIExcitation, PrescribedDualExcitation

export BallBerry, Medlyn, SteadyStateStomata, PrognosticStomata, NoStomatalDynamics,
    SeriesDiffusion, WaterVaporDiffusion, CoreOxygenExchange, NoOxygenExchange,
    EquilibriumSolver, Bisection, SafeguardedNewton

export NoFluorescence, RateCompetitionFluorescence, ConstantQuenching, SaturationQuenching,
    CytochromeFluorescence, PSIIRateParameters, PSIRateParameters,
    StateTransitionRule, PAMFluorescenceObserver

export photosynthesize, biochemical_turnover, oxygen_exchange, gas_exchange,
    leaf_tendencies, equilibrium_gas_exchange, assimilation_sensitivity,
    transpiration_sensitivity, psii_flux_inputs, photosystem_flux_inputs, fluorescence,
    fluorescence_valid, validate_fluorescence_coupling, state_transition_excitation,
    pam_fluorescence, pam_fluorescence_valid

# Supported API that is rarely typed: abstract types for new components, component-level
# functions, rate and input bundles, result types and status values. Qualify these as
# `LeafPhotosynthesis.name` or import them with `using LeafPhotosynthesis: name`. The
# `public` keyword parses from Julia 1.11 on.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("""
    public AbstractLeafScalar, AbstractPhotosynthesis, AbstractC3Photosynthesis,
        AbstractCompensationPoint, AbstractExcitation, AbstractFluorescence,
        AbstractQuenching, AbstractStomatalTarget, AbstractStomatalModel,
        AbstractOxygenExchange, AbstractLeafDiffusion, AbstractEquilibriumMethod,
        temperature_response, biochemical_rates, compensation_point, C3Rates, CytochromeRates,
        electron_transport, cytochrome_electron_supply, valid_excitation, psii_excitation,
        psi_excitation, cytochrome_cyclic_coefficients, cytochrome_carbon_coefficients,
        psi_to_psii_electron_ratio, colimit, carbon_assimilation,
        carbon_assimilation_derivative, StomatalDrivers, target_conductance,
        stomatal_tendency, OxygenFluxInputs,
        EquilibriumSolution, EquilibriumStatus, EquilibriumConverged,
        EquilibriumIterationLimit, EquilibriumPrecisionLimit, EquilibriumInvalidInput,
        EquilibriumInvalidBracket,
        PSIIFluxInputs, CoupledPhotosystemInputs, RateCompetitionResult,
        CytochromeFluorescenceResult, PAMFluorescenceResult, FluorescenceStatus,
        FluorescenceLit, FluorescenceDark, FluorescenceInvalidInput,
        FluorescenceIncompatibleRates, quenching_rate, maximum_photochemical_yield,
        psi_maximum_photochemical_yield, psii_maximum_photochemical_yield,
        fluorescence_low_light_yield, minimum_psii_partition,
        dark_balanced_psii_partition, state_transition_partition, dark_pam_levels
    """))
end

"""
    AbstractLeafScalar

Supertype of every component, parameter, input and result struct of LeafPhotosynthesis.

Its subtypes are immutable values that broadcast as scalars, so
`photosynthesize.(leaf, environments)` needs no `Ref(leaf)`. A new struct gets this by
subtyping `AbstractLeafScalar` or one of the package's abstract component types.
"""
abstract type AbstractLeafScalar end
Base.broadcastable(x::AbstractLeafScalar) = Ref(x)

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

end
