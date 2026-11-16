# Loaded first by every test group (see runtests.jl) and by gpu/runtests.jl, each time into a
# fresh module, so no group depends on another group having run.
using Test
using LeafPhotosynthesis
using LeafPhotosynthesis: AbstractPhotosynthesis, AbstractExcitation, AbstractQuenching,
    AbstractStomatalModel, AbstractOxygenExchange, temperature_response, biochemical_rates,
    compensation_point, CytochromeRates, electron_transport, valid_excitation,
    psii_excitation, cytochrome_carbon_coefficients, colimit, carbon_assimilation,
    carbon_assimilation_derivative, StomatalDrivers, target_conductance, stomatal_tendency,
    OxygenFluxInputs, EquilibriumConverged, EquilibriumIterationLimit,
    EquilibriumPrecisionLimit, EquilibriumInvalidInput, EquilibriumInvalidBracket,
    EquilibriumStomataClosed, PSIIFluxInputs,
    CoupledPhotosystemInputs, RateCompetitionResult, CytochromeFluorescenceResult,
    PAMFluorescenceResult, FluorescenceLit, FluorescenceDark, FluorescenceInvalidInput,
    FluorescenceIncompatibleRates, quenching_rate, maximum_photochemical_yield,
    PSIIOpenFractionInputs, ObservedFluorescenceInputs, LakeModelResult, lake_model_yields, C3Rates,
    psi_maximum_photochemical_yield, psii_maximum_photochemical_yield,
    minimum_psii_partition, dark_balanced_psii_partition, state_transition_partition,
    dark_pam_levels

# The example code is the test fixture. The physics tests build their leaves, parameter sets
# and response curves with the functions of examples/, the same code that produces the
# documentation gallery, so a broken example fails the tests. This is the only place where
# test code includes files from examples/; the parameter values stay in examples/parameters/.
const EXAMPLES = joinpath(dirname(@__DIR__), "examples")
# read_bonan_parameters, bonan_leaf, bonan_environment, write_tables, bonan_curves,
# core_leaf, core_environment, core_solver, equilibrium_example
include(joinpath(EXAMPLES, "core_exchange.jl"))
# setup_oxygen_response, oxygen_response_environments, with_compensation_point
include(joinpath(EXAMPLES, "setup", "oxygen_response.jl"))
# setup_fluorescence, fluorescence_row
include(joinpath(EXAMPLES, "setup", "fluorescence.jl"))
# setup_cytochrome_stress, evaluate_fluorescence_comparison
include(joinpath(EXAMPLES, "setup", "cytochrome_stress.jl"))
# flexas_2002_curves, setup_flexas_2002, evaluate_flexas_point and, through its setup,
# setup_johnson_berry_2021 and exact_dual_excitation
include(joinpath(EXAMPLES, "flexas_2002.jl"))
# gu_2019_curves, setup_gu_2019, measured_open_fraction, gross_assimilation_from_electron_flow
include(joinpath(EXAMPLES, "gu_2019.jl"))

# Settings of Bonan's supplemental program 11.1, used by several groups.
const parameter_path = joinpath(EXAMPLES, "parameters", "bonan_11_1.toml")

# The carbon fluxes of a flat result (`photosynthesize`, `gas_exchange` or the `exchange` of
# an `EquilibriumSolution`), in the order of `carbon_assimilation`.
carbon_fluxes(x) = NamedTuple{(:rubisco_limited, :electron_limited, :assimilation,
    :respiration, :net_assimilation)}(x)
