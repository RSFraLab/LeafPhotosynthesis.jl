module LeafPhotosynthesis

export AbstractPhotosynthesis, C3Photosynthesis, LeafModel, LeafEnvironment,
    Arrhenius, PeakedArrhenius, temperature_response, C3TemperatureResponse,
    C3Rates, biochemical_rates, ElectronTransportParameters, electron_transport,
    C3Stoichiometry, MinimumLimitation, QuadraticLimitation, colimit,
    SignedAssimilation, NonnegativeAssimilation, carbon_assimilation, photosynthesize

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
include("rate_limitation.jl")
include("electron_transport.jl")
include("carbon_assimilation.jl")
include("assimilation_derivatives.jl")
include("photosynthesis.jl")
include("stomatal_conductance.jl")
include("gas_diffusion.jl")
include("equilibrium_solver.jl")
include("gas_exchange_sensitivities.jl")
include("oxygen_exchange.jl")
include("leaf_model.jl")

end
