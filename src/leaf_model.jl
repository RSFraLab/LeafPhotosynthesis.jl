"""
    LeafModel(; photosynthesis, stomata=NoStomatalDynamics(), oxygen=NoOxygenExchange())

Small leaf composition boundary with concrete photosynthesis, optional stomatal,
and optional oxygen components. Explicit no-op types disable optional processes;
these structural defaults introduce no scientific parameter values. State,
timestepping, arrays, canopy transport, and soils remain external.
"""
struct LeafModel{P<:AbstractPhotosynthesis,S<:AbstractStomatalModel,O<:AbstractOxygenExchange}
    photosynthesis::P
    stomata::S
    oxygen::O
end
LeafModel(; photosynthesis, stomata = NoStomatalDynamics(), oxygen = NoOxygenExchange()) = LeafModel(photosynthesis, stomata, oxygen)

@inline photosynthesize(leaf::LeafModel, env::LeafEnvironment) =
    photosynthesize(leaf.photosynthesis, env)

@inline oxygen_exchange(leaf::LeafModel, env::LeafEnvironment) =
    oxygen_exchange(leaf.oxygen, leaf.photosynthesis, env)

@inline stomatal_tendency(leaf::LeafModel, state::StomatalState, drivers::StomatalDrivers) =
    stomatal_tendency(leaf.stomata, state, drivers)

@inline gas_exchange(leaf::LeafModel, env::GasExchangeEnvironment, state::StomatalState, diffusion::AbstractLeafDiffusion) =
    gas_exchange(leaf.photosynthesis, env, state, diffusion)

@inline equilibrium_gas_exchange(leaf::LeafModel{P,S}, env, diffusion, solver; initial_co2 = nothing) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(leaf.photosynthesis, leaf.stomata, env, diffusion, solver; initial_co2)

"""
    leaf_tendencies(leaf, environment, stomatal_state, diffusion)

Evaluate prescribed-conductance gas exchange, the prognostic stomatal tendency,
and optional oxygen diagnostics. Returns values only; the caller advances g_sw.
Requires PrognosticStomata or NoStomatalDynamics and the supported analytical C3 closure.
The instantaneous target depends on current assimilation; it is not the already
converged equilibrium of the fully coupled leaf system.
"""
@inline function leaf_tendencies(leaf::LeafModel, env::GasExchangeEnvironment, state::StomatalState,
    diffusion::AbstractLeafDiffusion)
    exchange = gas_exchange(leaf, env, state, diffusion)
    drivers = StomatalDrivers(exchange.carbon.net_assimilation, exchange.surface_co2,
        env.pressure, env.relative_humidity, env.vpd)
    site = LeafEnvironment(env.temperature, env.absorbed_par, exchange.biochemical_co2, env.oxygen)
    return (; gas_exchange = exchange,
        stomatal_outputs(leaf.stomata, state, drivers)...,
        oxygen = oxygen_exchange(leaf.oxygen, leaf.photosynthesis, site))
end
