"""
    LeafModel(; photosynthesis, stomata=NoStomatalDynamics(), oxygen=NoOxygenExchange())
    LeafModel(photosynthesis, stomata, oxygen)

Compose concrete `photosynthesis` ([`AbstractPhotosynthesis`](@ref)), `stomata`
([`AbstractStomatalModel`](@ref)), and `oxygen` ([`AbstractOxygenExchange`](@ref))
components into an immutable leaf model.

Explicit no-op defaults (`NoStomatalDynamics()`, `NoOxygenExchange()`) disable
optional processes without introducing hidden scientific parameter values.
"""
struct LeafModel{P<:AbstractPhotosynthesis,S<:AbstractStomatalModel,O<:AbstractOxygenExchange} <:
       AbstractLeafScalar
    photosynthesis::P
    stomata::S
    oxygen::O
end
LeafModel(; photosynthesis, stomata = NoStomatalDynamics(), oxygen = NoOxygenExchange()) =
    LeafModel(photosynthesis, stomata, oxygen)

@inline photosynthesize(
    leaf::LeafModel,
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = photosynthesize(leaf.photosynthesis, env, excitation)

@inline biochemical_turnover(
    leaf::LeafModel,
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = biochemical_turnover(leaf.photosynthesis, env, excitation)

@inline oxygen_exchange(
    leaf::LeafModel,
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = oxygen_exchange(leaf.oxygen, leaf.photosynthesis, env, excitation)

@inline stomatal_tendency(leaf::LeafModel, state::StomatalState, drivers::StomatalDrivers) =
    stomatal_tendency(leaf.stomata, state, drivers)

@inline gas_exchange(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = gas_exchange(leaf.photosynthesis, env, state, diffusion, excitation)

@inline equilibrium_gas_exchange(
    leaf::LeafModel{P,S}, env, diffusion, solver,
    excitation::AbstractExcitation = BroadbandExcitation(); initial_co2 = nothing,
) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(
        leaf.photosynthesis, leaf.stomata, env, diffusion, solver, excitation; initial_co2,
    )

"""
    leaf_tendencies(leaf, environment, stomatal_state, diffusion,
        excitation=BroadbandExcitation())

Evaluate prescribed-conductance [`gas_exchange`](@ref), prognostic stomatal
target and tendency, and optional [`oxygen_exchange`](@ref) diagnostics at
`stomatal_state` without mutating state.
"""
@inline function leaf_tendencies(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    exchange = gas_exchange(leaf, env, state, diffusion, excitation)
    drivers = StomatalDrivers(
        exchange.carbon.net_assimilation, exchange.surface_co2,
        env.pressure, env.relative_humidity, env.vpd,
    )
    site = LeafEnvironment(env.temperature, env.absorbed_par, exchange.biochemical_co2, env.oxygen)
    return (;
        gas_exchange = exchange,
        stomatal_outputs(leaf.stomata, state, drivers)...,
        oxygen = oxygen_exchange(leaf, site, excitation),
    )
end
