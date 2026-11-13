"""
    LeafModel(; photosynthesis, stomata=NoStomatalDynamics(), oxygen=NoOxygenExchange())
    LeafModel(photosynthesis, stomata, oxygen)

Compose concrete `photosynthesis` ([`AbstractPhotosynthesis`](@ref)), `stomata`
([`AbstractStomatalModel`](@ref)), and `oxygen` ([`AbstractOxygenExchange`](@ref))
components into an immutable leaf model.

Explicit no-op defaults (`NoStomatalDynamics()`, `NoOxygenExchange()`) disable
optional processes without introducing hidden scientific parameter values.
"""
struct LeafModel{P<:AbstractPhotosynthesis,S<:AbstractStomatalModel,O<:AbstractOxygenExchange}
    photosynthesis::P
    stomata::S
    oxygen::O
end
LeafModel(; photosynthesis, stomata = NoStomatalDynamics(), oxygen = NoOxygenExchange()) =
    LeafModel(photosynthesis, stomata, oxygen)

@inline photosynthesize(leaf::LeafModel, env::LeafEnvironment) =
    photosynthesize(leaf.photosynthesis, env)

@inline biochemical_turnover(leaf::LeafModel, env::LeafEnvironment) =
    biochemical_turnover(leaf.photosynthesis, env)

@inline oxygen_exchange(leaf::LeafModel, env::LeafEnvironment) =
    oxygen_exchange(leaf.oxygen, leaf.photosynthesis, env)

@inline photosynthesize(leaf::LeafModel, env::LeafEnvironment, excitation::AbstractExcitation) =
    photosynthesize(leaf.photosynthesis, env, excitation)
@inline biochemical_turnover(leaf::LeafModel, env::LeafEnvironment, excitation::AbstractExcitation) =
    biochemical_turnover(leaf.photosynthesis, env, excitation)
@inline oxygen_exchange(leaf::LeafModel, env::LeafEnvironment, excitation::AbstractExcitation) =
    oxygen_exchange(leaf.oxygen, leaf.photosynthesis, env, excitation)

# Broadband preserves legacy component entry points, including user-defined
# photosynthesis/oxygen implementations that only implement the two-argument arity.
@inline photosynthesize(leaf::LeafModel, env::LeafEnvironment, ::BroadbandExcitation) =
    photosynthesize(leaf, env)
@inline biochemical_turnover(leaf::LeafModel, env::LeafEnvironment, ::BroadbandExcitation) =
    biochemical_turnover(leaf, env)
@inline oxygen_exchange(leaf::LeafModel, env::LeafEnvironment, ::BroadbandExcitation) =
    oxygen_exchange(leaf, env)

@inline stomatal_tendency(leaf::LeafModel, state::StomatalState, drivers::StomatalDrivers) =
    stomatal_tendency(leaf.stomata, state, drivers)

@inline gas_exchange(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
) = gas_exchange(leaf.photosynthesis, env, state, diffusion)

@inline gas_exchange(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    excitation::AbstractExcitation,
) = gas_exchange(leaf.photosynthesis, env, state, diffusion, excitation)
@inline gas_exchange(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    ::BroadbandExcitation,
) = gas_exchange(leaf, env, state, diffusion)

@inline equilibrium_gas_exchange(
    leaf::LeafModel{P,S}, env, diffusion, solver; initial_co2 = nothing,
) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(leaf.photosynthesis, leaf.stomata, env, diffusion, solver; initial_co2)

@inline equilibrium_gas_exchange(
    leaf::LeafModel{P,S}, env, diffusion, solver,
    excitation::AbstractExcitation; initial_co2 = nothing,
) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(
        leaf.photosynthesis, leaf.stomata, env, diffusion, solver, excitation; initial_co2,
    )
@inline equilibrium_gas_exchange(
    leaf::LeafModel{P,S}, env, diffusion, solver,
    ::BroadbandExcitation; initial_co2 = nothing,
) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(leaf, env, diffusion, solver; initial_co2)

"""
    leaf_tendencies(leaf, environment, stomatal_state, diffusion[, excitation])

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
