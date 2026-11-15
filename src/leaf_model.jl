"""
    LeafModel(; photosynthesis, stomata=NoStomatalDynamics(), oxygen=NoOxygenExchange())
    LeafModel(photosynthesis, stomata, oxygen)

Compose concrete `photosynthesis` ([`AbstractPhotosynthesis`](@ref)), `stomata`
([`AbstractStomatalModel`](@ref)), and `oxygen` ([`AbstractOxygenExchange`](@ref))
components into an immutable leaf model.

Explicit no-op defaults (`NoStomatalDynamics()`, `NoOxygenExchange()`) disable
optional processes without introducing hidden scientific parameter values. Components of
different binary precisions are a `MethodError`.
"""
struct LeafModel{P<:AbstractPhotosynthesis,S<:AbstractStomatalModel,O<:AbstractOxygenExchange} <:
       AbstractLeafScalar
    photosynthesis::P
    stomata::S
    oxygen::O
    LeafModel{P,S,O}(fields::Vararg{Any,3}) where {P<:AbstractPhotosynthesis,S<:AbstractStomatalModel,O<:AbstractOxygenExchange} =
        new{P,S,O}(component_fields(LeafModel{P,S,O}, fields)...)
end
LeafModel(photosynthesis::AbstractPhotosynthesis, stomata::AbstractStomatalModel, oxygen::AbstractOxygenExchange) =
    LeafModel{typeof(photosynthesis),typeof(stomata),typeof(oxygen)}(photosynthesis, stomata, oxygen)
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

@inline fluorescence_inputs(
    leaf::LeafModel,
    site::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = fluorescence_inputs(leaf.photosynthesis, site, excitation)

@inline target_conductance(leaf::LeafModel, drivers::StomatalDrivers, diffusion::SeriesDiffusion) =
    target_conductance(leaf.stomata, drivers, diffusion)

@inline stomatal_tendency(leaf::LeafModel, state::StomatalState, drivers::StomatalDrivers, diffusion::SeriesDiffusion) =
    stomatal_tendency(leaf.stomata, state, drivers, diffusion)

@inline gas_exchange(
    leaf::LeafModel,
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
) = gas_exchange(leaf.photosynthesis, env, state, diffusion, excitation)

@inline equilibrium_gas_exchange(
    leaf::LeafModel{P,S}, env, diffusion, solver,
    excitation::AbstractExcitation = BroadbandExcitation(),
    initial_co2::Union{Nothing,Real} = nothing,
) where {P,S<:SteadyStateStomata} =
    equilibrium_gas_exchange(leaf.photosynthesis, leaf.stomata, env, diffusion, solver, excitation, initial_co2)

"""
    leaf_tendencies(leaf, environment, stomatal_state, diffusion,
        excitation=BroadbandExcitation())

Evaluate prescribed-conductance [`gas_exchange`](@ref), the stomatal target and tendency,
and the [`oxygen_exchange`](@ref) budget at `stomatal_state` without mutating state.

`leaf.stomata` is a [`PrognosticStomata`](@ref) or [`NoStomatalDynamics`](@ref); a
[`SteadyStateStomata`](@ref) leaf has no tendency (solve it with
[`equilibrium_gas_exchange`](@ref)) and is a `MethodError`.

Return `(; gas_exchange, conductance_target, conductance_tendency, oxygen)`:
- `gas_exchange`: the [`gas_exchange`](@ref) result at `stomatal_state`;
- `conductance_target` [mol H₂O m⁻² leaf s⁻¹] and `conductance_tendency`
  [mol H₂O m⁻² leaf s⁻²] of [`target_conductance`](@ref) and [`stomatal_tendency`](@ref),
  with drivers from the gas exchange (net assimilation, surface CO₂) and `environment`
  (pressure, relative humidity, VPD) and the diffusivity ratio of `diffusion`; for
  `NoStomatalDynamics` the target is the state conductance and the tendency zero;
- `oxygen`: the coupled [`oxygen_exchange`](@ref) at the solved `biochemical_co2`, or
  `nothing` for [`NoOxygenExchange`](@ref).

Each part follows its own domain: inputs outside the gas-exchange domain give `NaN` in every
field of `gas_exchange` and `oxygen` and, through the drivers, in the outputs of
`PrognosticStomata`; invalid stomatal drivers alone (for example relative humidity outside
`[0, 1]` with a Ball–Berry target, or a `NaN` VPD with a Medlyn target) give `NaN` stomatal
outputs only. The humidity driver the target does not read may be `NaN`.
"""
@inline function leaf_tendencies(
    leaf::LeafModel{<:AbstractPhotosynthesis,<:Union{PrognosticStomata,NoStomatalDynamics}},
    env::GasExchangeEnvironment,
    state::StomatalState,
    diffusion::AbstractLeafDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    require_compatible(leaf_tendencies, leaf, env, state, diffusion, excitation)
    exchange = gas_exchange(leaf, env, state, diffusion, excitation)
    drivers = StomatalDrivers(
        exchange.net_assimilation, exchange.surface_co2,
        env.pressure, env.relative_humidity, env.vpd,
    )
    site = LeafEnvironment(env.temperature, env.absorbed_par, exchange.biochemical_co2, env.oxygen)
    return (;
        gas_exchange = exchange,
        stomatal_outputs(leaf.stomata, state, drivers, diffusion)...,
        oxygen = oxygen_exchange(leaf, site, excitation),
    )
end
