"""
    AbstractPhotosynthesis

Supertype for leaf photosynthetic biochemistry components implementing [`photosynthesize`](@ref).
"""
abstract type AbstractPhotosynthesis end

"""
    AbstractC3Photosynthesis{L,C} <: AbstractPhotosynthesis

Supertype for C3 biochemistry that co-limits Rubisco carboxylation and a light-driven
electron supply with carbon co-limitation `L` and compensation treatment `C`.

[`C3Photosynthesis`](@ref) and [`CytochromePhotosynthesis`](@ref) are its subtypes. The
leaf entry points [`photosynthesize`](@ref), [`biochemical_turnover`](@ref),
[`gas_exchange`](@ref), [`equilibrium_gas_exchange`](@ref),
[`assimilation_sensitivity`](@ref) and the coupled [`oxygen_exchange`](@ref) are written
once for this supertype; a subtype supplies only how its electron supply limits carbon.

# Implementation

A subtype has the fields `temperature`, `electron_transport`, `stoichiometry`,
`limitation::L` and `compensation::C`, and implements these internal methods:

- `light_supply(p, environment, excitation) -> (rates, electrons, supply)`: kinetic
  rates at leaf temperature and O₂, the electron flux `electrons` that limits carbon
  [μmol e⁻ m⁻² leaf s⁻¹], and a NamedTuple `supply` of all available-supply scalars.
- `light_limited_kinetics(p, rates, electrons) -> (m, k)`: capacity `m`
  [μmol CO₂ m⁻² leaf s⁻¹] and half-saturation `k` [Pa] of the light-limited branch
  ``A_j = m (C_c - Γ^*) / (C_c + k)``.
- `electron_state(p, rates, supply, co2) -> (valid, electron_transport, extras)`:
  whether the electron state at finite biochemical-site `co2` lies in the model domain,
  the available linear electron transport reported to callers, and a NamedTuple
  of extra diagnostics appended to every leaf result.
- `valid_biochemistry(p, rates)` and `valid_equilibrium_biochemistry(p, rates)`: the
  parameter domains of the prescribed-CO₂ entry points and of the equilibrium solver.
- `electron_demand(rates, co2, stoichiometry)`, `capacity_and_slope` and the signed
  `biochemical_turnover` for the model's rate and stoichiometry types; the electron-
  limited capacity is `electrons / electron_demand` (see [`carbon_assimilation`](@ref)).
- `turnover_electron_flows(p, turnover)`: realized electron flows beyond
  `linear_electron_flow` that [`oxygen_exchange`](@ref) reports.
"""
abstract type AbstractC3Photosynthesis{L,C} <: AbstractPhotosynthesis end

"""
    C3Photosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation)

Compose Farquhar–von Caemmerer–Berry C3 biochemistry from concrete `temperature`,
`electron_transport`, `stoichiometry`, `limitation`, and `compensation` components.

`C3Photosynthesis` holds immutable parameters only, with no per-leaf arrays,
scratch buffers, or hidden empirical defaults.
"""
Base.@kwdef struct C3Photosynthesis{T,E,S,L,C} <: AbstractC3Photosynthesis{L,C}
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    compensation::C
end

"""
    LeafEnvironment{FT}(; temperature, absorbed_par, co2, oxygen)
    LeafEnvironment(temperature, absorbed_par, co2, oxygen)

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
and biochemical-site `co2` and `oxygen` partial pressures [Pa].
"""
Base.@kwdef struct LeafEnvironment{FT<:Real}
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
    function LeafEnvironment{FT}(temperature, absorbed_par, co2, oxygen) where {FT<:Real}
        require_continuous_scalar(
            FT, LeafEnvironment{FT}, (temperature, absorbed_par, co2, oxygen),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen)
    end
    function LeafEnvironment(
        temperature::FT,
        absorbed_par::FT,
        co2::FT,
        oxygen::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, LeafEnvironment, (temperature, absorbed_par, co2, oxygen),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen)
    end
end
LeafEnvironment(t::Real, par::Real, co2::Real, o2::Real) =
    LeafEnvironment(promote_compatible(LeafEnvironment, t, par, co2, o2)...)

# Physical domain of a biochemical site: positive finite temperature [K] and finite,
# nonnegative CO₂ and O₂ partial pressures [Pa]. Absorbed light is checked with the
# excitation by `valid_excitation`.
@inline function valid_site(temperature, co2, oxygen)
    t_p, co2_p, o2_p = map(primal_scalar, (temperature, co2, oxygen))
    return isfinite(t_p) && t_p > zero(t_p) &&
        isfinite(co2_p) && co2_p >= zero(co2_p) &&
        isfinite(o2_p) && o2_p >= zero(o2_p)
end
@inline valid_site(env) = valid_site(env.temperature, env.co2, env.oxygen)

# Domain of a leaf entry point at prescribed biochemical-site gases.
@inline valid_leaf(p, rates, env, excitation) =
    valid_excitation(excitation, env.absorbed_par) && valid_site(env) &&
    valid_biochemistry(p, rates)

@inline function light_supply(p::C3Photosynthesis, env, excitation)
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, rates.jmax, excitation)
    return rates, electrons, (; electron_transport = electrons)
end

@inline light_limited_kinetics(p::C3Photosynthesis, r, electrons) =
    (electrons / p.stoichiometry.electrons_per_carboxylation,
        p.stoichiometry.oxygenation_cost * r.gamma_star)

@inline electron_state(::C3Photosynthesis, r, supply, co2) =
    (true, supply.electron_transport, (;))

@inline valid_biochemistry(::C3Photosynthesis, r) = true

@inline function valid_equilibrium_biochemistry(p::C3Photosynthesis, r)
    s = p.stoichiometry
    return all(x -> isfinite(x) && x > zero(x), (s.electrons_per_carboxylation, s.oxygenation_cost)) &&
        all(isfinite, (r.vcmax, r.jmax, r.respiration, r.kc, r.ko, r.gamma_star)) &&
        min(primal_scalar(r.vcmax), primal_scalar(r.jmax), primal_scalar(r.respiration)) >= zero(primal_scalar(r.vcmax)) &&
        min(r.kc, r.ko, r.gamma_star) > zero(r.kc)
end

@inline turnover_electron_flows(::C3Photosynthesis, turnover) = (;)

# Output type of a leaf result: the promotion of parameters, forcing and electron supply.
@inline leaf_result_type(p, r, supply, forcing...) =
    typeof(r.vcmax + p.stoichiometry.electrons_per_carboxylation + +(forcing...) + sum(supply))

"""
    photosynthesize(photosynthesis, environment[, excitation])
    photosynthesize(leaf, environment[, excitation])

Evaluate leaf photosynthetic carbon assimilation and available electron supply
for `environment` using `BroadbandExcitation()` by default or an explicit
[`AbstractExcitation`](@ref).

Return the carbon fluxes of [`carbon_assimilation`](@ref) followed by
`electron_transport` [μmol e⁻ m⁻² leaf s⁻¹], the available linear electron supply
(use [`biochemical_turnover`](@ref) for realized `linear_electron_flow`);
[`CytochromePhotosynthesis`](@ref) appends `psi_electron_supply`, `psii_open_fraction`
and `psi_to_psii_ratio`. Inputs outside the physical domain (non-finite or negative
light, CO₂ or O₂, non-positive temperature, invalid excitation or model domain)
give `NaN` in every field.
"""
@inline function photosynthesize(
    p::AbstractC3Photosynthesis,
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    r, electrons, supply = light_supply(p, env, excitation)
    OutFT = leaf_result_type(p, r, supply, env.co2, env.oxygen)
    valid_co2, transport, extras = electron_state(p, r, supply, env.co2)
    carbon = carbon_assimilation(
        r, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation,
    )
    valid = valid_co2 && valid_leaf(p, r, env, excitation)
    return checked_result(OutFT, valid, (; carbon..., electron_transport = transport, extras...))
end

@inline function biochemical_turnover(
    p::AbstractC3Photosynthesis{L,SignedAssimilation},
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {L}
    r, electrons, supply = light_supply(p, env, excitation)
    OutFT = leaf_result_type(p, r, supply, env.co2, env.oxygen)
    valid_co2, transport, extras = electron_state(p, r, supply, env.co2)
    turnover = biochemical_turnover(
        r, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation,
    )
    valid = valid_co2 && valid_leaf(p, r, env, excitation)
    return checked_result(OutFT, valid, (; turnover..., electron_supply = transport, extras...))
end
