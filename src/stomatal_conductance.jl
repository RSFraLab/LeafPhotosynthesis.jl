"""Interface for an instantaneous stomatal target law, independent of time evolution."""
abstract type AbstractStomatalTarget end
"""Interface for equilibrium or prognostic stomatal behavior."""
abstract type AbstractStomatalModel end

"""Disable stomatal evolution; use the caller's prescribed conductance unchanged."""
struct NoStomatalDynamics <: AbstractStomatalModel end

"""
    StomatalDrivers{FT}(; assimilation, surface_co2, pressure, relative_humidity, vpd)

Net assimilation [μmol CO₂ m⁻² leaf s⁻¹], leaf-surface CO₂ [Pa], total pressure
[Pa], leaf-surface relative humidity [0,1], and leaf-surface VPD [Pa]. Supply
positive CO₂/pressure. Humidity is an explicit boundary input, not solved here.
"""
Base.@kwdef struct StomatalDrivers{FT}
    assimilation::FT
    surface_co2::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
end

"""
Ball–Berry target: g₀ + g₁ max(A,0) hₛ / χₛ.
g₀: mol H₂O m⁻² leaf s⁻¹; slope dimensionless; conversion: mol per μmol.
Supply nonnegative g₀/slope and a positive conversion explicitly.
"""
Base.@kwdef struct BallBerry{FT} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    micromole_to_mole::FT
end

"""
Medlyn target: g₀ + r (1 + g₁ / √D) max(A,0) / χₛ.
g₀: mol H₂O m⁻² leaf s⁻¹; slope: Pa^½ (not kPa^½); r is the H₂O/CO₂
diffusivity ratio. `vpd_floor` [Pa] is an explicit positive regularization.
The conversion is mol per μmol. All coefficients are required.
"""
Base.@kwdef struct Medlyn{FT} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    diffusivity_ratio::FT
    micromole_to_mole::FT
    vpd_floor::FT
end

"""Instantaneous target evaluation; a fully coupled equilibrium gas-exchange solve is separate."""
struct SteadyStateStomata{T<:AbstractStomatalTarget} <: AbstractStomatalModel
    target::T
end

"""
    PrognosticStomata(; target, opening_timescale, closing_timescale)

First-order conductance response from Liu et al. (2024), Eq. (3), with independently
specified opening and closing timescales [s]. Parameters are immutable; the caller
owns conductance state and integration. Both timescales must be positive and finite.
"""
struct PrognosticStomata{T<:AbstractStomatalTarget,FT} <: AbstractStomatalModel
    target::T
    opening_timescale::FT
    closing_timescale::FT
end
function PrognosticStomata(; target, opening_timescale::FT, closing_timescale::FT) where {FT}
    all(x -> isfinite(x) && x > zero(FT), (opening_timescale, closing_timescale)) ||
        throw(ArgumentError("stomatal timescales must be positive and finite"))
    return PrognosticStomata(target, opening_timescale, closing_timescale)
end

"""Prognostic H₂O conductance [mol m⁻² leaf s⁻¹]; separate from model parameters."""
struct StomatalState{FT}
    conductance::FT
end

@inline function target_conductance(p::BallBerry, d::StomatalDrivers)
    uptake = max(d.assimilation, zero(d.assimilation)) * p.micromole_to_mole
    return p.minimum_conductance + p.slope * uptake * d.relative_humidity * d.pressure / d.surface_co2
end

@inline function target_conductance(p::Medlyn, d::StomatalDrivers)
    uptake = max(d.assimilation, zero(d.assimilation)) * p.micromole_to_mole
    sensitivity = one(d.vpd) + p.slope / sqrt(max(d.vpd, p.vpd_floor))
    return p.minimum_conductance + p.diffusivity_ratio * sensitivity * uptake * d.pressure / d.surface_co2
end

@inline target_conductance(p::AbstractStomatalModel, d::StomatalDrivers) = target_conductance(p.target, d)
@inline target_conductance(::NoStomatalDynamics, ::StomatalDrivers) = nothing

"""Return dg_sw/dt [mol H₂O m⁻² leaf s⁻²]; never advance or mutate the state."""
@inline function stomatal_tendency(p::PrognosticStomata, state::StomatalState, d::StomatalDrivers)
    target = target_conductance(p, d)
    tau = target > state.conductance ? p.opening_timescale : p.closing_timescale
    return (target - state.conductance) / tau
end

@inline stomatal_tendency(::NoStomatalDynamics, state::StomatalState, ::StomatalDrivers) = zero(state.conductance)

@inline stomatal_outputs(p::PrognosticStomata, state, drivers) =
    (; conductance_target = target_conductance(p, drivers), conductance_tendency = stomatal_tendency(p, state, drivers))
@inline stomatal_outputs(::NoStomatalDynamics, state, drivers) =
    (; conductance_target = state.conductance, conductance_tendency = zero(state.conductance))
