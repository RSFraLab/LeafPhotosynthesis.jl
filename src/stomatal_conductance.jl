"""
    AbstractStomatalTarget

Supertype for instantaneous stomatal conductance target laws, kept separate from
prognostic time evolution.
"""
abstract type AbstractStomatalTarget <: AbstractLeafScalar end

"""
    AbstractStomatalModel

Supertype for equilibrium, prognostic, or disabled stomatal conductance models.
"""
abstract type AbstractStomatalModel <: AbstractLeafScalar end

"""
    NoStomatalDynamics()

Disable stomatal conductance evolution and retain the caller's prescribed conductance.
"""
struct NoStomatalDynamics <: AbstractStomatalModel end

"""
    StomatalDrivers{FT}(; assimilation, surface_co2, pressure, relative_humidity, vpd)
    StomatalDrivers(assimilation, surface_co2, pressure, relative_humidity, vpd)

Bundle leaf-surface drivers for stomatal target evaluation: net `assimilation`
[μmol CO₂ m⁻² leaf s⁻¹], `surface_co2` [Pa], total `pressure` [Pa],
`relative_humidity` (`∈ [0, 1]`), and leaf-surface `vpd` [Pa]. Each target law reads one
humidity driver ([`BallBerry`](@ref) `relative_humidity`, [`Medlyn`](@ref) `vpd`); the other
is not used and may be `NaN`.
"""
Base.@kwdef struct StomatalDrivers{FT<:Real} <: AbstractLeafScalar
    assimilation::FT
    surface_co2::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
    StomatalDrivers{FT}(fields::Vararg{Any,5}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, StomatalDrivers{FT}, fields)...)
end
StomatalDrivers(fields::Vararg{Real,5}) =
    StomatalDrivers{scalar_type(StomatalDrivers, fields...)}(fields...)

"""
    BallBerry{FT}(; minimum_conductance, slope, micromole_to_mole)
    BallBerry(minimum_conductance, slope, micromole_to_mole)

Parameterize the Ball–Berry stomatal H₂O conductance target
``g_0 + g_1 \\max(A_n, 0) h_s / χ_s`` [mol H₂O m⁻² leaf s⁻¹], where
``χ_s = C_s / P`` is the leaf-surface CO₂ mole fraction (`surface_co2 / pressure`
in [`StomatalDrivers`](@ref)), with `minimum_conductance` (`g₀` [mol H₂O m⁻² leaf s⁻¹]),
dimensionless `slope` (`g₁`), and `micromole_to_mole` [mol μmol⁻¹] converting `A_n`.

Throws `ArgumentError` unless `minimum_conductance ≥ 0`, `slope ≥ 0` and
`micromole_to_mole > 0` are finite; `minimum_conductance = 0` is valid (see
[`equilibrium_gas_exchange`](@ref) for the closed-stomata limit it allows). These are fixed
parameters, set once on the host, so an invalid value is reported where it is entered
rather than as `NaN` results later.
"""
Base.@kwdef struct BallBerry{FT<:Real} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    micromole_to_mole::FT
    function BallBerry{FT}(fields::Vararg{Any,3}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, BallBerry{FT}, fields)...)
        is_nonnegative(p.minimum_conductance) && is_nonnegative(p.slope) &&
            is_positive(p.micromole_to_mole) ||
            throw(ArgumentError("BallBerry requires finite minimum_conductance ≥ 0, slope ≥ 0 and micromole_to_mole > 0"))
        return p
    end
end
BallBerry(fields::Vararg{Real,3}) =
    BallBerry{scalar_type(BallBerry, fields...)}(fields...)

"""
    Medlyn{FT}(; minimum_conductance, slope, micromole_to_mole, vpd_floor)
    Medlyn(minimum_conductance, slope, micromole_to_mole, vpd_floor)

Parameterize the Medlyn stomatal H₂O conductance target
``g_0 + r (1 + g_1 / \\sqrt{\\max(D_s, D_{\\min})}) \\max(A_n, 0) / χ_s``
[mol H₂O m⁻² leaf s⁻¹], where ``χ_s = C_s / P`` is the leaf-surface CO₂ mole fraction
(`surface_co2 / pressure` in [`StomatalDrivers`](@ref)), with `minimum_conductance`
(`g₀`), `slope` (`g₁` [Pa¹ᐟ²]), `micromole_to_mole` [mol μmol⁻¹] converting `A_n`, and
positive regularization `vpd_floor` (`D_min` [Pa]).

The H₂O/CO₂ diffusivity ratio ``r`` is not a parameter of the target: it is the
`diffusivity_ratio` of the leaf's [`SeriesDiffusion`](@ref), passed to
[`target_conductance`](@ref) and read by [`equilibrium_gas_exchange`](@ref) and
[`leaf_tendencies`](@ref), so the stomatal law and the CO₂ diffusion path use one value.

Throws `ArgumentError` unless `minimum_conductance ≥ 0`, `slope ≥ 0`,
`micromole_to_mole > 0` and `vpd_floor > 0` are finite. These are fixed parameters, set
once on the host, so an invalid value is reported where it is entered rather than as
`NaN` results later.
"""
Base.@kwdef struct Medlyn{FT<:Real} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    micromole_to_mole::FT
    vpd_floor::FT
    function Medlyn{FT}(fields::Vararg{Any,4}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, Medlyn{FT}, fields)...)
        is_nonnegative(p.minimum_conductance) && is_nonnegative(p.slope) &&
            is_positive(p.micromole_to_mole) && is_positive(p.vpd_floor) ||
            throw(ArgumentError("Medlyn requires finite minimum_conductance ≥ 0, slope ≥ 0, and positive micromole_to_mole and vpd_floor"))
        return p
    end
end
Medlyn(fields::Vararg{Real,4}) =
    Medlyn{scalar_type(Medlyn, fields...)}(fields...)

"""
    SteadyStateStomata(target)
    SteadyStateStomata(; target)

Wrap an [`AbstractStomatalTarget`](@ref) for steady-state stomatal target
evaluation and coupled [`equilibrium_gas_exchange`](@ref) solves.
"""
Base.@kwdef struct SteadyStateStomata{T<:AbstractStomatalTarget} <: AbstractStomatalModel
    target::T
end

"""
    PrognosticStomata(; target, opening_timescale, closing_timescale)
    PrognosticStomata(target, opening_timescale, closing_timescale)

Parameterize first-order prognostic stomatal conductance relaxation
(Liu et al. 2024, Eq. 3) with independent positive finite `opening_timescale`
and `closing_timescale` [s].

Parameters are immutable; the caller owns [`StomatalState`](@ref) and time stepping.
"""
Base.@kwdef struct PrognosticStomata{T<:AbstractStomatalTarget,FT<:Real} <: AbstractStomatalModel
    target::T
    opening_timescale::FT
    closing_timescale::FT
    function PrognosticStomata{T,FT}(fields::Vararg{Any,3}) where {T<:AbstractStomatalTarget,FT<:Real}
        p = new{T,FT}(continuous_fields(FT, PrognosticStomata{T,FT}, fields)...)
        is_positive(p.opening_timescale) && is_positive(p.closing_timescale) ||
            throw(ArgumentError("stomatal timescales must be positive and finite"))
        return p
    end
end
PrognosticStomata(target::AbstractStomatalTarget, opening_timescale::Real, closing_timescale::Real) =
    PrognosticStomata{
        typeof(target),scalar_type(PrognosticStomata, opening_timescale, closing_timescale),
    }(target, opening_timescale, closing_timescale)

# Domain of the stomatal drivers for a target law, shared by the target, the tendency,
# `leaf_tendencies` and the equilibrium solver: finite assimilation of either sign, finite
# nonnegative surface CO₂ that is positive whenever assimilation is (uptake from zero CO₂
# has no steady state), positive finite pressure, and the humidity driver of the law
# (`valid_humidity`).
@inline valid_stomatal_drivers(p::AbstractStomatalTarget, d::StomatalDrivers) =
    is_finite(d.assimilation) && is_nonnegative(d.surface_co2) &&
    (is_positive(d.surface_co2) || !exceeds(d.assimilation, zero(d.assimilation))) &&
    is_positive(d.pressure) && valid_humidity(p, d.relative_humidity, d.vpd)

# Each law checks only the humidity driver it reads: Ball–Berry the relative humidity
# (in [0, 1]), Medlyn the vapour pressure deficit (finite; values below `vpd_floor`,
# including negative ones, are floored). The other driver is ignored and may be `NaN`.
@inline valid_humidity(::BallBerry, relative_humidity, vpd) = is_fraction(relative_humidity)
@inline valid_humidity(::Medlyn, relative_humidity, vpd) = is_finite(vpd)

"""
    target_conductance(target, drivers::StomatalDrivers, diffusion::SeriesDiffusion)
    target_conductance(stomatal_model, drivers::StomatalDrivers, diffusion::SeriesDiffusion)
    target_conductance(leaf, drivers::StomatalDrivers, diffusion::SeriesDiffusion)

Evaluate the steady-state stomatal H₂O conductance target [mol H₂O m⁻² leaf s⁻¹] of a
[`BallBerry`](@ref) or [`Medlyn`](@ref) target (directly, in a [`SteadyStateStomata`](@ref)
or [`PrognosticStomata`](@ref), or in a [`LeafModel`](@ref)) for `drivers`, or return
`nothing` for [`NoStomatalDynamics`](@ref).

`diffusion` is the leaf's CO₂ diffusion path; its `diffusivity_ratio` is the H₂O/CO₂ ratio
``r`` of the Medlyn law (the one value the stomatal law and the diffusion share). Ball–Berry
does not use it (its slope already includes the ratio).

Drivers outside their domain give `NaN`, never the minimum conductance: non-finite
assimilation, negative or non-finite surface CO₂ (zero only without net uptake), non-positive
or non-finite pressure, and an invalid humidity driver of the law: relative humidity outside
`[0, 1]` for Ball–Berry, a non-finite VPD for Medlyn. Each law checks only the driver it
uses; the other may be `NaN`. Net assimilation ≤ 0 gives the minimum conductance `g₀`.
"""
@inline function target_conductance(p::BallBerry, d::StomatalDrivers, diffusion::SeriesDiffusion)
    require_compatible(target_conductance, p, d, diffusion)
    pos_assim = exceeds(d.assimilation, zero(d.assimilation)) ? d.assimilation : zero(d.assimilation)
    uptake = pos_assim * p.micromole_to_mole
    S = typeof(p.minimum_conductance + p.slope * uptake * d.relative_humidity * d.pressure / d.surface_co2)
    valid_stomatal_drivers(p, d) || return S(NaN)
    (exact_zero_scalar(uptake) || (iszero(primal_scalar(uptake)) && iszero(primal_scalar(d.surface_co2)))) &&
        return S(p.minimum_conductance)
    return S(p.minimum_conductance + p.slope * uptake * d.relative_humidity * d.pressure / d.surface_co2)
end

@inline function target_conductance(p::Medlyn, d::StomatalDrivers, diffusion::SeriesDiffusion)
    require_compatible(target_conductance, p, d, diffusion)
    ratio = diffusion.diffusivity_ratio
    pos_assim = exceeds(d.assimilation, zero(d.assimilation)) ? d.assimilation : zero(d.assimilation)
    uptake = pos_assim * p.micromole_to_mole
    S = typeof(p.minimum_conductance + p.slope + ratio + p.vpd_floor +
        uptake * d.pressure / d.surface_co2 + d.vpd)
    valid_stomatal_drivers(p, d) && is_positive(ratio) || return S(NaN)
    (exact_zero_scalar(uptake) || (iszero(primal_scalar(uptake)) && iszero(primal_scalar(d.surface_co2)))) &&
        return S(p.minimum_conductance)
    sensitivity = one(d.vpd) + p.slope / sqrt(primal_max(d.vpd, p.vpd_floor))
    return S(p.minimum_conductance + ratio * sensitivity * uptake * d.pressure / d.surface_co2)
end

@inline target_conductance(p::Union{SteadyStateStomata,PrognosticStomata}, d::StomatalDrivers, diffusion::SeriesDiffusion) =
    target_conductance(p.target, d, diffusion)
@inline target_conductance(::NoStomatalDynamics, ::StomatalDrivers, ::SeriesDiffusion) = nothing

"""
    stomatal_tendency(model_or_leaf, state::StomatalState, drivers::StomatalDrivers,
        diffusion::SeriesDiffusion)

Compute the stomatal H₂O conductance rate of change ``\\mathrm{d}g_{sw}/\\mathrm{d}t``
[mol H₂O m⁻² leaf s⁻²] without mutating `state` (`0` for [`NoStomatalDynamics`](@ref)).
`diffusion` supplies the diffusivity ratio of the target, as in [`target_conductance`](@ref).
Invalid drivers (see [`target_conductance`](@ref)) or a negative or non-finite state
conductance give `NaN`; a closed state (`conductance = 0`) is valid.
"""
@inline function stomatal_tendency(p::PrognosticStomata, state::StomatalState, d::StomatalDrivers, diffusion::SeriesDiffusion)
    require_compatible(stomatal_tendency, p, state, d, diffusion)
    target = target_conductance(p, d, diffusion)
    tau = exceeds(target, state.conductance) ? p.opening_timescale : p.closing_timescale
    tendency = (target - state.conductance) / tau
    return is_nonnegative(state.conductance) ? tendency : oftype(tendency, NaN)
end

@inline stomatal_tendency(::NoStomatalDynamics, state::StomatalState, ::StomatalDrivers, ::SeriesDiffusion) =
    zero(state.conductance)

@inline stomatal_outputs(p::PrognosticStomata, state, drivers, diffusion) =
    (;
        conductance_target = target_conductance(p, drivers, diffusion),
        conductance_tendency = stomatal_tendency(p, state, drivers, diffusion),
    )
@inline function stomatal_outputs(::NoStomatalDynamics, state, drivers, diffusion)
    S = typeof(state.conductance + zero(drivers.assimilation))
    return (; conductance_target = S(state.conductance), conductance_tendency = zero(S))
end
