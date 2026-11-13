"""
    AbstractStomatalTarget

Supertype for instantaneous stomatal conductance target laws, kept separate from
prognostic time evolution.
"""
abstract type AbstractStomatalTarget end

"""
    AbstractStomatalModel

Supertype for equilibrium, prognostic, or disabled stomatal conductance models.
"""
abstract type AbstractStomatalModel end

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
`relative_humidity` (`∈ [0, 1]`), and leaf-surface `vpd` [Pa].
"""
Base.@kwdef struct StomatalDrivers{FT<:Real}
    assimilation::FT
    surface_co2::FT
    pressure::FT
    relative_humidity::FT
    vpd::FT
    function StomatalDrivers{FT}(
        assimilation,
        surface_co2,
        pressure,
        relative_humidity,
        vpd,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, StomatalDrivers{FT}, (assimilation, surface_co2, pressure, relative_humidity, vpd),
        )
        return new{FT}(assimilation, surface_co2, pressure, relative_humidity, vpd)
    end
    function StomatalDrivers(
        assimilation::FT,
        surface_co2::FT,
        pressure::FT,
        relative_humidity::FT,
        vpd::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, StomatalDrivers, (assimilation, surface_co2, pressure, relative_humidity, vpd),
        )
        return new{FT}(assimilation, surface_co2, pressure, relative_humidity, vpd)
    end
end
StomatalDrivers(a::Real, cs::Real, p::Real, rh::Real, vpd::Real) =
    StomatalDrivers(promote_compatible(StomatalDrivers, a, cs, p, rh, vpd)...)

"""
    BallBerry{FT}(; minimum_conductance, slope, micromole_to_mole)
    BallBerry(minimum_conductance, slope, micromole_to_mole)

Parameterize the Ball–Berry stomatal H₂O conductance target
``g_0 + g_1 \\max(A_n, 0) h_s / χ_s`` [mol H₂O m⁻² leaf s⁻¹], where
``χ_s = C_s / P`` is the leaf-surface CO₂ mole fraction (`surface_co2 / pressure`
in [`StomatalDrivers`](@ref)), with `minimum_conductance` (`g₀` [mol H₂O m⁻² leaf s⁻¹]),
dimensionless `slope` (`g₁`), and `micromole_to_mole` [mol μmol⁻¹] converting `A_n`.
"""
Base.@kwdef struct BallBerry{FT<:Real} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    micromole_to_mole::FT
    function BallBerry{FT}(minimum_conductance, slope, micromole_to_mole) where {FT<:Real}
        require_continuous_scalar(
            FT, BallBerry{FT}, (minimum_conductance, slope, micromole_to_mole),
        )
        return new{FT}(minimum_conductance, slope, micromole_to_mole)
    end
    function BallBerry(
        minimum_conductance::FT,
        slope::FT,
        micromole_to_mole::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, BallBerry, (minimum_conductance, slope, micromole_to_mole),
        )
        return new{FT}(minimum_conductance, slope, micromole_to_mole)
    end
end
BallBerry(g0::Real, g1::Real, conv::Real) =
    BallBerry(promote_compatible(BallBerry, g0, g1, conv)...)

"""
    Medlyn{FT}(; minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)
    Medlyn(minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)

Parameterize the Medlyn stomatal H₂O conductance target
``g_0 + r (1 + g_1 / \\sqrt{\\max(D_s, D_{\\min})}) \\max(A_n, 0) / χ_s``
[mol H₂O m⁻² leaf s⁻¹], where ``χ_s = C_s / P`` is the leaf-surface CO₂ mole fraction
(`surface_co2 / pressure` in [`StomatalDrivers`](@ref)), with `minimum_conductance`
(`g₀`), `slope` (`g₁` [Pa¹ᐟ²]), H₂O/CO₂ `diffusivity_ratio` (`r`),
`micromole_to_mole` [mol μmol⁻¹] converting `A_n`, and positive regularization
`vpd_floor` (`D_min` [Pa]).
"""
Base.@kwdef struct Medlyn{FT<:Real} <: AbstractStomatalTarget
    minimum_conductance::FT
    slope::FT
    diffusivity_ratio::FT
    micromole_to_mole::FT
    vpd_floor::FT
    function Medlyn{FT}(
        minimum_conductance,
        slope,
        diffusivity_ratio,
        micromole_to_mole,
        vpd_floor,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            Medlyn{FT},
            (minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor),
        )
        return new{FT}(minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)
    end
    function Medlyn(
        minimum_conductance::FT,
        slope::FT,
        diffusivity_ratio::FT,
        micromole_to_mole::FT,
        vpd_floor::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            Medlyn,
            (minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor),
        )
        return new{FT}(minimum_conductance, slope, diffusivity_ratio, micromole_to_mole, vpd_floor)
    end
end
Medlyn(g0::Real, g1::Real, r::Real, conv::Real, floor::Real) =
    Medlyn(promote_compatible(Medlyn, g0, g1, r, conv, floor)...)

"""
    SteadyStateStomata(target)

Wrap an [`AbstractStomatalTarget`](@ref) for steady-state stomatal target
evaluation and coupled [`equilibrium_gas_exchange`](@ref) solves.
"""
struct SteadyStateStomata{T<:AbstractStomatalTarget} <: AbstractStomatalModel
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
struct PrognosticStomata{T<:AbstractStomatalTarget,FT<:Real} <: AbstractStomatalModel
    target::T
    opening_timescale::FT
    closing_timescale::FT
    function PrognosticStomata(
        target::T,
        opening_timescale::FT,
        closing_timescale::FT,
    ) where {T<:AbstractStomatalTarget,FT<:Real}
        require_continuous_scalar(
            FT, PrognosticStomata, (target, opening_timescale, closing_timescale),
        )
        all(x -> isfinite(x) && x > zero(FT), (opening_timescale, closing_timescale)) ||
            throw(ArgumentError("stomatal timescales must be positive and finite"))
        return new{T,FT}(target, opening_timescale, closing_timescale)
    end
end
PrognosticStomata(target::AbstractStomatalTarget, opening_timescale::Real, closing_timescale::Real) =
    PrognosticStomata(
        target,
        promote_compatible(PrognosticStomata, opening_timescale, closing_timescale)...,
    )
PrognosticStomata(; target, opening_timescale, closing_timescale) =
    PrognosticStomata(target, opening_timescale, closing_timescale)

"""
    StomatalState(conductance)
    StomatalState{FT}(conductance)

Hold leaf stomatal H₂O `conductance` [mol H₂O m⁻² leaf s⁻¹] separately from
immutable model parameters.
"""
struct StomatalState{FT<:Real}
    conductance::FT
    function StomatalState{FT}(conductance) where {FT<:Real}
        require_continuous_scalar(FT, StomatalState{FT}, (conductance,))
        return new{FT}(conductance)
    end
    function StomatalState(conductance::FT) where {FT<:Real}
        require_continuous_scalar(FT, StomatalState, (conductance,))
        return new{FT}(conductance)
    end
end

"""
    target_conductance(target_or_model, drivers::StomatalDrivers)

Evaluate the steady-state stomatal H₂O conductance target [mol H₂O m⁻² leaf s⁻¹]
for `drivers`, or return `nothing` for [`NoStomatalDynamics`](@ref).
"""
@inline function target_conductance(p::BallBerry, d::StomatalDrivers)
    pos_assim = d.assimilation > zero(d.assimilation) ? d.assimilation : zero(d.assimilation)
    uptake = pos_assim * p.micromole_to_mole
    S = typeof(p.minimum_conductance + uptake * d.relative_humidity * d.pressure / d.surface_co2)
    (exact_zero_scalar(uptake) || (iszero(primal_scalar(uptake)) && iszero(primal_scalar(d.surface_co2)))) &&
        return S(p.minimum_conductance)
    return S(p.minimum_conductance + p.slope * uptake * d.relative_humidity * d.pressure / d.surface_co2)
end

@inline function target_conductance(p::Medlyn, d::StomatalDrivers)
    pos_assim = d.assimilation > zero(d.assimilation) ? d.assimilation : zero(d.assimilation)
    uptake = pos_assim * p.micromole_to_mole
    S = typeof(p.minimum_conductance + p.slope + uptake * d.pressure / d.surface_co2 + d.vpd)
    (exact_zero_scalar(uptake) || (iszero(primal_scalar(uptake)) && iszero(primal_scalar(d.surface_co2)))) &&
        return S(p.minimum_conductance)
    sensitivity = one(d.vpd) + p.slope / sqrt(max(d.vpd, p.vpd_floor))
    return S(p.minimum_conductance + p.diffusivity_ratio * sensitivity * uptake * d.pressure / d.surface_co2)
end

@inline target_conductance(p::AbstractStomatalModel, d::StomatalDrivers) =
    target_conductance(p.target, d)
@inline target_conductance(::NoStomatalDynamics, ::StomatalDrivers) = nothing

"""
    stomatal_tendency(model, state::StomatalState, drivers::StomatalDrivers)

Compute the stomatal H₂O conductance rate of change ``\\mathrm{d}g_{sw}/\\mathrm{d}t``
[mol H₂O m⁻² leaf s⁻²] without mutating `state` (`0` for [`NoStomatalDynamics`](@ref)).
"""
@inline function stomatal_tendency(p::PrognosticStomata, state::StomatalState, d::StomatalDrivers)
    target = target_conductance(p, d)
    tau = target > state.conductance ? p.opening_timescale : p.closing_timescale
    return (target - state.conductance) / tau
end

@inline stomatal_tendency(::NoStomatalDynamics, state::StomatalState, ::StomatalDrivers) =
    zero(state.conductance)

@inline stomatal_outputs(p::PrognosticStomata, state, drivers) =
    (;
        conductance_target = target_conductance(p, drivers),
        conductance_tendency = stomatal_tendency(p, state, drivers),
    )
@inline function stomatal_outputs(::NoStomatalDynamics, state, drivers)
    S = typeof(state.conductance + zero(drivers.assimilation))
    return (; conductance_target = S(state.conductance), conductance_tendency = zero(S))
end
