"""
    AbstractOxygenExchange

Supertype for leaf oxygen-exchange components; positive net flux denotes O₂ release.
"""
abstract type AbstractOxygenExchange <: AbstractLeafScalar end

"""
    NoOxygenExchange()

Disable oxygen-budget calculations explicitly; [`oxygen_exchange`](@ref) returns `nothing`.
"""
struct NoOxygenExchange <: AbstractOxygenExchange end
@inline oxygen_exchange(::NoOxygenExchange, args...) = nothing

"""
    CoreOxygenExchange{FT}(; oxygen_per_oxygenation, co2_per_oxygenation, respiratory_quotient)
    CoreOxygenExchange(oxygen_per_oxygenation, co2_per_oxygenation, respiratory_quotient)

Specify dimensionless stoichiometric coefficients of the leaf O₂ budget: O₂ taken up
(`oxygen_per_oxygenation`) and CO₂ released (`co2_per_oxygenation`) by the photorespiratory
cycle per Rubisco oxygenation, and the mitochondrial `respiratory_quotient` (CO₂ released per
O₂ consumed). Throws `ArgumentError` unless `co2_per_oxygenation` and `respiratory_quotient`
are finite and positive and `oxygen_per_oxygenation` is finite and nonnegative.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: OxygenFluxInputs

julia> ox = CoreOxygenExchange(1.5, 0.5, 1.0);

julia> res = oxygen_exchange(ox, OxygenFluxInputs(16.0, 4.0, 1.0));

julia> res.net_release
13.0
```
"""
Base.@kwdef struct CoreOxygenExchange{FT<:Real} <: AbstractOxygenExchange
    oxygen_per_oxygenation::FT
    co2_per_oxygenation::FT
    respiratory_quotient::FT
    function CoreOxygenExchange{FT}(fields::Vararg{Any,3}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, CoreOxygenExchange{FT}, fields)...)
        is_nonnegative(p.oxygen_per_oxygenation) && is_positive(p.co2_per_oxygenation) &&
            is_positive(p.respiratory_quotient) ||
            throw(ArgumentError("CoreOxygenExchange requires finite oxygen_per_oxygenation ≥ 0, co2_per_oxygenation > 0 and respiratory_quotient > 0"))
        return p
    end
end
CoreOxygenExchange(fields::Vararg{Real,3}) =
    CoreOxygenExchange{scalar_type(CoreOxygenExchange, fields...)}(fields...)

"""
    OxygenFluxInputs{FT}(; carboxylation, rubisco_oxygenation, respiration)
    OxygenFluxInputs(carboxylation, rubisco_oxygenation, respiration)

Bundle realized Rubisco `carboxylation` and `rubisco_oxygenation` [μmol reactions m⁻² leaf
s⁻¹] and mitochondrial CO₂ `respiration` [μmol CO₂ m⁻² leaf s⁻¹] for
[`oxygen_exchange`](@ref), the turnover the coupled form takes from
[`biochemical_turnover`](@ref).
"""
Base.@kwdef struct OxygenFluxInputs{FT<:Real} <: AbstractLeafScalar
    carboxylation::FT
    rubisco_oxygenation::FT
    respiration::FT
    OxygenFluxInputs{FT}(fields::Vararg{Any,3}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, OxygenFluxInputs{FT}, fields)...)
end
OxygenFluxInputs(fields::Vararg{Real,3}) =
    OxygenFluxInputs{scalar_type(OxygenFluxInputs, fields...)}(fields...)

"""
    oxygen_exchange(oxygen_model, inputs::OxygenFluxInputs)
    oxygen_exchange(oxygen_model, photosynthesis, environment,
        excitation=BroadbandExcitation())
    oxygen_exchange(leaf, environment, excitation=BroadbandExcitation())

Return `(; gross_production, photorespiratory_uptake, respiratory_uptake, net_release)`
[μmol O₂ m⁻² leaf s⁻¹]: gross O₂ evolution, photorespiratory and mitochondrial respiratory
O₂ uptake, and net O₂ release (`nothing` for [`NoOxygenExchange`](@ref)). The coupled forms
append the realized `carboxylation` and `rubisco_oxygenation` [μmol m⁻² leaf s⁻¹] and the
realized `linear_electron_flow` [μmol e⁻ m⁻² leaf s⁻¹]; [`CytochromePhotosynthesis`](@ref)
adds `psi_electron_flow` and `cyclic_electron_flow`. Input outside the physical domain gives
`NaN` in every field.

Both forms compute the budget from carboxylation ``v_c`` and oxygenation ``v_o``: with
[`OxygenFluxInputs`](@ref) the caller gives them; when coupled to a signed
[`C3Photosynthesis`](@ref) or [`CytochromePhotosynthesis`](@ref) model they are the realized
turnover (``v_o`` = photorespiratory CO₂ / `co2_per_oxygenation`), evaluated prior to net
carbon cancellation so the compensation point causes no `0/0` singularity. The same fluxes
give the same budget, bit for bit, in both forms. Gross O₂ evolution is the O₂ whose
reductant the carbon cycles consume: ``v_c + (o_{po} - c_{po}) v_o``, one O₂ per CO₂ reduced
to carbohydrate and, per oxygenation, the O₂ taken up (`oxygen_per_oxygenation`, ``o_{po}``)
less the CO₂ released (`co2_per_oxygenation`, ``c_{po}``) by the photorespiratory cycle;
with the canonical 1.5 and 0.5 this is ``v_c + v_o`` (two NADPH per carboxylation and per
oxygenation, Johnson & Berry 2021, Eqs. 7–9). It does not depend on the electron
stoichiometry of the photosynthesis model: electrons beyond the NADPH demand (an
ATP-limited FvCB stoichiometry such as ``e_c = 4.5``) are taken to return to O₂ without net
O₂ exchange. Net O₂ release is therefore ``v_c - c_{po} v_o - R_d / RQ``, equal to net CO₂
uptake when `respiratory_quotient` is one. The coupled forms still report the realized
`linear_electron_flow`.
"""
@inline function oxygen_exchange(p::CoreOxygenExchange, x::OxygenFluxInputs)
    require_compatible(oxygen_exchange, p, x)
    budget = oxygen_budget(p, x.carboxylation, x.rubisco_oxygenation, x.respiration)
    return typed_result(result_type(budget), budget)
end

# The O₂ budget of both `oxygen_exchange` forms, from carboxylation, oxygenation and
# mitochondrial CO₂ respiration.
@inline function oxygen_budget(p::CoreOxygenExchange, carboxylation, oxygenation, respiration_co2)
    gross = carboxylation + (p.oxygen_per_oxygenation - p.co2_per_oxygenation) * oxygenation
    photorespiration = p.oxygen_per_oxygenation * oxygenation
    respiration = respiration_co2 / p.respiratory_quotient
    return (;
        gross_production = gross,
        photorespiratory_uptake = photorespiration,
        respiratory_uptake = respiration,
        net_release = gross - photorespiration - respiration,
    )
end

@inline function oxygen_exchange(
    p::CoreOxygenExchange,
    model::AbstractC3Photosynthesis{L,SignedAssimilation},
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {L}
    require_compatible(oxygen_exchange, p, model, env, excitation)
    turnover = biochemical_turnover(model, env, excitation)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    result = (;
        oxygen_budget(p, turnover.carboxylation, oxygenation, turnover.respiration)...,
        turnover.carboxylation,
        rubisco_oxygenation = oxygenation,
        turnover.linear_electron_flow,
        turnover_electron_flows(model, turnover)...,
    )
    return typed_result(result_type(result), result)
end
