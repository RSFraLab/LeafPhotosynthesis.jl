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
    CoreOxygenExchange{FT}(; electrons_per_oxygen, oxygen_per_oxygenation,
        co2_per_oxygenation, respiratory_quotient)
    CoreOxygenExchange(electrons_per_oxygen, oxygen_per_oxygenation,
        co2_per_oxygenation, respiratory_quotient)

Specify dimensionless stoichiometric coefficients for PSII water-splitting O₂
evolution, Rubisco photorespiratory O₂ uptake, and mitochondrial respiratory O₂
uptake (`respiratory_quotient` is CO₂ released per O₂ consumed). Throws
`ArgumentError` unless `electrons_per_oxygen`, `co2_per_oxygenation` and
`respiratory_quotient` are finite and positive and `oxygen_per_oxygenation` is finite and
nonnegative.

# Examples
```jldoctest
julia> using LeafPhotosynthesis: OxygenFluxInputs

julia> ox = CoreOxygenExchange(4.0, 1.5, 0.5, 1.0);

julia> res = oxygen_exchange(ox, OxygenFluxInputs(80.0, 4.0, 1.0));

julia> res.net_release
13.0
```
"""
Base.@kwdef struct CoreOxygenExchange{FT<:Real} <: AbstractOxygenExchange
    electrons_per_oxygen::FT
    oxygen_per_oxygenation::FT
    co2_per_oxygenation::FT
    respiratory_quotient::FT
    function CoreOxygenExchange{FT}(fields::Vararg{Any,4}) where {FT<:Real}
        p = new{FT}(continuous_fields(FT, CoreOxygenExchange{FT}, fields)...)
        is_positive(p.electrons_per_oxygen) && is_nonnegative(p.oxygen_per_oxygenation) &&
            is_positive(p.co2_per_oxygenation) && is_positive(p.respiratory_quotient) ||
            throw(ArgumentError("CoreOxygenExchange requires finite electrons_per_oxygen > 0, oxygen_per_oxygenation ≥ 0, co2_per_oxygenation > 0 and respiratory_quotient > 0"))
        return p
    end
end
CoreOxygenExchange(fields::Vararg{Real,4}) =
    CoreOxygenExchange{scalar_type(CoreOxygenExchange, fields...)}(fields...)

"""
    OxygenFluxInputs{FT}(; linear_electron_flow, rubisco_oxygenation, respiration)
    OxygenFluxInputs(linear_electron_flow, rubisco_oxygenation, respiration)

Bundle realized `linear_electron_flow` [μmol electrons m⁻² leaf s⁻¹],
`rubisco_oxygenation` [μmol reactions m⁻² leaf s⁻¹], and mitochondrial CO₂
`respiration` [μmol CO₂ m⁻² leaf s⁻¹] for [`oxygen_exchange`](@ref).
"""
Base.@kwdef struct OxygenFluxInputs{FT<:Real} <: AbstractLeafScalar
    linear_electron_flow::FT
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

With [`OxygenFluxInputs`](@ref), gross O₂ evolution is the given linear electron flow divided
by `electrons_per_oxygen`. This agrees with the coupled form below only when the linear flow
is `electrons_per_oxygen` times `v_c + (o_po - c_po) v_o` (for example the classical FvCB
stoichiometry with `e_c = 4`, `a = 2` and four electrons per O₂); `electrons_per_oxygen` is
used only by this form.

When coupled to a signed [`C3Photosynthesis`](@ref) or
[`CytochromePhotosynthesis`](@ref) model, the budget is computed from the realized
carboxylation ``v_c`` and oxygenation ``v_o`` (= photorespiratory CO₂ / `co2_per_oxygenation`),
evaluated prior to net carbon cancellation so the compensation point causes no `0/0`
singularity. Gross O₂ evolution is the O₂ whose reductant the carbon cycles consume:
``v_c + (o_{po} - c_{po}) v_o``, one O₂ per CO₂ reduced to carbohydrate and, per
oxygenation, the O₂ taken up (`oxygen_per_oxygenation`, ``o_{po}``) less the CO₂ released
(`co2_per_oxygenation`, ``c_{po}``) by the photorespiratory cycle; with the canonical
1.5 and 0.5 this is ``v_c + v_o`` (two NADPH per carboxylation and per oxygenation,
Johnson & Berry 2021, Eqs. 7–9). It does not depend on the electron stoichiometry of the
photosynthesis model: electrons beyond the NADPH demand (an ATP-limited FvCB stoichiometry
such as ``e_c = 4.5``) are taken to return to O₂ without net O₂ exchange. Net O₂ release is
therefore ``v_c - c_{po} v_o - R_d / RQ``, equal to net CO₂ uptake when `respiratory_quotient`
is one. The realized `linear_electron_flow` is still reported.
"""
@inline function oxygen_exchange(p::CoreOxygenExchange, x::OxygenFluxInputs)
    require_compatible(oxygen_exchange, p, x)
    gross = x.linear_electron_flow / p.electrons_per_oxygen
    photorespiration = p.oxygen_per_oxygenation * x.rubisco_oxygenation
    respiration = x.respiration / p.respiratory_quotient
    net = gross - photorespiration - respiration
    return typed_result(typeof(gross + photorespiration + respiration), (;
        gross_production = gross,
        photorespiratory_uptake = photorespiration,
        respiratory_uptake = respiration,
        net_release = net,
    ))
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
    gross = turnover.carboxylation +
        (p.oxygen_per_oxygenation - p.co2_per_oxygenation) * oxygenation
    photorespiration = p.oxygen_per_oxygenation * oxygenation
    respiration = turnover.respiration / p.respiratory_quotient
    result = (;
        gross_production = gross,
        photorespiratory_uptake = photorespiration,
        respiratory_uptake = respiration,
        net_release = gross - photorespiration - respiration,
        turnover.carboxylation,
        rubisco_oxygenation = oxygenation,
        turnover.linear_electron_flow,
        turnover_electron_flows(model, turnover)...,
    )
    return typed_result(result_type(result), result)
end
