"""
    AbstractCompensationPoint

Supertype for photorespiratory CO₂ compensation point (`Γ*` [Pa]) parameterizations.
"""
abstract type AbstractCompensationPoint <: AbstractLeafScalar end

# A Γ* field of the rate bundles holds a temperature law or a compensation-point component.
const GammaStarLaw = Union{AbstractTemperatureResponse,AbstractCompensationPoint}

"""
    PrescribedCompensationPoint(response)
    PrescribedCompensationPoint(; response)

Wrap a temperature-only `Γ*` law `response` [Pa] calibrated at a fixed O₂ level.

The caller must supply an oxygen partial pressure compatible with `response`'s
calibration; this component does not infer oxygen scaling and is distinct from
the signed or clipped carbon-assimilation policy.
"""
Base.@kwdef struct PrescribedCompensationPoint{R<:AbstractTemperatureResponse} <: AbstractCompensationPoint
    response::R
end

"""
    OxygenScaledCompensationPoint(response, reference_oxygen)
    OxygenScaledCompensationPoint(; response, reference_oxygen)

Scale a reference `Γ*` temperature law linearly with biochemical-site oxygen:
``Γ^*(T, O) = Γ^*(T, O_{\\mathrm{ref}}) \\cdot O / O_{\\mathrm{ref}}``.

`response` evaluates `Γ*(T, O_ref)` [Pa] at the finite positive biochemical O₂
partial pressure `reference_oxygen` [Pa]. No ambient O₂ default or stoichiometric
recalibration is introduced, and evaluation without an explicit `oxygen` argument
is unsupported. The scalar law returns exact zero at `oxygen == 0`.
"""
Base.@kwdef struct OxygenScaledCompensationPoint{R<:AbstractTemperatureResponse,FT<:Real} <:
       AbstractCompensationPoint
    response::R
    reference_oxygen::FT
    function OxygenScaledCompensationPoint{R,FT}(
        fields::Vararg{Any,2},
    ) where {R<:AbstractTemperatureResponse,FT<:Real}
        p = new{R,FT}(continuous_fields(FT, OxygenScaledCompensationPoint{R,FT}, fields)...)
        is_positive(p.reference_oxygen) ||
            throw(ArgumentError("reference oxygen must be finite and positive [Pa]"))
        return p
    end
end
OxygenScaledCompensationPoint(response::AbstractTemperatureResponse, reference_oxygen::Real) =
    OxygenScaledCompensationPoint{
        typeof(response),scalar_type(OxygenScaledCompensationPoint, reference_oxygen),
    }(response, reference_oxygen)

"""
    compensation_point(response, temperature, oxygen)

Evaluate the photorespiratory CO₂ compensation point `Γ*` [Pa] at leaf
`temperature` [K] and biochemical-site `oxygen` partial pressure [Pa].

Plain temperature laws delegate to [`temperature_response`](@ref); implement a
method on `compensation_point` when adding a new [`AbstractCompensationPoint`](@ref).
"""
@inline compensation_point(response, temperature, oxygen) =
    temperature_response(response, temperature)
@inline function compensation_point(p::PrescribedCompensationPoint, temperature, oxygen)
    require_compatible(compensation_point, p, temperature, oxygen)
    return temperature_response(p.response, temperature)
end
@inline function compensation_point(p::OxygenScaledCompensationPoint, temperature, oxygen)
    require_compatible(compensation_point, p, temperature, oxygen)
    return temperature_response(p.response, temperature) * (oxygen / p.reference_oxygen)
end

# Only the prescribed component has an oxygen-free evaluation.
@inline temperature_response(p::PrescribedCompensationPoint, temperature) =
    temperature_response(p.response, temperature)
