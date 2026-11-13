"""
    AbstractCompensationPoint

Supertype for photorespiratory CO₂ compensation point (`Γ*` [Pa]) parameterizations.
"""
abstract type AbstractCompensationPoint end

"""
    PrescribedCompensationPoint(response)

Wrap a temperature-only `Γ*` law `response` [Pa] calibrated at a fixed O₂ level.

The caller must supply an oxygen partial pressure compatible with `response`'s
calibration; this component does not infer oxygen scaling and is distinct from
the signed or clipped carbon-assimilation policy.
"""
struct PrescribedCompensationPoint{R} <: AbstractCompensationPoint
    response::R
end

"""
    OxygenScaledCompensationPoint(response, reference_oxygen)

Scale a reference `Γ*` temperature law linearly with biochemical-site oxygen:
``Γ^*(T, O) = Γ^*(T, O_{\\mathrm{ref}}) \\cdot O / O_{\\mathrm{ref}}``.

`response` evaluates `Γ*(T, O_ref)` [Pa] at the finite positive biochemical O₂
partial pressure `reference_oxygen` [Pa]. No ambient O₂ default or stoichiometric
recalibration is introduced, and evaluation without an explicit `oxygen` argument
is unsupported. The scalar law returns exact zero at `oxygen == 0`.
"""
struct OxygenScaledCompensationPoint{R,FT<:Real} <: AbstractCompensationPoint
    response::R
    reference_oxygen::FT
    function OxygenScaledCompensationPoint(
        response::R,
        reference_oxygen::FT,
    ) where {R,FT<:Real}
        require_continuous_scalar(FT, OxygenScaledCompensationPoint, (response, reference_oxygen))
        isfinite(reference_oxygen) && reference_oxygen > zero(FT) ||
            throw(ArgumentError("reference oxygen must be finite and positive [Pa]"))
        return new{R,FT}(response, reference_oxygen)
    end
end

"""
    compensation_point(response, temperature, oxygen)

Evaluate the photorespiratory CO₂ compensation point `Γ*` [Pa] at leaf
`temperature` [K] and biochemical-site `oxygen` partial pressure [Pa].

Plain temperature laws delegate to [`temperature_response`](@ref); implement a
method on `compensation_point` when adding a new [`AbstractCompensationPoint`](@ref).
"""
@inline compensation_point(response, temperature, oxygen) =
    temperature_response(response, temperature)
@inline compensation_point(p::PrescribedCompensationPoint, temperature, oxygen) =
    temperature_response(p.response, temperature)
@inline compensation_point(p::OxygenScaledCompensationPoint, temperature, oxygen) =
    temperature_response(p.response, temperature) * (oxygen / p.reference_oxygen)

# Only the prescribed component has an oxygen-free evaluation.
@inline temperature_response(p::PrescribedCompensationPoint, temperature) =
    temperature_response(p.response, temperature)
