"""Interface for the photorespiratory CO₂ compensation point Γ* [Pa]."""
abstract type AbstractCompensationPoint end

"""
    PrescribedCompensationPoint(response)

Explicit temperature-only Γ* response. The caller must supply an oxygen level
compatible with its calibration. This component does not infer O₂ dependence.
It is distinct from the signed/clipped carbon-assimilation compensation policy.
"""
struct PrescribedCompensationPoint{R} <: AbstractCompensationPoint
    response::R
end

"""
    OxygenScaledCompensationPoint(response, reference_oxygen)

Γ*(T,O) = Γ*(T,Oref) O/Oref at unchanged Rubisco specificity and photorespiratory
stoichiometry. `response` is the entire Γ* temperature law [Pa] calibrated at the
explicit biochemical O₂ partial pressure `reference_oxygen` [Pa], finite and >0.
No ambient O₂ default, total-pressure multiplier or new stoichiometric constant
is introduced. All inputs/parameters use one FT. Evaluation without an oxygen
argument is unsupported.

The scalar law returns exact zero at O=0. Full gas exchange with this scaled
treatment currently requires positive Γ*/O₂; this law alone does not model anaerobic respiration or zero-gas
limits. Environmental oxygen must be finite and nonnegative.
"""
struct OxygenScaledCompensationPoint{R,FT<:AbstractFloat} <: AbstractCompensationPoint
    response::R
    reference_oxygen::FT
    function OxygenScaledCompensationPoint(
        response::R,
        reference_oxygen::FT,
    ) where {R,FT<:AbstractFloat}
        isfinite(reference_oxygen) && reference_oxygen > zero(FT) ||
            throw(ArgumentError("reference oxygen must be finite and positive [Pa]"))
        new{R,FT}(response, reference_oxygen)
    end
end

"""
    compensation_point(response, temperature, oxygen)

Evaluate Γ* [Pa] at temperature [K] and biochemical O₂ partial pressure [Pa].
Legacy `temperature_response(response,T)` extensions remain prescribed-Γ* laws.
Implement this function for a new AbstractCompensationPoint component.
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
