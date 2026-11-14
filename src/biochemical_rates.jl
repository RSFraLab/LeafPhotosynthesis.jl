"""
    C3TemperatureResponse(; vcmax, jmax, respiration, kc, ko, gamma_star)

Bundle the six C3 biochemical temperature responses (`vcmax`, `jmax`,
`respiration`, `kc`, `ko`, and `gamma_star`).

Every field is an [`AbstractTemperatureResponse`](@ref); `gamma_star` may also be an explicit
[`AbstractCompensationPoint`](@ref) component.
"""
Base.@kwdef struct C3TemperatureResponse{
    V<:AbstractTemperatureResponse,
    J<:AbstractTemperatureResponse,
    R<:AbstractTemperatureResponse,
    C<:AbstractTemperatureResponse,
    O<:AbstractTemperatureResponse,
    G<:GammaStarLaw,
} <: AbstractLeafScalar
    vcmax::V
    jmax::J
    respiration::R
    kc::C
    ko::O
    gamma_star::G
    C3TemperatureResponse{V,J,R,C,O,G}(fields::Vararg{Any,6}) where {V,J,R,C,O,G} =
        new{V,J,R,C,O,G}(component_fields(C3TemperatureResponse{V,J,R,C,O,G}, fields)...)
end
C3TemperatureResponse(
    vcmax::AbstractTemperatureResponse,
    jmax::AbstractTemperatureResponse,
    respiration::AbstractTemperatureResponse,
    kc::AbstractTemperatureResponse,
    ko::AbstractTemperatureResponse,
    gamma_star::GammaStarLaw,
) = C3TemperatureResponse{map(typeof, (vcmax, jmax, respiration, kc, ko, gamma_star))...}(
    vcmax, jmax, respiration, kc, ko, gamma_star,
)

"""
    C3Rates{FT}(; vcmax, jmax, respiration, kc, ko, gamma_star)
    C3Rates(vcmax, jmax, respiration, kc, ko, gamma_star)

Store evaluated C3 biochemical kinetic parameters at leaf temperature and, when
selected, biochemical-site O₂.

- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹.
- `jmax`: μmol electrons m⁻² leaf s⁻¹.
- `kc`, `ko`, `gamma_star`: Pa.
"""
Base.@kwdef struct C3Rates{FT<:Real} <: AbstractLeafScalar
    vcmax::FT
    jmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    C3Rates{FT}(fields::Vararg{Any,6}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, C3Rates{FT}, fields)...)
end
C3Rates(fields::Vararg{Real,6}) =
    C3Rates{scalar_type(C3Rates, fields...)}(fields...)

"""
    biochemical_rates(responses, temperature[, oxygen]) -> C3Rates or CytochromeRates

Evaluate the kinetic rate parameters of `responses` ([`C3TemperatureResponse`](@ref) or
[`CytochromeTemperatureResponse`](@ref)) at leaf `temperature` [K] and, when given,
biochemical-site `oxygen` partial pressure [Pa], returning a [`C3Rates`](@ref) or a
[`CytochromeRates`](@ref).

The three-argument form passes `oxygen` to [`compensation_point`](@ref) for `gamma_star`;
the two-argument form evaluates temperature-only laws, so it is a `MethodError` for an
[`OxygenScaledCompensationPoint`](@ref). A custom rate bundle needs only the two-argument
method: the three-argument call falls back to it.
"""
@inline function biochemical_rates(p::C3TemperatureResponse, temperature, oxygen = nothing)
    require_compatible(biochemical_rates, p, temperature, oxygen)
    return C3Rates(
        temperature_response(p.vcmax, temperature),
        temperature_response(p.jmax, temperature),
        temperature_response(p.respiration, temperature),
        temperature_response(p.kc, temperature),
        temperature_response(p.ko, temperature),
        gamma_star_response(p.gamma_star, temperature, oxygen),
    )
end

@inline biochemical_rates(responses, temperature, oxygen) = biochemical_rates(responses, temperature)

# Γ* of a rate bundle: temperature-only without oxygen, `compensation_point` with it.
@inline gamma_star_response(law, temperature, ::Nothing) = temperature_response(law, temperature)
@inline gamma_star_response(law, temperature, oxygen) = compensation_point(law, temperature, oxygen)
