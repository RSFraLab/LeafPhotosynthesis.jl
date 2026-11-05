"""
    Arrhenius{FT}(; value_ref, activation_energy, temperature_ref, gas_constant)

Temperature response normalized at `temperature_ref` [K]. `activation_energy`
is in J mol⁻¹ and `gas_constant` in J mol⁻¹ K⁻¹. `value_ref` carries the units
of the quantity being adjusted. All scientific parameters are supplied explicitly.
"""
Base.@kwdef struct Arrhenius{FT}
    value_ref::FT
    activation_energy::FT
    temperature_ref::FT
    gas_constant::FT
end

"""
    PeakedArrhenius(; arrhenius, deactivation_energy, entropy)

An Arrhenius response with high-temperature inhibition, normalized at the
same reference temperature. Deactivation energy: J mol⁻¹; entropy: J mol⁻¹ K⁻¹.
The three arguments must share their scalar type.
"""
struct PeakedArrhenius{FT}
    arrhenius::Arrhenius{FT}
    deactivation_energy::FT
    entropy::FT
end
PeakedArrhenius(; arrhenius::Arrhenius{FT}, deactivation_energy::FT,
    entropy::FT) where {FT} = PeakedArrhenius(arrhenius, deactivation_energy, entropy)

"""`temperature_response(response, temperature)` evaluates a response at temperature [K]."""
@inline function temperature_response(p::Arrhenius, temperature)
    exponent = p.activation_energy / p.gas_constant *
               (inv(p.temperature_ref) - inv(temperature))
    return p.value_ref * exp(exponent)
end

# Stable log(1 + exp(x)); 0 and 1 here are mathematical identities.
@inline softplus(x) = max(x, zero(x)) + log1p(exp(-abs(x)))

@inline function temperature_response(p::PeakedArrhenius, temperature)
    a = p.arrhenius
    inhibition_ref = (p.entropy - p.deactivation_energy / a.temperature_ref) / a.gas_constant
    inhibition = (p.entropy - p.deactivation_energy / temperature) / a.gas_constant
    exponent = a.activation_energy / a.gas_constant *
               (inv(a.temperature_ref) - inv(temperature))
    return a.value_ref * exp(exponent + softplus(inhibition_ref) - softplus(inhibition))
end

"""Six independent responses for C3 biochemical rates; see [`C3Rates`](@ref) for units."""
Base.@kwdef struct C3TemperatureResponse{V,J,R,C,O,G}
    vcmax::V
    jmax::J
    respiration::R
    kc::C
    ko::O
    gamma_star::G
end

"""
    C3Rates{FT}(; vcmax, jmax, respiration, kc, ko, gamma_star)

Parameters at the actual leaf temperature. Rates `vcmax`, `respiration`:
μmol CO₂ m⁻² leaf s⁻¹; `jmax`: μmol electrons m⁻² leaf s⁻¹.
`kc`, `ko`, `gamma_star`: Pa. Require nonnegative rates and gamma_star,
positive kc/ko, and positive assimilation denominators.
"""
Base.@kwdef struct C3Rates{FT}
    vcmax::FT
    jmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
end

"""Evaluate the six independent temperature responses, without mutating model parameters."""
@inline biochemical_rates(p::C3TemperatureResponse, temperature) = C3Rates(
    temperature_response(p.vcmax, temperature),
    temperature_response(p.jmax, temperature),
    temperature_response(p.respiration, temperature),
    temperature_response(p.kc, temperature),
    temperature_response(p.ko, temperature),
    temperature_response(p.gamma_star, temperature),
)
