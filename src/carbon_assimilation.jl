"""
    C3Stoichiometry{FT}(; electrons_per_carboxylation, oxygenation_cost)

Explicit stoichiometric parameters for the electron-limited denominator
`electrons_per_carboxylation * (co2 + oxygenation_cost * gamma_star)`.
Both are dimensionless and positive. Numerical values live in parameter files.
"""
Base.@kwdef struct C3Stoichiometry{FT}
    electrons_per_carboxylation::FT
    oxygenation_cost::FT
end

"""Retain the signed CO₂ response below the photorespiratory compensation point."""
struct SignedAssimilation end
"""Clip the CO₂ response at zero, as in Bonan supplemental program 11.1."""
struct NonnegativeAssimilation end
@inline co2_excess(co2, gamma_star, ::SignedAssimilation) = co2 - gamma_star
@inline co2_excess(co2, gamma_star, ::NonnegativeAssimilation) =
    max(co2 - gamma_star, zero(co2 - gamma_star))

"""
    carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)

C3 assimilation at prescribed biochemical-site CO₂ and O₂ partial pressures
[Pa]. `electrons` is μmol electrons m⁻² leaf s⁻¹. Returns Rubisco-limited,
electron-limited, combined pre-respiration, respiration, and net rates in
μmol CO₂ m⁻² leaf s⁻¹. Positive net means uptake. Pre-respiration rates already
include photorespiration; they are not raw Rubisco carboxylation fluxes.

No stomatal or mesophyll diffusion is solved here. Equating this site's CO₂
with intercellular CO₂ assumes negligible mesophyll resistance.
"""
@inline function carbon_assimilation(r::C3Rates, co2, oxygen, electrons,
    s::C3Stoichiometry, limitation, compensation)
    excess = co2_excess(co2, r.gamma_star, compensation)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons * excess /
            (s.electrons_per_carboxylation * (co2 + s.oxygenation_cost * r.gamma_star))
    assimilation = colimit(rubisco, light, limitation)
    return (; rubisco_limited = rubisco, electron_limited = light,
        assimilation, respiration = r.respiration, net_assimilation = assimilation - r.respiration)
end
