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

@inline function invalid_carbon(::Type{FT}) where {FT}
    n = FT(NaN)
    return (; rubisco_limited = n, electron_limited = n, assimilation = n,
        respiration = n, net_assimilation = n)
end
@inline co2_excess(co2, gamma_star, ::SignedAssimilation) = co2 - gamma_star
@inline co2_excess(co2, gamma_star, ::NonnegativeAssimilation) =
    max(co2 - gamma_star, zero(co2 - gamma_star))

# Co-limit positive biochemical capacities before multiplying by CO₂ excess.
# This preserves turnover at the compensation point, including photorespiration.
@inline function assimilation_capacity(r, co2, oxygen, electrons, s, limitation)
    rubisco = r.vcmax / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons / (s.electrons_per_carboxylation * (co2 + s.oxygenation_cost * r.gamma_star))
    return colimit(rubisco, light, limitation)
end

"""
    carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, compensation)

C3 assimilation at prescribed biochemical-site CO₂ and O₂ partial pressures
[Pa]. `electrons` is μmol electrons m⁻² leaf s⁻¹. Returns Rubisco-limited,
electron-limited, combined pre-respiration, respiration, and net rates in
μmol CO₂ m⁻² leaf s⁻¹. Positive net means uptake. Pre-respiration rates already
include photorespiration; they are not raw Rubisco carboxylation fluxes.

Limitation acts on nonnegative biochemical capacities before multiplying by CO₂
excess. For signed sub-compensation fluxes this is an explicit extension of the
usual positive-uptake minimum: it limits turnover, not the most negative rate.

No stomatal or mesophyll diffusion is solved here. Equating this site's CO₂
with intercellular CO₂ assumes negligible mesophyll resistance.
"""
@inline function carbon_assimilation(r::C3Rates, co2, oxygen, electrons,
    s::C3Stoichiometry, limitation, compensation)
    excess = co2_excess(co2, r.gamma_star, compensation)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons * excess /
            (s.electrons_per_carboxylation * (co2 + s.oxygenation_cost * r.gamma_star))
    assimilation = assimilation_capacity(r, co2, oxygen, electrons, s, limitation) * excess
    return (; rubisco_limited = rubisco, electron_limited = light,
        assimilation, respiration = r.respiration, net_assimilation = assimilation - r.respiration)
end

"""
    biochemical_turnover(rates, co2, oxygen, electrons, stoichiometry, limitation, SignedAssimilation())

Realized C3 turnover before cancellation of carboxylation and photorespiratory
CO₂ release. Carbon fluxes and respiration are μmol CO₂ m⁻² leaf s⁻¹;
`linear_electron_flow` is μmol electrons m⁻² leaf s⁻¹. The supplied `electrons`
is available capacity, not necessarily realized flow.

Assumes linear electron flow serves carboxylation/photorespiration only. No
alternative sinks, cyclic flow, oxygen limitation of respiration or fluorescence
yield is inferred. Gamma-star must match the supplied oxygen. Require finite,
nonnegative inputs and positive biochemical denominators, as in C3 assimilation.
Clipped sub-compensation carbon models do not define this turnover contract.

Computing the positive capacity before multiplying by `co2 - gamma_star` keeps
turnover finite at the photorespiratory compensation point; no carbon/efficiency
ratio or small denominator floor is used.
"""
@inline function biochemical_turnover(r::C3Rates, co2, oxygen, electrons,
    s::C3Stoichiometry, limitation, ::SignedAssimilation)
    capacity = assimilation_capacity(r, co2, oxygen, electrons, s, limitation)
    return (; carboxylation = capacity * co2,
        photorespiratory_co2 = capacity * r.gamma_star,
        respiration = r.respiration,
        linear_electron_flow = s.electrons_per_carboxylation * capacity *
            (co2 + s.oxygenation_cost * r.gamma_star))
end
