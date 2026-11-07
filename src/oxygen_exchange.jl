"""Interface for leaf oxygen budgets; positive net flux means O₂ release."""
abstract type AbstractOxygenExchange end

"""Disable oxygen calculations explicitly; `oxygen_exchange` returns `nothing`."""
struct NoOxygenExchange <: AbstractOxygenExchange end
@inline oxygen_exchange(::NoOxygenExchange, args...) = nothing

"""
    CoreOxygenExchange{FT}(; electrons_per_oxygen, oxygen_per_oxygenation,
        co2_per_oxygenation, respiratory_quotient)

Core PSII/photorespiration/mitochondrial respiration budget. Coefficients are
explicit and dimensionless. RQ is CO₂ released / O₂ consumed. All must be positive.
Alternative electron sinks, cyclic flow, and O₂ limitation of respiration are not
represented. This budget does not itself solve internal O₂ diffusion/storage.
"""
Base.@kwdef struct CoreOxygenExchange{FT} <: AbstractOxygenExchange
    electrons_per_oxygen::FT
    oxygen_per_oxygenation::FT
    co2_per_oxygenation::FT
    respiratory_quotient::FT
end

"""
    OxygenFluxInputs{FT}(; linear_electron_flow, rubisco_oxygenation, respiration)

Actual linear electron flow [μmol electrons m⁻² leaf s⁻¹], Rubisco oxygenation
[μmol reactions m⁻² leaf s⁻¹], and mitochondrial CO₂ respiration
[μmol CO₂ m⁻² leaf s⁻¹]. All are nonnegative. Electron flow is realized flux,
not a potential transport capacity, and excludes cyclic flow.
"""
Base.@kwdef struct OxygenFluxInputs{FT}
    linear_electron_flow::FT
    rubisco_oxygenation::FT
    respiration::FT
end

"""Return gross production, separate consumption terms, and net O₂ release [μmol m⁻² leaf s⁻¹]."""
@inline function oxygen_exchange(p::CoreOxygenExchange, x::OxygenFluxInputs)
    gross = x.linear_electron_flow / p.electrons_per_oxygen
    photorespiration = p.oxygen_per_oxygenation * x.rubisco_oxygenation
    respiration = x.respiration / p.respiratory_quotient
    return (; gross_production = gross, photorespiratory_uptake = photorespiration,
        respiratory_uptake = respiration, net_release = gross - photorespiration - respiration)
end

"""
    oxygen_exchange(oxygen_model, photosynthesis, environment)

Close the core O₂ budget with a signed C3 model, assuming linear electron flow
serves carboxylation/photorespiration only. Calculate realized biochemical turnover
before net cancellation, so the compensation point does not cause a 0/0 singularity.
Clipped sub-compensation models are intentionally unsupported for this budget.
Gamma-star must be appropriate for the supplied oxygen level; its O₂ response is
not inferred from a temperature law.
"""
@inline function oxygen_exchange(p::CoreOxygenExchange,
    model::C3Photosynthesis{T,E,S,L,SignedAssimilation}, env::LeafEnvironment) where {T,E,S,L}
    turnover = biochemical_turnover(model, env)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    inputs = OxygenFluxInputs(turnover.linear_electron_flow, oxygenation, turnover.respiration)
    return (; oxygen_exchange(p, inputs)...,
        carboxylation = turnover.carboxylation, rubisco_oxygenation = oxygenation,
        linear_electron_flow = turnover.linear_electron_flow)
end
