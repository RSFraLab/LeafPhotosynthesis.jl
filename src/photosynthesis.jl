"""Interface for a leaf photosynthesis component; concrete types implement `photosynthesize`."""
abstract type AbstractPhotosynthesis end

"""
    C3Photosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation)

Composable prescribed-CO₂ C3 biochemistry. Each physical/model choice is a
concrete field type. Contains parameters only, with no leaf arrays, scratch
cache, mutable iteration state, or hidden scientific defaults.
"""
Base.@kwdef struct C3Photosynthesis{T,E,S,L,C} <: AbstractPhotosynthesis
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    compensation::C
end

"""
    LeafEnvironment{FT}(; temperature, absorbed_par, co2, oxygen)

Prescribed leaf temperature [K], absorbed PAR [μmol photons m⁻² leaf s⁻¹],
and biochemical-site CO₂ and O₂ [Pa]. This is an input value, not prognostic
canopy air. Supply finite inputs, temperature > 0, other quantities ≥ 0.
Use one scalar type consistently throughout the model and environment.
"""
Base.@kwdef struct LeafEnvironment{FT}
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
end

"""
Evaluate one leaf environment. The `electron_transport` diagnostic is available
electron supply; `biochemical_turnover` separately reports realized `linear_electron_flow`.
"""
@inline function photosynthesize(p::C3Photosynthesis, env::LeafEnvironment)
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, rates.jmax)
    carbon = carbon_assimilation(rates, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation)
    return (; carbon..., electron_transport = electrons)
end

"""
    biochemical_turnover(photosynthesis, environment)

Signed-C3 carboxylation, photorespiratory CO₂ release, respiration and realized
linear electron flow, with `electron_supply` reported separately. Uses the
prescribed biochemical-site environment; does not solve stomatal diffusion.
Independent of whether optional oxygen exchange is enabled.
"""
@inline function biochemical_turnover(p::C3Photosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment) where {T,E,S,L}
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, r.jmax)
    turnover = biochemical_turnover(r, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation)
    return (; turnover..., electron_supply = electrons)
end
