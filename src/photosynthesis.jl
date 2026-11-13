"""
    AbstractPhotosynthesis

Supertype for leaf photosynthetic biochemistry components implementing [`photosynthesize`](@ref).
"""
abstract type AbstractPhotosynthesis end

"""
    C3Photosynthesis(; temperature, electron_transport, stoichiometry, limitation, compensation)

Compose Farquhar–von Caemmerer–Berry C3 biochemistry from concrete `temperature`,
`electron_transport`, `stoichiometry`, `limitation`, and `compensation` components.

`C3Photosynthesis` holds immutable parameters only, with no per-leaf arrays,
scratch buffers, or hidden empirical defaults.
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
    LeafEnvironment(temperature, absorbed_par, co2, oxygen)

Represent prescribed leaf `temperature` [K], `absorbed_par` [μmol photons m⁻² leaf s⁻¹],
and biochemical-site `co2` and `oxygen` partial pressures [Pa].
"""
Base.@kwdef struct LeafEnvironment{FT<:Real}
    temperature::FT
    absorbed_par::FT
    co2::FT
    oxygen::FT
    function LeafEnvironment{FT}(temperature, absorbed_par, co2, oxygen) where {FT<:Real}
        require_continuous_scalar(
            FT, LeafEnvironment{FT}, (temperature, absorbed_par, co2, oxygen),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen)
    end
    function LeafEnvironment(
        temperature::FT,
        absorbed_par::FT,
        co2::FT,
        oxygen::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, LeafEnvironment, (temperature, absorbed_par, co2, oxygen),
        )
        return new{FT}(temperature, absorbed_par, co2, oxygen)
    end
end
LeafEnvironment(t::Real, par::Real, co2::Real, o2::Real) =
    LeafEnvironment(promote_compatible(LeafEnvironment, t, par, co2, o2)...)

"""
    photosynthesize(photosynthesis, environment[, excitation])
    photosynthesize(leaf, environment[, excitation])

Evaluate leaf photosynthetic carbon assimilation and available electron supply
for `environment` using `BroadbandExcitation()` by default or an explicit
[`AbstractExcitation`](@ref).

The returned `electron_transport` field reports available electron supply; use
[`biochemical_turnover`](@ref) for realized `linear_electron_flow`.
"""
@inline function photosynthesize(
    p::C3Photosynthesis,
    env::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {FT<:Real}
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, rates.jmax, excitation)
    S = typeof(rates.vcmax + p.stoichiometry.electrons_per_carboxylation + env.co2 + env.oxygen + electrons)
    valid_excitation(excitation, env.absorbed_par) ||
        return (; invalid_carbon(S)..., electron_transport = S(NaN))
    carbon = carbon_assimilation(
        rates, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation,
    )
    return (;
        rubisco_limited = S(carbon.rubisco_limited),
        electron_limited = S(carbon.electron_limited),
        assimilation = S(carbon.assimilation),
        respiration = S(carbon.respiration),
        net_assimilation = S(carbon.net_assimilation),
        electron_transport = S(electrons),
    )
end

@inline function biochemical_turnover(
    p::C3Photosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L,FT<:Real}
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    electrons = electron_transport(p.electron_transport, env.absorbed_par, r.jmax, excitation)
    OutFT = typeof(r.vcmax + p.stoichiometry.electrons_per_carboxylation + env.co2 + env.oxygen + electrons)
    if !valid_excitation(excitation, env.absorbed_par)
        n = OutFT(NaN)
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, electron_supply = n)
    end
    turnover = biochemical_turnover(
        r, env.co2, env.oxygen, electrons,
        p.stoichiometry, p.limitation, p.compensation,
    )
    return (;
        carboxylation = OutFT(turnover.carboxylation),
        photorespiratory_co2 = OutFT(turnover.photorespiratory_co2),
        respiration = OutFT(turnover.respiration),
        linear_electron_flow = OutFT(turnover.linear_electron_flow),
        electron_supply = OutFT(electrons),
    )
end
