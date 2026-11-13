"""
    AbstractOxygenExchange

Supertype for leaf oxygen-exchange components; positive net flux denotes O₂ release.
"""
abstract type AbstractOxygenExchange end

"""
    NoOxygenExchange()

Disable oxygen-budget calculations explicitly; [`oxygen_exchange`](@ref) returns `nothing`.
"""
struct NoOxygenExchange <: AbstractOxygenExchange end
@inline oxygen_exchange(::NoOxygenExchange, args...) = nothing

"""
    CoreOxygenExchange{FT}(; electrons_per_oxygen, oxygen_per_oxygenation,
        co2_per_oxygenation, respiratory_quotient)
    CoreOxygenExchange(electrons_per_oxygen, oxygen_per_oxygenation,
        co2_per_oxygenation, respiratory_quotient)

Specify dimensionless stoichiometric coefficients for PSII water-splitting O₂
evolution, Rubisco photorespiratory O₂ uptake, and mitochondrial respiratory O₂
uptake (`respiratory_quotient` is CO₂ released per O₂ consumed).

# Examples
```jldoctest
julia> ox = CoreOxygenExchange(4.0, 1.5, 0.5, 1.0);

julia> res = oxygen_exchange(ox, OxygenFluxInputs(80.0, 4.0, 1.0));

julia> res.net_release
13.0
```
"""
Base.@kwdef struct CoreOxygenExchange{FT<:Real} <: AbstractOxygenExchange
    electrons_per_oxygen::FT
    oxygen_per_oxygenation::FT
    co2_per_oxygenation::FT
    respiratory_quotient::FT
    function CoreOxygenExchange{FT}(
        electrons_per_oxygen,
        oxygen_per_oxygenation,
        co2_per_oxygenation,
        respiratory_quotient,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CoreOxygenExchange{FT},
            (
                electrons_per_oxygen,
                oxygen_per_oxygenation,
                co2_per_oxygenation,
                respiratory_quotient,
            ),
        )
        return new{FT}(
            electrons_per_oxygen,
            oxygen_per_oxygenation,
            co2_per_oxygenation,
            respiratory_quotient,
        )
    end
    function CoreOxygenExchange(
        electrons_per_oxygen::FT,
        oxygen_per_oxygenation::FT,
        co2_per_oxygenation::FT,
        respiratory_quotient::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CoreOxygenExchange,
            (
                electrons_per_oxygen,
                oxygen_per_oxygenation,
                co2_per_oxygenation,
                respiratory_quotient,
            ),
        )
        return new{FT}(
            electrons_per_oxygen,
            oxygen_per_oxygenation,
            co2_per_oxygenation,
            respiratory_quotient,
        )
    end
end
CoreOxygenExchange(epo::Real, opo::Real, cpo::Real, rq::Real) =
    CoreOxygenExchange(promote_compatible(CoreOxygenExchange, epo, opo, cpo, rq)...)

"""
    OxygenFluxInputs{FT}(; linear_electron_flow, rubisco_oxygenation, respiration)
    OxygenFluxInputs(linear_electron_flow, rubisco_oxygenation, respiration)

Bundle realized `linear_electron_flow` [μmol electrons m⁻² leaf s⁻¹],
`rubisco_oxygenation` [μmol reactions m⁻² leaf s⁻¹], and mitochondrial CO₂
`respiration` [μmol CO₂ m⁻² leaf s⁻¹] for [`oxygen_exchange`](@ref).
"""
Base.@kwdef struct OxygenFluxInputs{FT<:Real}
    linear_electron_flow::FT
    rubisco_oxygenation::FT
    respiration::FT
    function OxygenFluxInputs{FT}(
        linear_electron_flow,
        rubisco_oxygenation,
        respiration,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            OxygenFluxInputs{FT},
            (linear_electron_flow, rubisco_oxygenation, respiration),
        )
        return new{FT}(linear_electron_flow, rubisco_oxygenation, respiration)
    end
    function OxygenFluxInputs(
        linear_electron_flow::FT,
        rubisco_oxygenation::FT,
        respiration::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            OxygenFluxInputs,
            (linear_electron_flow, rubisco_oxygenation, respiration),
        )
        return new{FT}(linear_electron_flow, rubisco_oxygenation, respiration)
    end
end
OxygenFluxInputs(lef::Real, vo::Real, rd::Real) =
    OxygenFluxInputs(promote_compatible(OxygenFluxInputs, lef, vo, rd)...)

"""
    oxygen_exchange(oxygen_model, inputs::OxygenFluxInputs)
    oxygen_exchange(oxygen_model, photosynthesis, environment[, excitation])
    oxygen_exchange(leaf, environment[, excitation])

Compute gross O₂ evolution, photorespiratory O₂ uptake, mitochondrial
respiratory O₂ uptake, and net O₂ release [μmol O₂ m⁻² leaf s⁻¹] (`nothing` for
[`NoOxygenExchange`](@ref)).

When coupled to a signed [`C3Photosynthesis`](@ref) or
[`CytochromePhotosynthesis`](@ref) model, realized linear electron flow and
photorespiratory turnover are evaluated prior to net carbon cancellation so the
compensation point causes no `0/0` singularity.
"""
@inline function oxygen_exchange(p::CoreOxygenExchange, x::OxygenFluxInputs)
    gross = x.linear_electron_flow / p.electrons_per_oxygen
    photorespiration = p.oxygen_per_oxygenation * x.rubisco_oxygenation
    respiration = x.respiration / p.respiratory_quotient
    net = gross - photorespiration - respiration
    S = typeof(gross + photorespiration + respiration)
    return (;
        gross_production = S(gross),
        photorespiratory_uptake = S(photorespiration),
        respiratory_uptake = S(respiration),
        net_release = S(net),
    )
end

@inline function oxygen_exchange(
    p::CoreOxygenExchange,
    model::C3Photosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L}
    turnover = biochemical_turnover(model, env, excitation)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    inputs = OxygenFluxInputs(turnover.linear_electron_flow, oxygenation, turnover.respiration)
    ox = oxygen_exchange(p, inputs)
    OutFT = typeof(ox.net_release + turnover.carboxylation)
    return (;
        gross_production = OutFT(ox.gross_production),
        photorespiratory_uptake = OutFT(ox.photorespiratory_uptake),
        respiratory_uptake = OutFT(ox.respiratory_uptake),
        net_release = OutFT(ox.net_release),
        carboxylation = OutFT(turnover.carboxylation),
        rubisco_oxygenation = OutFT(oxygenation),
        linear_electron_flow = OutFT(turnover.linear_electron_flow),
    )
end
