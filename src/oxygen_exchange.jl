"""
    AbstractOxygenExchange

Supertype for leaf oxygen-exchange components; positive net flux denotes O₂ release.
"""
abstract type AbstractOxygenExchange <: AbstractLeafScalar end

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
julia> using LeafPhotosynthesis: OxygenFluxInputs

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
    CoreOxygenExchange{FT}(fields::Vararg{Any,4}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, CoreOxygenExchange{FT}, fields)...)
end
CoreOxygenExchange(fields::Vararg{Real,4}) =
    CoreOxygenExchange{scalar_type(CoreOxygenExchange, fields...)}(fields...)

"""
    OxygenFluxInputs{FT}(; linear_electron_flow, rubisco_oxygenation, respiration)
    OxygenFluxInputs(linear_electron_flow, rubisco_oxygenation, respiration)

Bundle realized `linear_electron_flow` [μmol electrons m⁻² leaf s⁻¹],
`rubisco_oxygenation` [μmol reactions m⁻² leaf s⁻¹], and mitochondrial CO₂
`respiration` [μmol CO₂ m⁻² leaf s⁻¹] for [`oxygen_exchange`](@ref).
"""
Base.@kwdef struct OxygenFluxInputs{FT<:Real} <: AbstractLeafScalar
    linear_electron_flow::FT
    rubisco_oxygenation::FT
    respiration::FT
    OxygenFluxInputs{FT}(fields::Vararg{Any,3}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, OxygenFluxInputs{FT}, fields)...)
end
OxygenFluxInputs(fields::Vararg{Real,3}) =
    OxygenFluxInputs{scalar_type(OxygenFluxInputs, fields...)}(fields...)

"""
    oxygen_exchange(oxygen_model, inputs::OxygenFluxInputs)
    oxygen_exchange(oxygen_model, photosynthesis, environment,
        excitation=BroadbandExcitation())
    oxygen_exchange(leaf, environment, excitation=BroadbandExcitation())

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
    return typed_result(typeof(gross + photorespiration + respiration), (;
        gross_production = gross,
        photorespiratory_uptake = photorespiration,
        respiratory_uptake = respiration,
        net_release = net,
    ))
end

@inline function oxygen_exchange(
    p::CoreOxygenExchange,
    model::AbstractC3Photosynthesis{L,SignedAssimilation},
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {L}
    # Gross O₂ evolution comes from PSII water oxidation (realized linear_electron_flow).
    turnover = biochemical_turnover(model, env, excitation)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    inputs = OxygenFluxInputs(turnover.linear_electron_flow, oxygenation, turnover.respiration)
    ox = oxygen_exchange(p, inputs)
    OutFT = typeof(ox.net_release + turnover.carboxylation)
    return typed_result(OutFT, (;
        ox...,
        turnover.carboxylation,
        rubisco_oxygenation = oxygenation,
        turnover.linear_electron_flow,
        turnover_electron_flows(model, turnover)...,
    ))
end
