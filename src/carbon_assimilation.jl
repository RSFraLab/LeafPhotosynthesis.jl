"""
    C3Stoichiometry{FT}(; electrons_per_carboxylation, oxygenation_cost)
    C3Stoichiometry(electrons_per_carboxylation, oxygenation_cost)

Specify dimensionless C3 stoichiometric coefficients in the electron-limited
denominator ``e_c (C_c + c_o Γ^*)``.

Both parameters must be supplied explicitly from parameter files. Throws
`ArgumentError` unless `electrons_per_carboxylation > 0` and `oxygenation_cost ≥ 0` are
finite.
"""
Base.@kwdef struct C3Stoichiometry{FT<:Real} <: AbstractLeafScalar
    electrons_per_carboxylation::FT
    oxygenation_cost::FT
    function C3Stoichiometry{FT}(fields::Vararg{Any,2}) where {FT<:Real}
        s = new{FT}(continuous_fields(FT, C3Stoichiometry{FT}, fields)...)
        is_positive(s.electrons_per_carboxylation) && is_nonnegative(s.oxygenation_cost) ||
            throw(ArgumentError("C3Stoichiometry requires finite electrons_per_carboxylation > 0 and oxygenation_cost ≥ 0"))
        return s
    end
end
C3Stoichiometry(fields::Vararg{Real,2}) =
    C3Stoichiometry{scalar_type(C3Stoichiometry, fields...)}(fields...)

"""
    AbstractAssimilationSign

Supertype of the treatments of CO₂ assimilation below the photorespiratory compensation
point `Γ*`: [`SignedAssimilation`](@ref) keeps its sign, [`NonnegativeAssimilation`](@ref)
clips it at zero.
"""
abstract type AbstractAssimilationSign <: AbstractLeafScalar end

"""
    SignedAssimilation()

Retain signed CO₂ exchange below the photorespiratory compensation point `Γ*`.
"""
struct SignedAssimilation <: AbstractAssimilationSign end

"""
    NonnegativeAssimilation()

Clip pre-respiration CO₂ assimilation at zero below `Γ*`, as in Bonan
supplemental program 11.1.
"""
struct NonnegativeAssimilation <: AbstractAssimilationSign end

@inline function invalid_carbon(::Type{FT}) where {FT}
    n = FT(NaN)
    return (; rubisco_limited = n, electron_limited = n, assimilation = n,
        respiration = n, net_assimilation = n)
end
@inline co2_excess(co2, gamma_star, ::SignedAssimilation) = co2 - gamma_star
@inline co2_excess(co2, gamma_star, ::NonnegativeAssimilation) =
    primal_max(co2 - gamma_star, zero(co2 - gamma_star))

# Electrons consumed per unit of carboxylation capacity [e⁻ Pa]: the light-limited
# capacity is electrons / electron_demand, e (C + a Γ*) for C3 linear electron flow.
@inline electron_demand(r::C3Rates, co2, s::C3Stoichiometry) =
    s.electrons_per_carboxylation * (co2 + s.oxygenation_cost * r.gamma_star)

# Co-limit positive biochemical capacities before multiplying by CO₂ excess.
# This preserves turnover at the compensation point, including photorespiration.
@inline function assimilation_capacity(r, co2, oxygen, electrons, s, limitation)
    rubisco = r.vcmax / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons / electron_demand(r, co2, s)
    return colimit(rubisco, light, limitation)
end

"""
    carbon_assimilation(rates, co2, oxygen, electrons, stoichiometry, limitation, assimilation_sign)

Compute C3 carbon assimilation fluxes [μmol CO₂ m⁻² leaf s⁻¹] at prescribed
biochemical-site `co2` and `oxygen` partial pressures [Pa] and available
electron flux `electrons` [μmol electrons m⁻² leaf s⁻¹].

The electron-limited flux is ``J (C_c - Γ^*) / D``, with electron demand
``D = e_c (C_c + a Γ^*)`` for [`C3Rates`](@ref) and ``D = c_1 C_c + c_2 Γ^*`` (PSI
electrons, see [`cytochrome_carbon_coefficients`](@ref)) for
[`CytochromeRates`](@ref).

Return `(; rubisco_limited, electron_limited, assimilation, respiration, net_assimilation)`,
where positive `net_assimilation` denotes net CO₂ uptake and `assimilation`
already accounts for photorespiratory CO₂ release prior to subtracting
mitochondrial `respiration`. Co-limitation acts on nonnegative biochemical
capacities before multiplying by CO₂ excess so realized turnover remains
well-defined across `Γ*`.
"""
@inline function carbon_assimilation(r, co2, oxygen, electrons, s, limitation, assimilation_sign)
    require_compatible(carbon_assimilation, r, co2, oxygen, electrons, s, limitation)
    excess = co2_excess(co2, r.gamma_star, assimilation_sign)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons * excess / electron_demand(r, co2, s)
    assimilation = assimilation_capacity(r, co2, oxygen, electrons, s, limitation) * excess
    net = assimilation - r.respiration
    return typed_result(typeof(rubisco + light + net), (;
        rubisco_limited = rubisco,
        electron_limited = light,
        assimilation,
        r.respiration,
        net_assimilation = net,
    ))
end

"""
    biochemical_turnover(rates, co2, oxygen, electrons, stoichiometry, limitation, SignedAssimilation())
    biochemical_turnover(photosynthesis, environment, excitation=BroadbandExcitation())
    biochemical_turnover(leaf, environment, excitation=BroadbandExcitation())

Compute realized C3 gross biochemical turnover prior to cancellation between
Rubisco carboxylation and photorespiratory CO₂ release.

Return `(; carboxylation, photorespiratory_co2, respiration, linear_electron_flow)`,
with carbon fluxes in μmol CO₂ m⁻² leaf s⁻¹ and electron fluxes in μmol electrons
m⁻² leaf s⁻¹. [`CytochromePhotosynthesis`](@ref) adds the realized `psi_electron_flow`
and `cyclic_electron_flow`. Environment entry points append the available linear
`electron_transport` and, for [`CytochromePhotosynthesis`](@ref), `psi_electron_supply`,
`psii_open_fraction` and `psi_to_psii_ratio`, as [`photosynthesize`](@ref) does; input
outside the physical domain gives `NaN` in every field. Only
[`SignedAssimilation`](@ref) defines this turnover contract.
"""
@inline function biochemical_turnover(
    r::C3Rates,
    co2,
    oxygen,
    electrons,
    s::C3Stoichiometry,
    limitation,
    ::SignedAssimilation,
)
    capacity = assimilation_capacity(r, co2, oxygen, electrons, s, limitation)
    carboxylation = capacity * co2
    photorespiratory_co2 = capacity * r.gamma_star
    linear_flow = s.electrons_per_carboxylation * capacity *
                  (co2 + s.oxygenation_cost * r.gamma_star)
    return typed_result(typeof(carboxylation + photorespiratory_co2 + linear_flow + r.respiration), (;
        carboxylation,
        photorespiratory_co2,
        r.respiration,
        linear_electron_flow = linear_flow,
    ))
end
