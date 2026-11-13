"""
    C3Stoichiometry{FT}(; electrons_per_carboxylation, oxygenation_cost)
    C3Stoichiometry(electrons_per_carboxylation, oxygenation_cost)

Specify dimensionless C3 stoichiometric coefficients in the electron-limited
denominator ``e_c (C_c + c_o Γ^*)``.

Both parameters must be supplied explicitly from parameter files.
"""
Base.@kwdef struct C3Stoichiometry{FT<:Real}
    electrons_per_carboxylation::FT
    oxygenation_cost::FT
    function C3Stoichiometry{FT}(electrons_per_carboxylation, oxygenation_cost) where {FT<:Real}
        require_continuous_scalar(
            FT, C3Stoichiometry{FT}, (electrons_per_carboxylation, oxygenation_cost),
        )
        return new{FT}(electrons_per_carboxylation, oxygenation_cost)
    end
    function C3Stoichiometry(
        electrons_per_carboxylation::FT,
        oxygenation_cost::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, C3Stoichiometry, (electrons_per_carboxylation, oxygenation_cost),
        )
        return new{FT}(electrons_per_carboxylation, oxygenation_cost)
    end
end
C3Stoichiometry(ec::Real, oc::Real) =
    C3Stoichiometry(promote_compatible(C3Stoichiometry, ec, oc)...)

"""
    SignedAssimilation()

Retain signed CO₂ exchange below the photorespiratory compensation point `Γ*`.
"""
struct SignedAssimilation end

"""
    NonnegativeAssimilation()

Clip pre-respiration CO₂ assimilation at zero below `Γ*`, as in Bonan
supplemental program 11.1.
"""
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

Compute C3 carbon assimilation fluxes [μmol CO₂ m⁻² leaf s⁻¹] at prescribed
biochemical-site `co2` and `oxygen` partial pressures [Pa] and available
electron flux `electrons` [μmol electrons m⁻² leaf s⁻¹].

Return `(; rubisco_limited, electron_limited, assimilation, respiration, net_assimilation)`,
where positive `net_assimilation` denotes net CO₂ uptake and `assimilation`
already accounts for photorespiratory CO₂ release prior to subtracting
mitochondrial `respiration`. Co-limitation acts on nonnegative biochemical
capacities before multiplying by CO₂ excess so realized turnover remains
well-defined across `Γ*`.
"""
@inline function carbon_assimilation(
    r::C3Rates,
    co2,
    oxygen,
    electrons,
    s::C3Stoichiometry,
    limitation,
    compensation,
)
    excess = co2_excess(co2, r.gamma_star, compensation)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    light = electrons * excess /
            (s.electrons_per_carboxylation * (co2 + s.oxygenation_cost * r.gamma_star))
    assimilation = assimilation_capacity(r, co2, oxygen, electrons, s, limitation) * excess
    net = assimilation - r.respiration
    S = typeof(rubisco + light + net)
    return (;
        rubisco_limited = S(rubisco),
        electron_limited = S(light),
        assimilation = S(assimilation),
        respiration = S(r.respiration),
        net_assimilation = S(net),
    )
end

"""
    biochemical_turnover(rates, co2, oxygen, electrons, stoichiometry, limitation, SignedAssimilation())
    biochemical_turnover(photosynthesis, environment[, excitation])
    biochemical_turnover(leaf, environment[, excitation])

Compute realized C3 gross biochemical turnover prior to cancellation between
Rubisco carboxylation and photorespiratory CO₂ release.

Return `(; carboxylation, photorespiratory_co2, respiration, linear_electron_flow)`
(plus `electron_supply` on environment entry points, and `psi_electron_flow`,
`cyclic_electron_flow`, `psi_electron_supply`, `psii_open_fraction` on
[`CytochromePhotosynthesis`](@ref)), with carbon fluxes in μmol CO₂ m⁻² leaf s⁻¹
and electron fluxes in μmol electrons m⁻² leaf s⁻¹. Only [`SignedAssimilation`](@ref)
defines this turnover contract.
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
    S = typeof(carboxylation + photorespiratory_co2 + linear_flow + r.respiration)
    return (;
        carboxylation = S(carboxylation),
        photorespiratory_co2 = S(photorespiratory_co2),
        respiration = S(r.respiration),
        linear_electron_flow = S(linear_flow),
    )
end
