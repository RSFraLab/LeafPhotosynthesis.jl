"""
    StateTransitionRule{FT}(; minimum_psii_fraction = zero(FT), maximum_psii_fraction = one(FT))

Closed-form state-transition antenna partitioning rule (`qT`; Johnson & Berry 2021,
Eqs. 25a, 28, 30a, Cases 10–12).

Under Cytochrome b₆f limitation (`J_P700,j ≤ J_P700,c`), solves analytically at each
point for the PSII absorbed-photon fraction `β = Q₂ / Q_abs` that balances linear
and PSI cyclic electron transport with zero dynamic feedback quenching (`K_N = 0`,
`Φ_P2 = Φ_P2^(K_N=0)(q₂)`). At and above the light-saturation crossover
(`J_P700,j ≥ J_P700,c`, where `q₂,sat = max(Vqmax - J_P700,c, 0) / Vqmax`), freezes
`β` at its crossover value `β_sat` so dynamic regulated heat dissipation (`K_N ≥ 0`)
absorbs excess excitation. That crossover freeze is exact for `MinimumLimitation()`
and remains `K_N ≥ 0` compatible for any concave `Ac`/`Aj` colimitation
(`QuadraticLimitation`, `SerialLimitation`).

Optional physiological LHCII mobility bounds clamp `β` to
`[minimum_psii_fraction, maximum_psii_fraction]` with `0 ≤ min < max ≤ 1`. Keep
`maximum_psii_fraction ≥ dark_balanced_psii_partition(...)` if `K_N ≥ 0` is required
down to `Q_abs → 0`; clamping `β` below the light-limited balance starves PSII and
causes `CytochromeFluorescence` to report `FluorescenceIncompatibleRates`.
"""
struct StateTransitionRule{FT<:AbstractFloat}
    minimum_psii_fraction::FT
    maximum_psii_fraction::FT
    function StateTransitionRule(
        minimum_psii_fraction::FT,
        maximum_psii_fraction::FT,
    ) where {FT<:AbstractFloat}
        isfinite(minimum_psii_fraction) && isfinite(maximum_psii_fraction) &&
            zero(FT) <= minimum_psii_fraction < maximum_psii_fraction <= one(FT) ||
            throw(ArgumentError("StateTransitionRule requires finite bounds 0 ≤ minimum_psii_fraction < maximum_psii_fraction ≤ 1"))
        new{FT}(minimum_psii_fraction, maximum_psii_fraction)
    end
end

StateTransitionRule{FT}(;
    minimum_psii_fraction::FT = zero(FT),
    maximum_psii_fraction::FT = one(FT),
) where {FT<:AbstractFloat} = StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)

StateTransitionRule(;
    minimum_psii_fraction::FT,
    maximum_psii_fraction::FT,
) where {FT<:AbstractFloat} = StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)

@inline compatible_state_scalar(::Type{FTR}, ::Type{FT}) where {FTR<:AbstractFloat,FT<:AbstractFloat} = FTR === FT
@inline compatible_state_scalar(::Type{FTR}, ::Type{FT}) where {FTR<:AbstractFloat,FT<:Real} = true

"""
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, temperature, co2, oxygen)
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, site::LeafEnvironment)

Low-light (`Q_abs → 0`, `q₂ → 1`, `K_N = 0`) PSII absorbed-photon fraction
`β₀ = Φ_P1,max / (Φ_P1,max + η(T, Cc, O₂) * Φ_P2,max(1))` that balances open-center
PSII and PSI photochemistry at `(temperature, co2, oxygen)` (Johnson & Berry 2021,
Eqs. 28, 30, Table 2). Nonthrowing on GPU; returns `NaN` for invalid inputs.
"""
@inline function dark_balanced_psii_partition(
    fl::CytochromeFluorescence{FTR},
    p::CytochromePhotosynthesis,
    temperature::FT,
    co2::FT,
    oxygen::FT,
) where {FTR<:AbstractFloat,FT<:Real}
    rates = biochemical_rates(p.temperature, temperature, oxygen)
    T = typeof(zero(FTR) + rates.vcmax + rates.vqmax + rates.gamma_star + temperature + co2 + oxygen)
    nan = T(NaN)
    compatible_state_scalar(FTR, FT) &&
        compatible_state_scalar(FTR, typeof(p.electron_transport.psi_maximum_yield)) &&
        compatible_state_scalar(FTR, typeof(p.stoichiometry.electrons_per_carboxylation)) &&
        compatible_state_scalar(FTR, typeof(rates.vcmax)) || return nan
    p.electron_transport.limitation isa SerialLimitation || return nan
    phi1_max = psi_maximum_photochemical_yield(fl)
    abs(p.electron_transport.psi_maximum_yield - phi1_max) <= fl.yield_tolerance || return nan
    isfinite(temperature) && temperature > zero(FT) &&
        isfinite(co2) && co2 >= zero(FT) &&
        isfinite(oxygen) && oxygen >= zero(FT) || return nan
    valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) &&
        !(iszero(co2) && iszero(rates.gamma_star)) || return nan
    eta = psi_to_psii_electron_ratio(rates, co2, p.stoichiometry)
    phi2_max = psii_maximum_photochemical_yield(fl, one(FTR))
    return T(phi1_max / (phi1_max + eta * phi2_max))
end

@inline function dark_balanced_psii_partition(
    fl::CytochromeFluorescence,
    p::CytochromePhotosynthesis,
    site::LeafEnvironment,
)
    val = dark_balanced_psii_partition(fl, p, site.temperature, site.co2, site.oxygen)
    return isfinite(site.absorbed_par) && site.absorbed_par >= zero(site.absorbed_par) ?
        val : oftype(val, NaN)
end

"""
    state_transition_partition(rule::StateTransitionRule,
        fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis,
        site::LeafEnvironment)

Pointwise, closed-form PSII absorbed-photon fraction `β ∈ (0, 1)` under
`StateTransitionRule` (Johnson & Berry 2021, Cases 10–12). Uses scale-normalized
dimensionless rate and irradiance ratios so `Float32`, `Float64`, GPU kernels, and
automatic-differentiation dual inputs evaluate without iteration, array scans,
overflow, or low-light subtractive cancellation. Returns `NaN` for invalid inputs
or mismatched `psi_maximum_yield`.
"""
@inline function state_transition_partition(
    rule::StateTransitionRule{FTR},
    fl::CytochromeFluorescence{FTR},
    p::CytochromePhotosynthesis,
    site::LeafEnvironment{FT},
) where {FTR<:AbstractFloat,FT<:Real}
    rates = biochemical_rates(p.temperature, site.temperature, site.oxygen)
    T = typeof(zero(FTR) + rates.vcmax + rates.vqmax + rates.gamma_star +
               site.temperature + site.absorbed_par + site.co2 + site.oxygen)
    nan = T(NaN)
    compatible_state_scalar(FTR, FT) &&
        compatible_state_scalar(FTR, typeof(p.electron_transport.psi_maximum_yield)) &&
        compatible_state_scalar(FTR, typeof(p.stoichiometry.electrons_per_carboxylation)) &&
        compatible_state_scalar(FTR, typeof(rates.vcmax)) || return nan
    p.electron_transport.limitation isa SerialLimitation || return nan
    phi1_max = psi_maximum_photochemical_yield(fl)
    abs(p.electron_transport.psi_maximum_yield - phi1_max) <= fl.yield_tolerance || return nan
    isfinite(site.temperature) && site.temperature > zero(FT) &&
        isfinite(site.absorbed_par) && site.absorbed_par >= zero(FT) &&
        isfinite(site.co2) && site.co2 >= zero(FT) &&
        isfinite(site.oxygen) && site.oxygen >= zero(FT) || return nan

    valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) &&
        !(iszero(site.co2) && iszero(rates.gamma_star)) || return nan

    # Dimensionless PSII rate ratios bounded in [0, 1] (zero rate^2 overflow in Float32).
    psii = fl.psii
    d0 = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
    s_base = d0 + psii.exciton_transfer_rate
    s_tot = s_base + psii.maximum_photochemical_rate
    d0_hat = d0 / s_tot
    ku_hat = psii.exciton_transfer_rate / s_tot
    kp_hat = psii.maximum_photochemical_rate / s_tot
    sb_hat = s_base / s_tot

    eta = psi_to_psii_electron_ratio(rates, site.co2, p.stoichiometry)
    m1 = phi1_max * d0_hat
    mu = phi1_max * kp_hat * ku_hat
    m2 = eta * kp_hat * sb_hat

    # Scale-normalized quadratic for light-limited α₁ = 1 - β (K_N = 0):
    # w_q (m₁ + m₂) α₁² + [w_v (m₁ + m_u + m₂) - w_q m₂] α₁ - w_v m₂ = 0.
    q_light = phi1_max * site.absorbed_par
    s_vq = max(rates.vqmax, q_light)
    alpha1_opt = if iszero(s_vq)
        m2 / (m1 + mu + m2)
    else
        w_v = rates.vqmax / s_vq
        w_q = q_light / s_vq
        a_coef = w_q * (m1 + m2)
        b_coef = w_v * (m1 + mu + m2) - w_q * m2
        c_pos = w_v * m2
        disc = sqrt(b_coef * b_coef + 4 * a_coef * c_pos)
        b_coef >= zero(b_coef) ? (2 * c_pos) / (b_coef + disc) : (-b_coef + disc) / (2 * a_coef)
    end
    beta_opt = one(T) - alpha1_opt

    # Closed-form light-saturation floor β_sat at J_P700,j = J_P700,c (q₂,sat = max(Vqmax - J_P700,c, 0) / Vqmax).
    c1, c2 = cytochrome_carbon_coefficients(rates, p.stoichiometry)
    km = site.co2 + rates.kc * (one(T) + site.oxygen / rates.ko)
    j_psi_c = (rates.vcmax / km) * (c1 * site.co2 + c2 * rates.gamma_star)
    q2_sat = rates.vqmax > zero(rates.vqmax) ?
        max(rates.vqmax - j_psi_c, zero(T)) / rates.vqmax : zero(T)
    num_sat = m1 + q2_sat * mu
    beta_sat = num_sat / (num_sat + m2)

    beta_raw = max(beta_opt, beta_sat)
    return T(clamp(beta_raw, T(rule.minimum_psii_fraction), T(rule.maximum_psii_fraction)))
end

@inline state_transition_partition(
    rule::StateTransitionRule,
    fl::CytochromeFluorescence,
    leaf::LeafModel,
    site::LeafEnvironment,
) = state_transition_partition(rule, fl, leaf.photosynthesis, site)

"""
    state_transition_excitation(rule::StateTransitionRule,
        fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis,
        site::LeafEnvironment)

Return `PrescribedDualExcitation(Q₂, Q₁)` implementing Johnson & Berry (2021)
dynamic state transitions (`qT`) at `site`. Uses a Sterbenz-exact split of
`site.absorbed_par` so `Q₂ + Q₁ == site.absorbed_par` holds identically across all
IEEE-754 binades without clipping or automatic-differentiation tangent distortion.
"""
@inline function state_transition_excitation(
    rule::StateTransitionRule{FTR},
    fl::CytochromeFluorescence{FTR},
    p::CytochromePhotosynthesis,
    site::LeafEnvironment{FT},
) where {FTR<:AbstractFloat,FT<:Real}
    beta = state_transition_partition(rule, fl, p, site)
    T_q = typeof(beta * site.absorbed_par)
    q_abs = T_q(site.absorbed_par)
    if !isfinite(beta)
        nan = T_q(NaN)
        return PrescribedDualExcitation(nan, nan)
    end
    half = inv(oftype(beta, 2))
    if beta >= half
        q2 = beta * q_abs
        q1 = q_abs - q2
    else
        q1 = (one(beta) - beta) * q_abs
        q2 = q_abs - q1
    end
    return PrescribedDualExcitation(q2, q1)
end

@inline state_transition_excitation(
    rule::StateTransitionRule,
    fl::CytochromeFluorescence,
    leaf::LeafModel,
    site::LeafEnvironment,
) = state_transition_excitation(rule, fl, leaf.photosynthesis, site)
