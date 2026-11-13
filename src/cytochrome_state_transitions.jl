"""
    StateTransitionRule{FT}(; minimum_psii_fraction=zero(FT), maximum_psii_fraction=one(FT))
    StateTransitionRule(; minimum_psii_fraction, maximum_psii_fraction)
    StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)

Configure closed-form state-transition antenna partitioning (`qT`; Johnson &
Berry 2021, Eqs. 25a, 28, 30a, Cases 10–12) with optional LHCII mobility bounds
`0 ≤ minimum_psii_fraction < maximum_psii_fraction ≤ 1`.

Under Cytochrome b₆f limitation (`J_P700,j ≤ J_P700,c`), solves analytically for
the PSII absorbed-photon fraction `β = Q₂ / Q_abs` that balances linear and
cyclic electron transport at `K_N = 0`, freezing `β` at its light-saturation
crossover value `β_sat` above `J_P700,j ≥ J_P700,c`.
"""
struct StateTransitionRule{FT<:Real}
    minimum_psii_fraction::FT
    maximum_psii_fraction::FT
    function StateTransitionRule(
        minimum_psii_fraction::FT,
        maximum_psii_fraction::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, StateTransitionRule, (minimum_psii_fraction, maximum_psii_fraction),
        )
        lo_p = primal_scalar(minimum_psii_fraction)
        hi_p = primal_scalar(maximum_psii_fraction)
        isfinite(lo_p) && isfinite(hi_p) &&
            zero(lo_p) <= lo_p < hi_p <= one(hi_p) ||
            throw(ArgumentError("StateTransitionRule requires finite bounds 0 ≤ minimum_psii_fraction < maximum_psii_fraction ≤ 1"))
        return new{FT}(minimum_psii_fraction, maximum_psii_fraction)
    end
end

StateTransitionRule(minimum_psii_fraction::Real, maximum_psii_fraction::Real) =
    StateTransitionRule(
        promote_compatible(StateTransitionRule, minimum_psii_fraction, maximum_psii_fraction)...,
    )

StateTransitionRule{FT}(;
    minimum_psii_fraction::Real = zero(FT),
    maximum_psii_fraction::Real = one(FT),
) where {FT<:Real} = StateTransitionRule(FT(minimum_psii_fraction), FT(maximum_psii_fraction))

StateTransitionRule(; minimum_psii_fraction, maximum_psii_fraction) =
    StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)

"""
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, temperature, co2, oxygen)
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, site::LeafEnvironment)

Compute the low-light (`Q_abs → 0`, `q₂ → 1`, `K_N = 0`) PSII absorbed-photon
fraction ``β_0 = Φ_{P1,\\max} / (Φ_{P1,\\max} + η(T, C_c, O_2) Φ_{P2,\\max}(1))``
that balances open-center PSII and PSI photochemistry (Johnson & Berry 2021,
Eqs. 28, 30, Table 2), returning `NaN` on invalid inputs without throwing.
"""
@inline function dark_balanced_psii_partition(
    fl::CytochromeFluorescence{FTR},
    p::CytochromePhotosynthesis,
    temperature::FTT,
    co2::FTC,
    oxygen::FTO,
) where {FTR<:Real,FTT<:Real,FTC<:Real,FTO<:Real}
    rates = biochemical_rates(p.temperature, temperature, oxygen)
    T = typeof(zero(FTR) + rates.vcmax + rates.vqmax + rates.gamma_star + temperature + co2 + oxygen)
    nan = T(NaN)
    compatible_scalar(FTR, FTT) &&
        compatible_scalar(FTR, FTC) &&
        compatible_scalar(FTR, FTO) &&
        compatible_scalar(FTR, typeof(p.electron_transport.psi_maximum_yield)) &&
        compatible_scalar(FTR, typeof(p.stoichiometry.electrons_per_carboxylation)) &&
        compatible_scalar(FTR, typeof(rates.vcmax)) || return nan
    p.electron_transport.limitation isa SerialLimitation || return nan
    phi1_max = psi_maximum_photochemical_yield(fl)
    primal_scalar(abs(p.electron_transport.psi_maximum_yield - phi1_max)) <= primal_scalar(fl.yield_tolerance) || return nan
    t_p, co2_p, o2_p = map(primal_scalar, (temperature, co2, oxygen))
    isfinite(t_p) && t_p > zero(t_p) &&
        isfinite(co2_p) && co2_p >= zero(co2_p) &&
        isfinite(o2_p) && o2_p >= zero(o2_p) || return nan
    valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) &&
        !(iszero(co2_p) && iszero(primal_scalar(rates.gamma_star))) || return nan
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
    par_p = primal_scalar(site.absorbed_par)
    return isfinite(par_p) && par_p >= zero(par_p) ?
        val : oftype(val, NaN)
end

"""
    state_transition_partition(rule::StateTransitionRule,
        fluorescence_model::CytochromeFluorescence,
        photosynthesis_or_leaf, site::LeafEnvironment)

Evaluate the pointwise, closed-form PSII absorbed-photon fraction `β ∈ (0, 1)`
under `rule` at `site` (Johnson & Berry 2021, Cases 10–12), returning `NaN` on
invalid or incompatible inputs without throwing.
"""
@inline function state_transition_partition(
    rule::StateTransitionRule{FTR},
    fl::CytochromeFluorescence{FTF},
    p::CytochromePhotosynthesis,
    site::LeafEnvironment{FT},
) where {FTR<:Real,FTF<:Real,FT<:Real}
    rates = biochemical_rates(p.temperature, site.temperature, site.oxygen)
    T = typeof(
        zero(FTR) + zero(FTF) + rates.vcmax + rates.vqmax + rates.gamma_star +
        site.temperature + site.absorbed_par + site.co2 + site.oxygen,
    )
    nan = T(NaN)
    compatible_scalar(FTR, FTF) &&
        compatible_scalar(FTR, FT) &&
        compatible_scalar(FTR, typeof(p.electron_transport.psi_maximum_yield)) &&
        compatible_scalar(FTR, typeof(p.stoichiometry.electrons_per_carboxylation)) &&
        compatible_scalar(FTR, typeof(rates.vcmax)) || return nan
    p.electron_transport.limitation isa SerialLimitation || return nan
    phi1_max = psi_maximum_photochemical_yield(fl)
    primal_scalar(abs(p.electron_transport.psi_maximum_yield - phi1_max)) <= primal_scalar(fl.yield_tolerance) || return nan
    t_p, par_p, co2_p, o2_p = map(
        primal_scalar,
        (site.temperature, site.absorbed_par, site.co2, site.oxygen),
    )
    isfinite(t_p) && t_p > zero(t_p) &&
        isfinite(par_p) && par_p >= zero(par_p) &&
        isfinite(co2_p) && co2_p >= zero(co2_p) &&
        isfinite(o2_p) && o2_p >= zero(o2_p) || return nan

    valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) &&
        !(iszero(co2_p) && iszero(primal_scalar(rates.gamma_star))) || return nan

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
    alpha1_opt = if iszero(primal_scalar(s_vq))
        m2 / (m1 + mu + m2)
    else
        w_v = rates.vqmax / s_vq
        w_q = q_light / s_vq
        a_coef = w_q * (m1 + m2)
        b_coef = w_v * (m1 + mu + m2) - w_q * m2
        c_pos = w_v * m2
        rad2 = b_coef * b_coef + 4 * a_coef * c_pos
        disc = rad2 > zero(rad2) ? sqrt(rad2) : zero(rad2)
        b_coef >= zero(b_coef) ? (2 * c_pos) / (b_coef + disc) : (-b_coef + disc) / (2 * a_coef)
    end
    beta_opt = one(T) - alpha1_opt

    # Closed-form light-saturation floor β_sat at J_P700,j = J_P700,c.
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
        photosynthesis_or_leaf, site::LeafEnvironment)

Return a [`PrescribedDualExcitation`](@ref) `(Q₂, Q₁)` implementing Johnson &
Berry (2021) dynamic state transitions (`qT`) at `site` via a Sterbenz-exact
split of `site.absorbed_par` so `Q₂ + Q₁ == site.absorbed_par` holds identically
across all IEEE-754 binades.
"""
@inline function state_transition_excitation(
    rule::StateTransitionRule{FTR},
    fl::CytochromeFluorescence{FTF},
    p::CytochromePhotosynthesis,
    site::LeafEnvironment{FT},
) where {FTR<:Real,FTF<:Real,FT<:Real}
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
