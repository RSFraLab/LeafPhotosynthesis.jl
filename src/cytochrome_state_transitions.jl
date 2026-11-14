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
struct StateTransitionRule{FT<:Real} <: AbstractLeafScalar
    minimum_psii_fraction::FT
    maximum_psii_fraction::FT
    function StateTransitionRule{FT}(fields::Vararg{Any,2}) where {FT<:Real}
        r = new{FT}(continuous_fields(FT, StateTransitionRule{FT}, fields)...)
        lo, hi = r.minimum_psii_fraction, r.maximum_psii_fraction
        is_fraction(lo) && is_fraction(hi) && exceeds(hi, lo) ||
            throw(ArgumentError("StateTransitionRule requires finite bounds 0 ≤ minimum_psii_fraction < maximum_psii_fraction ≤ 1"))
        return r
    end
end
StateTransitionRule(fields::Vararg{Real,2}) =
    StateTransitionRule{scalar_type(StateTransitionRule, fields...)}(fields...)

StateTransitionRule{FT}(;
    minimum_psii_fraction::Real = zero(FT),
    maximum_psii_fraction::Real = one(FT),
) where {FT<:Real} = StateTransitionRule{FT}(minimum_psii_fraction, maximum_psii_fraction)

StateTransitionRule(; minimum_psii_fraction, maximum_psii_fraction) =
    StateTransitionRule(minimum_psii_fraction, maximum_psii_fraction)

# Domain shared by the partition closed forms: PSI maximum yields of the photosynthesis and
# fluorescence models agree within the tolerance, a valid site and rate domain, and a
# defined PSI/PSII electron ratio (not C = Γ* = 0).
@inline valid_partition_inputs(fl, p, rates, phi1_max, temperature, co2, oxygen) =
    is_at_most(abs(p.electron_transport.psi_maximum_yield - phi1_max), fl.yield_tolerance) &&
    valid_site(temperature, co2, oxygen) && valid_biochemistry(p, rates) &&
    !(is_zero(co2) && is_zero(rates.gamma_star))

"""
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, temperature, co2, oxygen)
    dark_balanced_psii_partition(fluorescence_model::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, site::LeafEnvironment)

Compute the low-light (`Q_abs → 0`, `q₂ → 1`, `K_N = 0`) PSII absorbed-photon
fraction ``β_0 = Φ_{P1,\\max} / (Φ_{P1,\\max} + η(T, C_c, O_2) Φ_{P2,\\max}(1))``
that balances open-center PSII and PSI photochemistry (Johnson & Berry 2021,
Eqs. 28, 30, Table 2), returning `NaN` on invalid inputs without throwing. The closed form
holds for [`SerialLimitation`](@ref) in `photosynthesis.electron_transport`; another
co-limitation is a `MethodError`.
"""
@inline function dark_balanced_psii_partition(
    fl::CytochromeFluorescence{FTR},
    p::SerialCytochromePhotosynthesis,
    temperature::FTT,
    co2::FTC,
    oxygen::FTO,
) where {FTR<:Real,FTT<:Real,FTC<:Real,FTO<:Real}
    require_compatible(dark_balanced_psii_partition, fl, p, temperature, co2, oxygen)
    rates = biochemical_rates(p.temperature, temperature, oxygen)
    phi1_max = psi_maximum_photochemical_yield(fl)
    eta = psi_to_psii_electron_ratio(rates, co2, p.stoichiometry)
    phi2_max = psii_maximum_photochemical_yield(fl, one(FTR))
    beta = phi1_max / (phi1_max + eta * phi2_max)
    T = typeof(beta + rates.vcmax + rates.vqmax + rates.gamma_star + temperature + co2 + oxygen)
    valid = valid_partition_inputs(fl, p, rates, phi1_max, temperature, co2, oxygen)
    return valid ? T(beta) : T(NaN)
end

@inline function dark_balanced_psii_partition(
    fl::CytochromeFluorescence,
    p::SerialCytochromePhotosynthesis,
    site::LeafEnvironment,
)
    val = dark_balanced_psii_partition(fl, p, site.temperature, site.biochemical_co2, site.oxygen)
    return is_nonnegative(site.absorbed_par) ? val : oftype(val, NaN)
end

"""
    state_transition_partition(rule::StateTransitionRule,
        fluorescence_model::CytochromeFluorescence,
        photosynthesis_or_leaf, site::LeafEnvironment)

Evaluate the pointwise, closed-form PSII absorbed-photon fraction `β ∈ (0, 1)`
under `rule` at `site` (Johnson & Berry 2021, Cases 10–12), returning `NaN` on
invalid or incompatible inputs without throwing. As for
[`dark_balanced_psii_partition`](@ref), the model must use [`SerialLimitation`](@ref) in its
electron transport; another co-limitation is a `MethodError`.
"""
@inline function state_transition_partition(
    rule::StateTransitionRule{FTR},
    fl::CytochromeFluorescence{FTF},
    p::SerialCytochromePhotosynthesis,
    site::LeafEnvironment{FT},
) where {FTR<:Real,FTF<:Real,FT<:Real}
    require_compatible(state_transition_partition, rule, fl, p, site)
    rates = biochemical_rates(p.temperature, site.temperature, site.oxygen)
    phi1_max = psi_maximum_photochemical_yield(fl)
    valid = valid_partition_inputs(fl, p, rates, phi1_max, site.temperature, site.biochemical_co2, site.oxygen) &&
        is_nonnegative(site.absorbed_par)
    T = typeof(
        zero(FTR) + zero(FTF) + rates.vcmax + rates.vqmax + rates.gamma_star +
        site.temperature + site.absorbed_par + site.biochemical_co2 + site.oxygen,
    )

    # Dimensionless PSII rate ratios bounded in [0, 1] (zero rate^2 overflow in Float32).
    psii = fl.psii
    d0 = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
    s_base = d0 + psii.exciton_transfer_rate
    s_tot = s_base + psii.maximum_photochemical_rate
    d0_hat = d0 / s_tot
    ku_hat = psii.exciton_transfer_rate / s_tot
    kp_hat = psii.maximum_photochemical_rate / s_tot
    sb_hat = s_base / s_tot

    eta = psi_to_psii_electron_ratio(rates, site.biochemical_co2, p.stoichiometry)
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
        disc = exceeds(rad2, zero(rad2)) ? sqrt(rad2) : zero(rad2)
        is_at_most(zero(b_coef), b_coef) ? (2 * c_pos) / (b_coef + disc) : (-b_coef + disc) / (2 * a_coef)
    end
    beta_opt = one(T) - alpha1_opt

    # Closed-form light-saturation floor β_sat at J_P700,j = J_P700,c.
    c1, c2 = cytochrome_carbon_coefficients(rates, p.stoichiometry)
    km = site.biochemical_co2 + rates.kc * (one(T) + site.oxygen / rates.ko)
    j_psi_c = (rates.vcmax / km) * (c1 * site.biochemical_co2 + c2 * rates.gamma_star)
    q2_sat = exceeds(rates.vqmax, zero(rates.vqmax)) ?
        primal_max(rates.vqmax - j_psi_c, zero(T)) / rates.vqmax : zero(T)
    num_sat = m1 + q2_sat * mu
    beta_sat = num_sat / (num_sat + m2)

    beta_raw = primal_max(beta_opt, beta_sat)
    R = promote_type(T, typeof(beta_raw))
    valid || return R(NaN)
    return R(primal_clamp(beta_raw, R(rule.minimum_psii_fraction), R(rule.maximum_psii_fraction)))
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
    p::SerialCytochromePhotosynthesis,
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
    if is_at_most(half, beta)
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
