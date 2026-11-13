"""
    PAMFluorescenceObserver(; psii_sensitivity, psi_sensitivity)
    PAMFluorescenceObserver(psii_sensitivity, psi_sensitivity)

Specify pulse-amplitude modulated (PAM) fluorometer detector weights
`(psii_sensitivity, psi_sensitivity)` (`ε_F2 > 0`, `ε_F1 ≥ 0`) mapping internal
PSII and PSI fluorescence yields to apparent detector levels `F / S`
(Johnson & Berry 2021, Eqs. 38–42).
"""
struct PAMFluorescenceObserver{FT<:Real}
    psii_sensitivity::FT
    psi_sensitivity::FT
    function PAMFluorescenceObserver(
        psii_sensitivity::FT,
        psi_sensitivity::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT, PAMFluorescenceObserver, (psii_sensitivity, psi_sensitivity),
        )
        e2p = primal_scalar(psii_sensitivity)
        e1p = primal_scalar(psi_sensitivity)
        isfinite(e2p) && isfinite(e1p) &&
            e2p > zero(e2p) && e1p >= zero(e1p) ||
            throw(ArgumentError("PAMFluorescenceObserver requires finite psii_sensitivity > 0 and psi_sensitivity ≥ 0"))
        return new{FT}(psii_sensitivity, psi_sensitivity)
    end
end

PAMFluorescenceObserver(psii_sensitivity::Real, psi_sensitivity::Real) =
    PAMFluorescenceObserver(
        promote_compatible(PAMFluorescenceObserver, psii_sensitivity, psi_sensitivity)...,
    )

PAMFluorescenceObserver(; psii_sensitivity, psi_sensitivity) =
    PAMFluorescenceObserver(psii_sensitivity, psi_sensitivity)

"""
    PAMFluorescenceResult{FT}

Store apparent PAM fluorometer levels (`fs`, `fm_dark`, `fo_dark`, `fm_prime`,
`fo_prime` in `F / S` units; Johnson & Berry 2021, Eqs. 38–42) and derived
pulse-ratio indices (`apparent_photochemical_yield`, `apparent_regulated_heat_yield`,
`apparent_residual_yield`, `apparent_npq`, `qp`, `ql`).

Check [`pam_fluorescence_valid`](@ref) before consuming fields; non-lit states
carry `NaN` fields without throwing on GPU devices.
"""
struct PAMFluorescenceResult{FT<:Real}
    status::FluorescenceStatus
    fs::FT
    fm_dark::FT
    fo_dark::FT
    fm_prime::FT
    fo_prime::FT
    apparent_photochemical_yield::FT
    apparent_regulated_heat_yield::FT
    apparent_residual_yield::FT
    apparent_npq::FT
    qp::FT
    ql::FT
    function PAMFluorescenceResult{FT}(
        status::FluorescenceStatus,
        fs,
        fm_dark,
        fo_dark,
        fm_prime,
        fo_prime,
        apparent_photochemical_yield,
        apparent_regulated_heat_yield,
        apparent_residual_yield,
        apparent_npq,
        qp,
        ql,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            PAMFluorescenceResult{FT},
            (
                status,
                fs,
                fm_dark,
                fo_dark,
                fm_prime,
                fo_prime,
                apparent_photochemical_yield,
                apparent_regulated_heat_yield,
                apparent_residual_yield,
                apparent_npq,
                qp,
                ql,
            ),
        )
        return new{FT}(
            status,
            fs,
            fm_dark,
            fo_dark,
            fm_prime,
            fo_prime,
            apparent_photochemical_yield,
            apparent_regulated_heat_yield,
            apparent_residual_yield,
            apparent_npq,
            qp,
            ql,
        )
    end
    function PAMFluorescenceResult(
        status::FluorescenceStatus,
        fs::FT,
        fm_dark::FT,
        fo_dark::FT,
        fm_prime::FT,
        fo_prime::FT,
        apparent_photochemical_yield::FT,
        apparent_regulated_heat_yield::FT,
        apparent_residual_yield::FT,
        apparent_npq::FT,
        qp::FT,
        ql::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            PAMFluorescenceResult,
            (
                status,
                fs,
                fm_dark,
                fo_dark,
                fm_prime,
                fo_prime,
                apparent_photochemical_yield,
                apparent_regulated_heat_yield,
                apparent_residual_yield,
                apparent_npq,
                qp,
                ql,
            ),
        )
        return new{FT}(
            status,
            fs,
            fm_dark,
            fo_dark,
            fm_prime,
            fo_prime,
            apparent_photochemical_yield,
            apparent_regulated_heat_yield,
            apparent_residual_yield,
            apparent_npq,
            qp,
            ql,
        )
    end
end

"""
    pam_fluorescence_valid(result::PAMFluorescenceResult) -> Bool

Return `true` if `result.status` is [`FluorescenceLit`](@ref).
"""
@inline pam_fluorescence_valid(r::PAMFluorescenceResult) = r.status == FluorescenceLit

@inline function unavailable_pam_fluorescence(::Type{FT}, status::FluorescenceStatus) where {FT}
    nan = FT(NaN)
    return PAMFluorescenceResult(
        status,
        nan, nan, nan, nan, nan,
        nan, nan, nan, nan, nan, nan,
    )
end

"""
    dark_pam_levels(observer::PAMFluorescenceObserver,
        fluorescence_model::CytochromeFluorescence,
        alpha_psii_dark, alpha_psi_dark)

Compute dark-adapted (`K_N = 0`) PAM reference levels `(; fo, fm, fv_fm)`
(`F_o / S`, `F_m / S`, and `(F_m - F_o) / F_m`; Johnson & Berry 2021, Eqs. 39, 41),
returning `NaN` fields on invalid cross-sections without throwing.
"""
@inline function dark_pam_levels(
    observer::PAMFluorescenceObserver{FTO},
    fl::CytochromeFluorescence{FTF},
    alpha_psii_dark::FT2,
    alpha_psi_dark::FT1,
) where {FTO<:Real,FTF<:Real,FT2<:Real,FT1<:Real}
    _, _, a2_d, a1_d = promote_compatible(
        dark_pam_levels, zero(FTO), zero(FTF), alpha_psii_dark, alpha_psi_dark,
    )
    S = typeof(a2_d)
    nan = S(NaN)
    a2d_p, a1d_p = primal_scalar(a2_d), primal_scalar(a1_d)
    isfinite(a2d_p) && isfinite(a1d_p) &&
        a2d_p > zero(a2d_p) && a1d_p >= zero(a1d_p) ||
        return (; fo = nan, fm = nan, fv_fm = nan)

    psii = fl.psii
    psi = fl.psi
    eps2 = S(observer.psii_sensitivity)
    eps1 = S(observer.psi_sensitivity)

    s2_closed = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
    s2_open = s2_closed + psii.maximum_photochemical_rate
    s1_open = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.maximum_photochemical_rate
    s1_closed = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.closed_quenching_rate

    phi_f2_fm = S(psii.fluorescence_rate / s2_closed)
    phi_f2_fo = S(psii.fluorescence_rate / s2_open)
    phi_f1_closed = S(psi.fluorescence_rate / s1_closed)
    phi_f1_open = S(psi.fluorescence_rate / s1_open)

    fm = S((a2_d * phi_f2_fm) * eps2 + (a1_d * phi_f1_closed) * eps1)
    fo = S((a2_d * phi_f2_fo) * eps2 + (a1_d * phi_f1_open) * eps1)
    psi_diff = S(phi_f1_closed * ((psi.maximum_photochemical_rate - psi.closed_quenching_rate) / s1_open))
    fv = S(
        (a2_d * (phi_f2_fm * S(psii.maximum_photochemical_rate / s2_open))) * eps2 +
        (a1_d * psi_diff) * eps1,
    )
    fv_fm = S(fv / fm)
    all(isfinite, (fo, fm, fv_fm)) || return (; fo = nan, fm = nan, fv_fm = nan)
    return (; fo, fm, fv_fm)
end

"""
    pam_fluorescence(observer::PAMFluorescenceObserver,
        fluorescence_model::CytochromeFluorescence,
        result::CytochromeFluorescenceResult,
        alpha_psii, alpha_psi[, alpha_psii_dark, alpha_psi_dark])

Evaluate apparent PAM fluorometer observables (Johnson & Berry 2021, Eqs. 38–42)
from coupled internal fluorescence `result`, light-adapted absorption
cross-sections `(alpha_psii, alpha_psi)`, and dark-reference cross-sections
`(alpha_psii_dark, alpha_psi_dark)` (defaulting to `(alpha_psii, alpha_psi)`).
"""
@inline function pam_fluorescence(
    observer::PAMFluorescenceObserver{FTO},
    fl::CytochromeFluorescence{FTF},
    result::CytochromeFluorescenceResult{FTR},
    alpha_psii::FT2,
    alpha_psi::FT1,
    alpha_psii_dark::FT2D = alpha_psii,
    alpha_psi_dark::FT1D = alpha_psi,
) where {FTO<:Real,FTF<:Real,FTR<:Real,FT2<:Real,FT1<:Real,FT2D<:Real,FT1D<:Real}
    _, _, _, a2, a1, a2_d, a1_d = promote_compatible(
        pam_fluorescence,
        zero(FTO), zero(FTF), zero(FTR),
        alpha_psii, alpha_psi, alpha_psii_dark, alpha_psi_dark,
    )
    S = typeof(a2)
    a2_p, a1_p, a2d_p, a1d_p = map(primal_scalar, (a2, a1, a2_d, a1_d))
    isfinite(a2_p) && isfinite(a1_p) &&
        isfinite(a2d_p) && isfinite(a1d_p) &&
        a2_p > zero(a2_p) && a1_p >= zero(a1_p) &&
        a2d_p > zero(a2d_p) && a1d_p >= zero(a1d_p) ||
        return unavailable_pam_fluorescence(S, FluorescenceInvalidInput)

    result.status == FluorescenceLit ||
        return unavailable_pam_fluorescence(S, result.status)

    psii = fl.psii
    psi = fl.psi
    eps2 = S(observer.psii_sensitivity)
    eps1 = S(observer.psi_sensitivity)
    kn = S(result.psii_quenching_rate)

    s2_closed_dark = S(psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate)
    s2_open_dark = S(s2_closed_dark + psii.maximum_photochemical_rate)
    s2_closed_light = s2_closed_dark + kn
    s_scale = max(s2_closed_light, S(psii.maximum_photochemical_rate))
    denom_open_light = (s2_closed_light / s_scale) + (S(psii.maximum_photochemical_rate) / s_scale)

    s1_open = S(psi.fluorescence_rate + psi.constitutive_heat_rate + psi.maximum_photochemical_rate)
    s1_closed = S(psi.fluorescence_rate + psi.constitutive_heat_rate + psi.closed_quenching_rate)

    phi_f2_fm_dark = S(psii.fluorescence_rate) / s2_closed_dark
    phi_f2_fo_dark = S(psii.fluorescence_rate) / s2_open_dark
    phi_f2_fmp = S(psii.fluorescence_rate) / s2_closed_light
    phi_f2_fop = (S(psii.fluorescence_rate) / s_scale) / denom_open_light
    kp2_frac_light = (S(psii.maximum_photochemical_rate) / s_scale) / denom_open_light
    phi_f1_closed = S(psi.fluorescence_rate) / s1_closed
    phi_f1_open = S(psi.fluorescence_rate) / s1_open

    fm_dark = (a2_d * phi_f2_fm_dark) * eps2 + (a1_d * phi_f1_closed) * eps1
    fo_dark = (a2_d * phi_f2_fo_dark) * eps2 + (a1_d * phi_f1_open) * eps1
    fm_prime = (a2 * phi_f2_fmp) * eps2 + (a1 * phi_f1_closed) * eps1
    fo_prime = (a2 * phi_f2_fop) * eps2 + (a1 * phi_f1_open) * eps1
    fs = (a2 * S(result.psii_fluorescence_yield)) * eps2 +
         (a1 * S(result.psi_fluorescence_yield)) * eps1

    psi_closed_minus_open = phi_f1_closed *
        ((S(psi.maximum_photochemical_rate) - S(psi.closed_quenching_rate)) / s1_open)
    fmp_minus_fop = (a2 * (phi_f2_fmp * kp2_frac_light)) * eps2 +
        (a1 * psi_closed_minus_open) * eps1
    fmp_minus_fs = (a2 * (S(result.psii_photochemical_yield) * phi_f2_fmp)) * eps2 +
        (a1 * (S(result.psi_open_fraction) * psi_closed_minus_open)) * eps1
    fm_minus_fmp = ((a2_d - a2) * phi_f2_fm_dark +
        a2 * (phi_f2_fm_dark * (kn / s2_closed_light))) * eps2 +
        ((a1_d - a1) * phi_f1_closed) * eps1

    phi_p_app = fmp_minus_fs / fm_prime
    phi_df_app = fs / fm_dark
    npq_app = fm_minus_fmp / fm_prime
    phi_n_app = phi_df_app * npq_app
    qp = fmp_minus_fs / fmp_minus_fop
    ql = qp * (fo_prime / fs)

    all(isfinite, (fs, fm_dark, fo_dark, fm_prime, fo_prime,
                   phi_p_app, phi_n_app, phi_df_app, npq_app, qp, ql)) ||
        return unavailable_pam_fluorescence(S, FluorescenceIncompatibleRates)

    return PAMFluorescenceResult{S}(
        FluorescenceLit,
        fs,
        fm_dark,
        fo_dark,
        fm_prime,
        fo_prime,
        phi_p_app,
        phi_n_app,
        phi_df_app,
        npq_app,
        qp,
        ql,
    )
end
