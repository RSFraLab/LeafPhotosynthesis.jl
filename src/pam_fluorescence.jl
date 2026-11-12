"""
    PAMFluorescenceObserver(; psii_sensitivity, psi_sensitivity)

Pulse-amplitude modulated (PAM) fluorometer observation weights `(ε_F2, ε_F1)`
mapping internal PSII and PSI fluorescence yields to apparent detector levels
`F / S` (Johnson & Berry 2021, Eqs. 38–42):
- `psii_sensitivity` (`ε_F2 > 0`): detector weight for PSII fluorescence
- `psi_sensitivity` (`ε_F1 ≥ 0`): detector weight for PSI fluorescence (`0` for an
  ideal PSII-only detector; `ε_F1 / ε_F2 ≈ 1.6–2.5` for standard long-pass >700 nm
  leaf-chamber fluorometers where PSI contributes ~39% of `F_o` and ~7% of `F_m`)
"""
struct PAMFluorescenceObserver{FT<:AbstractFloat}
    psii_sensitivity::FT
    psi_sensitivity::FT
    function PAMFluorescenceObserver(
        psii_sensitivity::FT,
        psi_sensitivity::FT,
    ) where {FT<:AbstractFloat}
        isfinite(psii_sensitivity) && isfinite(psi_sensitivity) &&
            psii_sensitivity > zero(FT) && psi_sensitivity >= zero(FT) ||
            throw(ArgumentError("PAMFluorescenceObserver requires finite psii_sensitivity > 0 and psi_sensitivity ≥ 0"))
        new{FT}(psii_sensitivity, psi_sensitivity)
    end
end

PAMFluorescenceObserver(; psii_sensitivity, psi_sensitivity) =
    PAMFluorescenceObserver(psii_sensitivity, psi_sensitivity)

"""
    PAMFluorescenceResult{FT}

Apparent PAM fluorometer levels (`F / S`, Johnson & Berry 2021, Eqs. 38–42) and
derived pulse-ratio indices:
- `status`: `FluorescenceStatus` (`FluorescenceLit` when populated; non-lit states
  carry `NaN` fields without throwing on GPU)
- `fs`, `fm_dark`, `fo_dark`, `fm_prime`, `fo_prime`: detector-weighted steady-state,
  dark-reference, and light-adapted fluorescence levels (`F / S`)
- `apparent_photochemical_yield`: Genty et al. (1989) `1 - F_s / F_m'`
- `apparent_regulated_heat_yield`: Hendrickson et al. (2004) `F_s * (1/F_m' - 1/F_m)`
- `apparent_residual_yield`: Hendrickson et al. (2004) `F_s / F_m`
- `apparent_npq`: Bilger & Björkman (1990) `F_m / F_m' - 1`
- `qp`: Schreiber et al. (1986) puddle-model open fraction `(F_m' - F_s) / (F_m' - F_o')`
- `ql`: Kramer et al. (2004) lake-model open fraction `q_P * (F_o' / F_s)`
"""
struct PAMFluorescenceResult{FT}
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
end

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

Dark-adapted (`K_N = 0`) PAM reference levels `(; fo, fm, fv_fm)` (`F_o / S`,
`F_m / S`, and `(F_m - F_o) / F_m`; Johnson & Berry 2021, Eqs. 39, 41). Nonthrowing
on GPU; returns `NaN` fields for invalid cross-sections.
"""
@inline function dark_pam_levels(
    observer::PAMFluorescenceObserver{FT},
    fl::CytochromeFluorescence{FT},
    alpha_psii_dark::FT,
    alpha_psi_dark::FT,
) where {FT<:AbstractFloat}
    nan = FT(NaN)
    isfinite(alpha_psii_dark) && isfinite(alpha_psi_dark) &&
        alpha_psii_dark > zero(FT) && alpha_psi_dark >= zero(FT) ||
        return (; fo = nan, fm = nan, fv_fm = nan)

    psii = fl.psii
    psi = fl.psi
    eps2 = observer.psii_sensitivity
    eps1 = observer.psi_sensitivity

    s2_closed = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
    s2_open = s2_closed + psii.maximum_photochemical_rate
    s1_open = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.maximum_photochemical_rate
    s1_closed = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.closed_quenching_rate

    phi_f2_fm = psii.fluorescence_rate / s2_closed
    phi_f2_fo = psii.fluorescence_rate / s2_open
    phi_f1_closed = psi.fluorescence_rate / s1_closed
    phi_f1_open = psi.fluorescence_rate / s1_open

    fm = (alpha_psii_dark * phi_f2_fm) * eps2 + (alpha_psi_dark * phi_f1_closed) * eps1
    fo = (alpha_psii_dark * phi_f2_fo) * eps2 + (alpha_psi_dark * phi_f1_open) * eps1
    psi_diff = phi_f1_closed * ((psi.maximum_photochemical_rate - psi.closed_quenching_rate) / s1_open)
    fv = (alpha_psii_dark * (phi_f2_fm * (psii.maximum_photochemical_rate / s2_open))) * eps2 +
         (alpha_psi_dark * psi_diff) * eps1
    fv_fm = fv / fm
    all(isfinite, (fo, fm, fv_fm)) || return (; fo = nan, fm = nan, fv_fm = nan)
    return (; fo, fm, fv_fm)
end

"""
    pam_fluorescence(observer::PAMFluorescenceObserver,
        fluorescence_model::CytochromeFluorescence,
        result::CytochromeFluorescenceResult,
        alpha_psii, alpha_psi[, alpha_psii_dark, alpha_psi_dark])

Evaluate apparent PAM fluorometer observables (Johnson & Berry 2021, Eqs. 38–42)
from coupled internal fluorescence `result`, light-adapted absorption cross-sections
`(alpha_psii, alpha_psi)`, and dark-reference cross-sections
`(alpha_psii_dark, alpha_psi_dark)` (defaulting to `(alpha_psii, alpha_psi)` when
pigment cross-sections are static).

Ratio-based indices (`apparent_photochemical_yield`, `apparent_regulated_heat_yield`,
`apparent_residual_yield`, `apparent_npq`, `qp`, `ql`) are invariant to any common
positive scaling of `(alpha_psii, alpha_psi, alpha_psii_dark, alpha_psi_dark)` (such
as constant leaf PAR absorptance). Normalized rate ratios and factored differences
avoid both `rate²` overflow and `Float32` PSI-offset cancellation in `F_m' - F_o'`,
`F_m' - F_s`, and `F_m - F_m'`.
"""
@inline function pam_fluorescence(
    observer::PAMFluorescenceObserver{FT},
    fl::CytochromeFluorescence{FT},
    result::CytochromeFluorescenceResult{FT},
    alpha_psii::FT,
    alpha_psi::FT,
    alpha_psii_dark::FT = alpha_psii,
    alpha_psi_dark::FT = alpha_psi,
) where {FT<:AbstractFloat}
    isfinite(alpha_psii) && isfinite(alpha_psi) &&
        isfinite(alpha_psii_dark) && isfinite(alpha_psi_dark) &&
        alpha_psii > zero(FT) && alpha_psi >= zero(FT) &&
        alpha_psii_dark > zero(FT) && alpha_psi_dark >= zero(FT) ||
        return unavailable_pam_fluorescence(FT, FluorescenceInvalidInput)

    result.status == FluorescenceLit ||
        return unavailable_pam_fluorescence(FT, result.status)

    psii = fl.psii
    psi = fl.psi
    eps2 = observer.psii_sensitivity
    eps1 = observer.psi_sensitivity
    kn = result.psii_quenching_rate

    s2_closed_dark = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
    s2_open_dark = s2_closed_dark + psii.maximum_photochemical_rate
    s2_closed_light = s2_closed_dark + kn
    s_scale = max(s2_closed_light, psii.maximum_photochemical_rate)
    denom_open_light = (s2_closed_light / s_scale) + (psii.maximum_photochemical_rate / s_scale)

    s1_open = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.maximum_photochemical_rate
    s1_closed = psi.fluorescence_rate + psi.constitutive_heat_rate + psi.closed_quenching_rate

    phi_f2_fm_dark = psii.fluorescence_rate / s2_closed_dark
    phi_f2_fo_dark = psii.fluorescence_rate / s2_open_dark
    phi_f2_fmp = psii.fluorescence_rate / s2_closed_light
    phi_f2_fop = (psii.fluorescence_rate / s_scale) / denom_open_light
    kp2_frac_light = (psii.maximum_photochemical_rate / s_scale) / denom_open_light
    phi_f1_closed = psi.fluorescence_rate / s1_closed
    phi_f1_open = psi.fluorescence_rate / s1_open

    fm_dark = (alpha_psii_dark * phi_f2_fm_dark) * eps2 + (alpha_psi_dark * phi_f1_closed) * eps1
    fo_dark = (alpha_psii_dark * phi_f2_fo_dark) * eps2 + (alpha_psi_dark * phi_f1_open) * eps1
    fm_prime = (alpha_psii * phi_f2_fmp) * eps2 + (alpha_psi * phi_f1_closed) * eps1
    fo_prime = (alpha_psii * phi_f2_fop) * eps2 + (alpha_psi * phi_f1_open) * eps1
    fs = (alpha_psii * result.psii_fluorescence_yield) * eps2 +
         (alpha_psi * result.psi_fluorescence_yield) * eps1

    # Overflow- and cancellation-free differences:
    # Φ_F2,m' - Φ_F2 = Φ_P2 * Φ_F2,m' avoids 1 - Φ_P2 cancellation at low Φ_P2,
    # and normalized ratios avoid rate^2 and (alpha * eps)^2 overflow/underflow in Float32.
    psi_closed_minus_open = phi_f1_closed *
        ((psi.maximum_photochemical_rate - psi.closed_quenching_rate) / s1_open)
    fmp_minus_fop = (alpha_psii * (phi_f2_fmp * kp2_frac_light)) * eps2 +
        (alpha_psi * psi_closed_minus_open) * eps1
    fmp_minus_fs = (alpha_psii * (result.psii_photochemical_yield * phi_f2_fmp)) * eps2 +
        (alpha_psi * (result.psi_open_fraction * psi_closed_minus_open)) * eps1
    fm_minus_fmp = ((alpha_psii_dark - alpha_psii) * phi_f2_fm_dark +
        alpha_psii * (phi_f2_fm_dark * (kn / s2_closed_light))) * eps2 +
        ((alpha_psi_dark - alpha_psi) * phi_f1_closed) * eps1

    phi_p_app = fmp_minus_fs / fm_prime
    phi_df_app = fs / fm_dark
    npq_app = fm_minus_fmp / fm_prime
    phi_n_app = phi_df_app * npq_app
    qp = fmp_minus_fs / fmp_minus_fop
    ql = qp * (fo_prime / fs)

    all(isfinite, (fs, fm_dark, fo_dark, fm_prime, fo_prime,
                   phi_p_app, phi_n_app, phi_df_app, npq_app, qp, ql)) ||
        return unavailable_pam_fluorescence(FT, FluorescenceIncompatibleRates)

    return PAMFluorescenceResult(
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
