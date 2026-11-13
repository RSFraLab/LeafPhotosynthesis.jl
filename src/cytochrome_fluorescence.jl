"""
    PSIIRateParameters(; fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)
    PSIIRateParameters(fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)

Specify first-order PSII excitation decay and inter-unit transfer rate constants
sharing a common normalization (conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F2 > 0`).
- `constitutive_heat_rate` (`K_D2 ≥ 0`).
- `maximum_photochemical_rate` (`K_P2 > 0`).
- `exciton_transfer_rate` (`K_U2 ≥ 0`; `0` selects the separate-units puddle limit,
  `> 0` selects connected-units exciton sharing).
- `sustained_quenching_rate` (`K_N,sus ≥ 0`; dark-sustained thermal dissipation).
"""
struct PSIIRateParameters{FT<:Real}
    fluorescence_rate::FT
    constitutive_heat_rate::FT
    maximum_photochemical_rate::FT
    exciton_transfer_rate::FT
    sustained_quenching_rate::FT
    function PSIIRateParameters(
        fluorescence_rate::FT,
        constitutive_heat_rate::FT,
        maximum_photochemical_rate::FT,
        exciton_transfer_rate::FT,
        sustained_quenching_rate::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            PSIIRateParameters,
            (
                fluorescence_rate,
                constitutive_heat_rate,
                maximum_photochemical_rate,
                exciton_transfer_rate,
                sustained_quenching_rate,
            ),
        )
        kf, kd, kp, ku, kns = map(
            primal_scalar,
            (
                fluorescence_rate,
                constitutive_heat_rate,
                maximum_photochemical_rate,
                exciton_transfer_rate,
                sustained_quenching_rate,
            ),
        )
        d0 = kf + kd + kns
        s_base = d0 + ku
        total_open = d0 + kp
        phi_open = kp / total_open
        d0_hat = d0 / s_base
        ku_hat = ku / s_base
        kp_hat = kp / s_base
        phi_conn = kp_hat / (d0_hat + kp_hat * (d0_hat + ku_hat))
        all(isfinite, (kf, kd, kp, ku, kns, d0, s_base, total_open, s_base + kp, phi_conn)) &&
            kf > zero(kf) && kd >= zero(kd) && kp > zero(kp) &&
            ku >= zero(ku) && kns >= zero(kns) &&
            zero(phi_open) < phi_open < one(phi_open) && zero(phi_conn) < phi_conn < one(phi_conn) ||
            throw(ArgumentError("Invalid or unrepresentable PSII rate parameters"))
        return new{FT}(
            fluorescence_rate,
            constitutive_heat_rate,
            maximum_photochemical_rate,
            exciton_transfer_rate,
            sustained_quenching_rate,
        )
    end
end
PSIIRateParameters(kf::Real, kd::Real, kp::Real, ku::Real, kns::Real) =
    PSIIRateParameters(promote_compatible(PSIIRateParameters, kf, kd, kp, ku, kns)...)
PSIIRateParameters(;
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    exciton_transfer_rate,
    sustained_quenching_rate,
) = PSIIRateParameters(
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    exciton_transfer_rate,
    sustained_quenching_rate,
)

"""
    PSIRateParameters(; fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, closed_quenching_rate)
    PSIRateParameters(fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, closed_quenching_rate)

Specify first-order PSI excitation decay rate constants sharing one normalization
(conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F1 > 0`).
- `constitutive_heat_rate` (`K_D1 ≥ 0`).
- `maximum_photochemical_rate` (`K_P1 > 0` in open P700 centers).
- `closed_quenching_rate` (`K_X1 > 0` radiationless decay in oxidized P700⁺ centers).
"""
struct PSIRateParameters{FT<:Real}
    fluorescence_rate::FT
    constitutive_heat_rate::FT
    maximum_photochemical_rate::FT
    closed_quenching_rate::FT
    function PSIRateParameters(
        fluorescence_rate::FT,
        constitutive_heat_rate::FT,
        maximum_photochemical_rate::FT,
        closed_quenching_rate::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            PSIRateParameters,
            (
                fluorescence_rate,
                constitutive_heat_rate,
                maximum_photochemical_rate,
                closed_quenching_rate,
            ),
        )
        kf, kd, kp, kx = map(
            primal_scalar,
            (
                fluorescence_rate,
                constitutive_heat_rate,
                maximum_photochemical_rate,
                closed_quenching_rate,
            ),
        )
        s_open = kf + kd + kp
        s_closed = kf + kd + kx
        phi_open = kp / s_open
        phi_closed = kx / s_closed
        all(isfinite, (kf, kd, kp, kx, s_open, s_closed)) &&
            kf > zero(kf) && kd >= zero(kd) && kp > zero(kp) && kx > zero(kx) &&
            zero(phi_open) < phi_open < one(phi_open) && zero(phi_closed) < phi_closed < one(phi_closed) ||
            throw(ArgumentError("Invalid or unrepresentable PSI rate parameters"))
        return new{FT}(
            fluorescence_rate,
            constitutive_heat_rate,
            maximum_photochemical_rate,
            closed_quenching_rate,
        )
    end
end
PSIRateParameters(kf::Real, kd::Real, kp::Real, kx::Real) =
    PSIRateParameters(promote_compatible(PSIRateParameters, kf, kd, kp, kx)...)
PSIRateParameters(;
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    closed_quenching_rate,
) = PSIRateParameters(
    fluorescence_rate,
    constitutive_heat_rate,
    maximum_photochemical_rate,
    closed_quenching_rate,
)

@inline _convert_psii(::Type{FT}, p::PSIIRateParameters{FT}) where {FT<:Real} = p
@inline _convert_psii(::Type{FT}, p::PSIIRateParameters) where {FT<:Real} =
    PSIIRateParameters(
        FT(p.fluorescence_rate),
        FT(p.constitutive_heat_rate),
        FT(p.maximum_photochemical_rate),
        FT(p.exciton_transfer_rate),
        FT(p.sustained_quenching_rate),
    )

@inline _convert_psi(::Type{FT}, p::PSIRateParameters{FT}) where {FT<:Real} = p
@inline _convert_psi(::Type{FT}, p::PSIRateParameters) where {FT<:Real} =
    PSIRateParameters(
        FT(p.fluorescence_rate),
        FT(p.constitutive_heat_rate),
        FT(p.maximum_photochemical_rate),
        FT(p.closed_quenching_rate),
    )

"""
    CytochromeFluorescence(; psii, psi, yield_tolerance)
    CytochromeFluorescence(psii, psi, yield_tolerance)

Parameterize the mechanistic steady-state PSII + PSI fluorescence model coupled
through Cytochrome b₆f redox poise and energy-dependent feedback quenching `K_N`
(Johnson & Berry 2021, Eqs. 7–30).
"""
struct CytochromeFluorescence{FT<:Real} <: AbstractFluorescence
    psii::PSIIRateParameters{FT}
    psi::PSIRateParameters{FT}
    yield_tolerance::FT
    function CytochromeFluorescence(
        psii::PSIIRateParameters{FT},
        psi::PSIRateParameters{FT},
        yield_tolerance::FT,
    ) where {FT<:Real}
        require_continuous_scalar(FT, CytochromeFluorescence, (psii, psi, yield_tolerance))
        phi2_max = psii_maximum_photochemical_yield(psii, one(FT))
        phi1_max = psi_maximum_photochemical_yield(psi)
        headroom = primal_scalar(min(one(FT) - phi2_max, one(FT) - phi1_max))
        ytp = primal_scalar(yield_tolerance)
        isfinite(ytp) && zero(ytp) <= ytp < headroom ||
            throw(ArgumentError("Invalid Cytochrome fluorescence yield tolerance"))
        return new{FT}(psii, psi, yield_tolerance)
    end
end
function CytochromeFluorescence(
    psii::PSIIRateParameters{FT2},
    psi::PSIRateParameters{FT1},
    yield_tolerance::Real,
) where {FT2<:Real,FT1<:Real}
    _, _, yt = promote_compatible(CytochromeFluorescence, zero(FT2), zero(FT1), yield_tolerance)
    FT = typeof(yt)
    return CytochromeFluorescence(_convert_psii(FT, psii), _convert_psi(FT, psi), yt)
end
CytochromeFluorescence(; psii, psi, yield_tolerance) =
    CytochromeFluorescence(psii, psi, yield_tolerance)

"""
    psi_maximum_photochemical_yield(psi_or_model)

Compute the open-center PSI maximum photochemical quantum yield
``K_{P1} / (K_{F1} + K_{D1} + K_{P1})``.
"""
@inline psi_maximum_photochemical_yield(p::PSIRateParameters) =
    p.maximum_photochemical_rate /
    (p.fluorescence_rate + p.constitutive_heat_rate + p.maximum_photochemical_rate)
@inline psi_maximum_photochemical_yield(m::CytochromeFluorescence) =
    psi_maximum_photochemical_yield(m.psi)

"""
    psii_maximum_photochemical_yield(psii_or_model[, psii_open_fraction = 1])

Compute the PSII photochemical quantum yield at zero dynamic feedback quenching
(`K_N = 0`) and open reaction-center fraction `psii_open_fraction` (`∈ [0, 1]`),
accounting for connected-units exciton sharing (`K_U2`) and sustained quenching (`K_N,sus`).
"""
@inline function psii_maximum_photochemical_yield(
    p::PSIIRateParameters{FT},
    psii_open_fraction::FTQ = one(FT),
) where {FT<:Real,FTQ<:Real}
    compatible_scalar(FT, FTQ) ||
        throw(MethodError(psii_maximum_photochemical_yield, (p, psii_open_fraction)))
    S = promote_type(FT, FTQ)
    d0 = p.fluorescence_rate + p.constitutive_heat_rate + p.sustained_quenching_rate
    s_base = d0 + p.exciton_transfer_rate
    d0_hat = d0 / s_base
    ku_hat = p.exciton_transfer_rate / s_base
    kp_hat = p.maximum_photochemical_rate / s_base
    m0_hat = d0_hat + kp_hat * (d0_hat + psii_open_fraction * ku_hat)
    return S((psii_open_fraction * kp_hat) / m0_hat)
end
@inline psii_maximum_photochemical_yield(
    m::CytochromeFluorescence{FT},
    psii_open_fraction::Real = one(FT),
) where {FT<:Real} = psii_maximum_photochemical_yield(m.psii, psii_open_fraction)

"""
    CoupledPhotosystemInputs{FT}(psii_excitation, psi_excitation,
        linear_electron_flow, psi_electron_flow, psii_open_fraction)
    CoupledPhotosystemInputs(psii_excitation, psi_excitation,
        linear_electron_flow, psi_electron_flow, psii_open_fraction)
    CoupledPhotosystemInputs(; psii_excitation, psi_excitation,
        linear_electron_flow, psi_electron_flow, psii_open_fraction)

Bundle scalar inputs for [`CytochromeFluorescence`](@ref):
- `psii_excitation`, `psi_excitation`: absorbed photon fluxes delivered to PSII
  and PSI [μmol photons m⁻² one-sided leaf s⁻¹].
- `linear_electron_flow`: realized PSII linear electron flux `J_P680` [μmol e⁻ m⁻² leaf s⁻¹].
- `psi_electron_flow`: realized PSI electron flux `J_P700 = J_P680 + J_CEF` [μmol e⁻ m⁻² leaf s⁻¹].
- `psii_open_fraction`: open PSII reaction-center fraction `q₂ ∈ [0, 1]`.
"""
struct CoupledPhotosystemInputs{FT<:Real}
    psii_excitation::FT
    psi_excitation::FT
    linear_electron_flow::FT
    psi_electron_flow::FT
    psii_open_fraction::FT
    function CoupledPhotosystemInputs{FT}(
        q2::FT, q1::FT, j2::FT, j1::FT, q2_open::FT,
    ) where {FT<:Real}
        require_continuous_scalar(FT, CoupledPhotosystemInputs, (q2, q1, j2, j1, q2_open))
        return new{FT}(q2, q1, j2, j1, q2_open)
    end
end
CoupledPhotosystemInputs(q2::FT, q1::FT, j2::FT, j1::FT, q2_open::FT) where {FT<:Real} =
    CoupledPhotosystemInputs{FT}(q2, q1, j2, j1, q2_open)
CoupledPhotosystemInputs(q2::Real, q1::Real, j2::Real, j1::Real, q2_open::Real) =
    CoupledPhotosystemInputs(
        promote_compatible(CoupledPhotosystemInputs, q2, q1, j2, j1, q2_open)...,
    )
CoupledPhotosystemInputs(;
    psii_excitation,
    psi_excitation,
    linear_electron_flow,
    psi_electron_flow,
    psii_open_fraction,
) = CoupledPhotosystemInputs(
    psii_excitation,
    psi_excitation,
    linear_electron_flow,
    psi_electron_flow,
    psii_open_fraction,
)

"""
    CytochromeFluorescenceResult{FT}

Store scalar PSII and PSI fluorescence fluxes, quantum yields, open-center
fractions, dynamic quenching `psii_quenching_rate` (`K_N`), and internal
Stern–Volmer `npq` for [`CytochromeFluorescence`](@ref).

`psii_photon_flux`, `psi_photon_flux`, and `total_photon_flux` represent
internal leaf emission rates [μmol photons m⁻² one-sided leaf s⁻¹] prior to
within-leaf optical escape and reabsorption. Check [`fluorescence_valid`](@ref)
before consuming fluxes or yields.
"""
struct CytochromeFluorescenceResult{FT<:Real}
    status::FluorescenceStatus
    boundary_roundoff::Bool
    psii_photon_flux::FT
    psi_photon_flux::FT
    total_photon_flux::FT
    psii_photochemical_yield::FT
    psii_fluorescence_yield::FT
    psii_constitutive_heat_yield::FT
    psii_regulated_heat_yield::FT
    psi_photochemical_yield::FT
    psi_fluorescence_yield::FT
    psi_constitutive_heat_yield::FT
    psi_regulated_heat_yield::FT
    total_fluorescence_yield::FT
    psii_open_fraction::FT
    psi_open_fraction::FT
    psii_quenching_rate::FT
    npq::FT
    function CytochromeFluorescenceResult{FT}(
        status::FluorescenceStatus,
        boundary_roundoff::Bool,
        psii_photon_flux,
        psi_photon_flux,
        total_photon_flux,
        psii_photochemical_yield,
        psii_fluorescence_yield,
        psii_constitutive_heat_yield,
        psii_regulated_heat_yield,
        psi_photochemical_yield,
        psi_fluorescence_yield,
        psi_constitutive_heat_yield,
        psi_regulated_heat_yield,
        total_fluorescence_yield,
        psii_open_fraction,
        psi_open_fraction,
        psii_quenching_rate,
        npq,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeFluorescenceResult{FT},
            (
                status,
                boundary_roundoff,
                psii_photon_flux,
                psi_photon_flux,
                total_photon_flux,
                psii_photochemical_yield,
                psii_fluorescence_yield,
                psii_constitutive_heat_yield,
                psii_regulated_heat_yield,
                psi_photochemical_yield,
                psi_fluorescence_yield,
                psi_constitutive_heat_yield,
                psi_regulated_heat_yield,
                total_fluorescence_yield,
                psii_open_fraction,
                psi_open_fraction,
                psii_quenching_rate,
                npq,
            ),
        )
        return new{FT}(
            status,
            boundary_roundoff,
            psii_photon_flux,
            psi_photon_flux,
            total_photon_flux,
            psii_photochemical_yield,
            psii_fluorescence_yield,
            psii_constitutive_heat_yield,
            psii_regulated_heat_yield,
            psi_photochemical_yield,
            psi_fluorescence_yield,
            psi_constitutive_heat_yield,
            psi_regulated_heat_yield,
            total_fluorescence_yield,
            psii_open_fraction,
            psi_open_fraction,
            psii_quenching_rate,
            npq,
        )
    end
    function CytochromeFluorescenceResult(
        status::FluorescenceStatus,
        boundary_roundoff::Bool,
        psii_photon_flux::FT,
        psi_photon_flux::FT,
        total_photon_flux::FT,
        psii_photochemical_yield::FT,
        psii_fluorescence_yield::FT,
        psii_constitutive_heat_yield::FT,
        psii_regulated_heat_yield::FT,
        psi_photochemical_yield::FT,
        psi_fluorescence_yield::FT,
        psi_constitutive_heat_yield::FT,
        psi_regulated_heat_yield::FT,
        total_fluorescence_yield::FT,
        psii_open_fraction::FT,
        psi_open_fraction::FT,
        psii_quenching_rate::FT,
        npq::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeFluorescenceResult,
            (
                status,
                boundary_roundoff,
                psii_photon_flux,
                psi_photon_flux,
                total_photon_flux,
                psii_photochemical_yield,
                psii_fluorescence_yield,
                psii_constitutive_heat_yield,
                psii_regulated_heat_yield,
                psi_photochemical_yield,
                psi_fluorescence_yield,
                psi_constitutive_heat_yield,
                psi_regulated_heat_yield,
                total_fluorescence_yield,
                psii_open_fraction,
                psi_open_fraction,
                psii_quenching_rate,
                npq,
            ),
        )
        return new{FT}(
            status,
            boundary_roundoff,
            psii_photon_flux,
            psi_photon_flux,
            total_photon_flux,
            psii_photochemical_yield,
            psii_fluorescence_yield,
            psii_constitutive_heat_yield,
            psii_regulated_heat_yield,
            psi_photochemical_yield,
            psi_fluorescence_yield,
            psi_constitutive_heat_yield,
            psi_regulated_heat_yield,
            total_fluorescence_yield,
            psii_open_fraction,
            psi_open_fraction,
            psii_quenching_rate,
            npq,
        )
    end
end

@inline fluorescence_valid(r::CytochromeFluorescenceResult) =
    r.status == FluorescenceLit || r.status == FluorescenceDark

@inline function unavailable_cytochrome_fluorescence(::Type{FT}, status) where {FT}
    nan = FT(NaN)
    flux = status == FluorescenceDark ? zero(FT) : nan
    return CytochromeFluorescenceResult(
        status, false,
        flux, flux, flux,
        nan, nan, nan, nan,
        nan, nan, nan, nan,
        nan, nan, nan, nan, nan,
    )
end

@inline function fluorescence(
    model::CytochromeFluorescence{FTM},
    x::CoupledPhotosystemInputs{FTX},
) where {FTM<:Real,FTX<:Real}
    compatible_scalar(FTM, FTX) || throw(MethodError(fluorescence, (model, x)))
    S = promote_type(FTM, FTX)
    q2_exc, q1_exc = S(x.psii_excitation), S(x.psi_excitation)
    j2, j1, q2_open = S(x.linear_electron_flow), S(x.psi_electron_flow), S(x.psii_open_fraction)
    q2p, q1p, j2p, j1p, q2op = map(primal_scalar, (q2_exc, q1_exc, j2, j1, q2_open))
    all(isfinite, (q2p, q1p, j2p, j1p, q2op)) &&
        q2p >= zero(q2p) && q1p >= zero(q1p) &&
        j2p >= zero(j2p) && j1p >= zero(j1p) &&
        zero(q2op) <= q2op <= one(q2op) ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceInvalidInput)

    if iszero(q2p) && iszero(q1p)
        return unavailable_cytochrome_fluorescence(
            S,
            iszero(j2p) && iszero(j1p) ? FluorescenceDark : FluorescenceInvalidInput,
        )
    elseif iszero(q2p)
        return unavailable_cytochrome_fluorescence(
            S,
            iszero(j2p) ? FluorescenceIncompatibleRates : FluorescenceInvalidInput,
        )
    elseif iszero(q1p)
        return unavailable_cytochrome_fluorescence(
            S,
            iszero(j1p) ? FluorescenceIncompatibleRates : FluorescenceInvalidInput,
        )
    end

    # 1. PSI excitation partitioning (oxidized P700⁺ quenches closed traps at K_X1).
    phi_p1 = j1 / q1_exc
    phi_p1_max = S(psi_maximum_photochemical_yield(model.psi))
    isfinite(phi_p1) && phi_p1 < one(S) && phi_p1 - phi_p1_max <= model.yield_tolerance ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceIncompatibleRates)
    boundary1 = phi_p1 > phi_p1_max
    kf1 = S(model.psi.fluorescence_rate)
    kd1 = S(model.psi.constitutive_heat_rate)
    kp1 = S(model.psi.maximum_photochemical_rate)
    kx1 = S(model.psi.closed_quenching_rate)
    s_open1 = kf1 + kd1 + kp1
    s_closed1 = kf1 + kd1 + kx1
    if boundary1
        q1_open = one(S)
        rem1 = one(S) - phi_p1
        base1 = kf1 + kd1
        phi_f1 = rem1 * (kf1 / base1)
        phi_d1 = rem1 * (kd1 / base1)
        phi_n1 = zero(S)
    else
        q1_open = phi_p1 / phi_p1_max
        closed1 = one(S) - q1_open
        phi_f1 = q1_open * (kf1 / s_open1) + closed1 * (kf1 / s_closed1)
        phi_d1 = q1_open * (kd1 / s_open1) + closed1 * (kd1 / s_closed1)
        phi_n1 = closed1 * (kx1 / s_closed1)
    end

    # 2. PSII connected-units feedback quenching K_N >= 0 and excitation partitioning.
    phi_p2 = j2 / q2_exc
    q2_open > zero(S) && isfinite(phi_p2) && zero(S) < phi_p2 < one(S) ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceIncompatibleRates)

    kf2 = S(model.psii.fluorescence_rate)
    kd2 = S(model.psii.constitutive_heat_rate)
    kp2 = S(model.psii.maximum_photochemical_rate)
    ku2 = S(model.psii.exciton_transfer_rate)
    kns = S(model.psii.sustained_quenching_rate)

    d0 = kf2 + kd2 + kns
    s_base = d0 + ku2
    d0_hat = d0 / s_base
    ku_hat = ku2 / s_base
    kp_hat = kp2 / s_base
    m0_hat = d0_hat + kp_hat * (d0_hat + q2_open * ku_hat)
    phi_p2_max = (q2_open * kp_hat) / m0_hat

    phi_p2 - phi_p2_max <= model.yield_tolerance ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceIncompatibleRates)
    boundary2 = phi_p2 > phi_p2_max

    kn = if boundary2
        zero(S)
    else
        # Solve phi_p2 * kn_hat^2 + b_hat * kn_hat - p0_hat = 0 directly for kn_hat >= 0,
        # avoiding S_1 - S_base cancellation near kn = 0, rate^2 overflow in Float32,
        # and sqrt(p0_hat) derivative blow-up at p0_hat = 0 (since b_hat >= phi_p2 > 0).
        raw_p0 = q2_open * kp_hat - phi_p2 * m0_hat
        p0_hat = raw_p0 > zero(raw_p0) ? raw_p0 : zero(S)
        kn_hat = if exact_zero_scalar(p0_hat)
            zero(S)
        else
            b_hat = phi_p2 * (one(S) + d0_hat + kp_hat) - q2_open * kp_hat
            two = one(S) + one(S)
            four = two + two
            disc = b_hat * b_hat + four * (phi_p2 * p0_hat)
            radical = disc > zero(disc) ? sqrt(disc) : zero(S)
            b_hat >= zero(S) ? (two * p0_hat) / (b_hat + radical) :
                               (radical - b_hat) / (two * phi_p2)
        end
        s_base * kn_hat
    end

    kdiss2 = d0 + kn
    isfinite(kn) && kn >= zero(S) && isfinite(kdiss2) && kdiss2 > zero(S) ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceIncompatibleRates)

    rem2 = one(S) - phi_p2
    phi_f2 = rem2 * (kf2 / kdiss2)
    phi_d2 = rem2 * (kd2 / kdiss2)
    phi_n2 = rem2 * ((kns + kn) / kdiss2)
    npq = (kn + kns) / (kf2 + kd2)

    psii_flux = q2_exc * phi_f2
    psi_flux = q1_exc * phi_f1
    total_flux = psii_flux + psi_flux
    q_total = q2_exc + q1_exc
    total_yield = total_flux / q_total

    all(isfinite, (psii_flux, psi_flux, total_flux, total_yield, npq)) ||
        return unavailable_cytochrome_fluorescence(S, FluorescenceIncompatibleRates)

    return CytochromeFluorescenceResult{S}(
        FluorescenceLit,
        boundary1 || boundary2,
        psii_flux,
        psi_flux,
        total_flux,
        phi_p2,
        phi_f2,
        phi_d2,
        phi_n2,
        phi_p1,
        phi_f1,
        phi_d1,
        phi_n1,
        total_yield,
        q2_open,
        q1_open,
        kn,
        npq,
    )
end

"""
    photosystem_flux_inputs(photosynthesis::CytochromePhotosynthesis,
        site::LeafEnvironment[, excitation])

Extract coupled PSII and PSI excitations, realized linear (`J_P680`) and PSI
(`J_P700`) electron flows, and open PSII fraction (`q₂`) from signed Cytochrome
turnover at biochemical-site `site.co2` (`Cc`).
"""
@inline function photosystem_flux_inputs(
    p::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation},
    site::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L,FT<:Real}
    turnover = biochemical_turnover(p, site, excitation)
    q2 = psii_excitation(p.electron_transport, site.absorbed_par, excitation)
    q1 = psi_excitation(p.electron_transport, site.absorbed_par, excitation)
    OutFT = promote_type(typeof(q2), typeof(q1), typeof(turnover.linear_electron_flow))
    nan = OutFT(NaN)
    t_p, par_p, co2_p, o2_p = map(
        primal_scalar,
        (site.temperature, site.absorbed_par, site.co2, site.oxygen),
    )
    isfinite(t_p) && t_p > zero(t_p) &&
        isfinite(par_p) && par_p >= zero(par_p) &&
        isfinite(co2_p) && co2_p >= zero(co2_p) &&
        isfinite(o2_p) && o2_p >= zero(o2_p) &&
        valid_excitation(excitation, site.absorbed_par) ||
        return CoupledPhotosystemInputs(nan, nan, nan, nan, nan)
    return CoupledPhotosystemInputs(
        OutFT(q2),
        OutFT(q1),
        OutFT(turnover.linear_electron_flow),
        OutFT(turnover.psi_electron_flow),
        OutFT(turnover.psii_open_fraction),
    )
end

"""
    minimum_psii_partition(fluorescence_model::CytochromeFluorescence,
        stoichiometry::CytochromeStoichiometry;
        linear_atp_per_electron, cyclic_atp_per_electron)

Compute the minimum broadband PSII photon fraction `β_min` ensuring
`Φ_P2 ≤ Φ_P2,max(q₂)` (`K_N ≥ 0`) across all positive irradiances, Cytochrome
capacities, and CO₂ partial pressures under [`SerialLimitation`](@ref)
(Johnson & Berry 2021, Eq. 20).
"""
function minimum_psii_partition(
    fl::CytochromeFluorescence{FTF},
    s::CytochromeStoichiometry{FTS};
    linear_atp_per_electron::Real,
    cyclic_atp_per_electron::Real,
) where {FTF<:Real,FTS<:Real}
    _, _, nl, nc = promote_compatible(
        minimum_psii_partition, zero(FTF), zero(FTS), linear_atp_per_electron, cyclic_atp_per_electron,
    )
    FT = typeof(nl)
    nlp, ncp = primal_scalar(nl), primal_scalar(nc)
    isfinite(nlp) && nlp >= zero(nlp) &&
        isfinite(ncp) && ncp > zero(ncp) ||
        throw(ArgumentError("Invalid ATP/electron coupling ratios"))
    all(x -> (xp = primal_scalar(x); isfinite(xp) && xp > zero(xp)), (s.electrons_per_carboxylation,
        s.electrons_per_photorespiratory_co2, s.atp_per_carboxylation,
        s.atp_per_photorespiratory_co2)) ||
        throw(ArgumentError("Invalid Cytochrome stoichiometry"))
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(nl, nc, s)
    c1p, c2p = primal_scalar(cef_c1), primal_scalar(cef_c2)
    isfinite(c1p) && isfinite(c2p) && c1p >= zero(c1p) && c2p >= zero(c2p) ||
        throw(ArgumentError("Implied minimum PSI/PSII electron ratio η_min < 1 implies negative cyclic flow"))
    eta_min = one(FT) + min(
        cef_c1 / s.electrons_per_carboxylation,
        cef_c2 / s.electrons_per_photorespiratory_co2,
    )
    phi1_max = FT(psi_maximum_photochemical_yield(fl))
    phi2_max = FT(psii_maximum_photochemical_yield(fl, one(FTF)))
    return FT(phi1_max / (phi1_max + eta_min * phi2_max))
end

function validate_fluorescence_coupling(
    fl::CytochromeFluorescence{FTF},
    et::CytochromeElectronParameters{FTE},
    s::CytochromeStoichiometry{FTS};
    linear_atp_per_electron::Real,
    cyclic_atp_per_electron::Real,
) where {FTF<:Real,FTE<:Real,FTS<:Real}
    _, _, _, nl, nc = promote_compatible(
        validate_fluorescence_coupling,
        zero(FTF), zero(FTE), zero(FTS), linear_atp_per_electron, cyclic_atp_per_electron,
    )
    FT = typeof(nl)
    et.limitation isa SerialLimitation ||
        throw(ArgumentError("Coupled Cytochrome fluorescence requires SerialLimitation() in CytochromeElectronParameters"))
    expected_psi = FT(psi_maximum_photochemical_yield(fl))
    primal_scalar(abs(et.psi_maximum_yield - expected_psi)) <= primal_scalar(fl.yield_tolerance) ||
        throw(ArgumentError("Cytochrome electron transport and PSI fluorescence maximum yields disagree"))
    beta_min = minimum_psii_partition(fl, s;
        linear_atp_per_electron = nl, cyclic_atp_per_electron = nc)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(nl, nc, s)
    eta_min = one(FT) + min(
        cef_c1 / s.electrons_per_carboxylation,
        cef_c2 / s.electrons_per_photorespiratory_co2,
    )
    phi2_max = FT(psii_maximum_photochemical_yield(fl, one(FTF)))
    primal_scalar(et.psii_fraction) > zero(primal_scalar(et.psii_fraction)) ||
        throw(ArgumentError("Broadband PSII fraction must exceed β_min = $(beta_min)"))
    low_light_phi2 = ((one(FT) - et.psii_fraction) * et.psi_maximum_yield) /
                     (et.psii_fraction * eta_min)
    primal_scalar(low_light_phi2 - phi2_max) <= primal_scalar(fl.yield_tolerance) ||
        throw(ArgumentError("Broadband PSII fraction $(et.psii_fraction) < β_min ($(beta_min)) implies super-maximal PSII yield"))
    return nothing
end

function validate_fluorescence_coupling(
    fl::CytochromeFluorescence{FTF},
    photo::CytochromePhotosynthesis,
    temperature::Real,
    oxygen::Real,
) where {FTF<:Real}
    _, t_s, o_s = promote_compatible(validate_fluorescence_coupling, zero(FTF), temperature, oxygen)
    r = biochemical_rates(photo.temperature, t_s, o_s)
    return validate_fluorescence_coupling(
        fl,
        photo.electron_transport,
        photo.stoichiometry;
        linear_atp_per_electron = r.linear_atp_per_electron,
        cyclic_atp_per_electron = r.cyclic_atp_per_electron,
    )
end
