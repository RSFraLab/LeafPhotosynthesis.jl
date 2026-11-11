"""
    PSIIRateParameters(; fluorescence_rate, constitutive_heat_rate,
        maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate)

First-order PSII excitation decay and inter-unit transfer rate constants sharing
a common rate normalization (conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F2 > 0`)
- `constitutive_heat_rate` (`K_D2 ≥ 0`)
- `maximum_photochemical_rate` (`K_P2 > 0`)
- `exciton_transfer_rate` (`K_U2 ≥ 0`; `0` gives the separate-units puddle limit,
  `> 0` gives the connected-units exciton-sharing model)
- `sustained_quenching_rate` (`K_N,sus ≥ 0`; dark-sustained thermal dissipation)
"""
struct PSIIRateParameters{FT<:AbstractFloat}
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
    ) where {FT<:AbstractFloat}
        kf, kd, kp, ku, kns = fluorescence_rate, constitutive_heat_rate,
            maximum_photochemical_rate, exciton_transfer_rate, sustained_quenching_rate
        d0 = kf + kd + kns
        s_base = d0 + ku
        total_open = d0 + kp
        phi_open = kp / total_open
        d0_hat = d0 / s_base
        ku_hat = ku / s_base
        kp_hat = kp / s_base
        phi_conn = kp_hat / (d0_hat + kp_hat * (d0_hat + ku_hat))
        all(isfinite, (kf, kd, kp, ku, kns, d0, s_base, total_open, s_base + kp, phi_conn)) &&
            kf > zero(FT) && kd >= zero(FT) && kp > zero(FT) &&
            ku >= zero(FT) && kns >= zero(FT) &&
            zero(FT) < phi_open < one(FT) && zero(FT) < phi_conn < one(FT) ||
            throw(ArgumentError("Invalid or unrepresentable PSII rate parameters"))
        new{FT}(kf, kd, kp, ku, kns)
    end
end
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

First-order PSI excitation decay rate constants sharing one normalization
(conventionally ns⁻¹, Johnson & Berry 2021, Table 1):
- `fluorescence_rate` (`K_F1 > 0`)
- `constitutive_heat_rate` (`K_D1 ≥ 0`)
- `maximum_photochemical_rate` (`K_P1 > 0` in open P700 centers)
- `closed_quenching_rate` (`K_X1 > 0` radiationless decay in oxidized P700⁺ / closed PSI centers)
"""
struct PSIRateParameters{FT<:AbstractFloat}
    fluorescence_rate::FT
    constitutive_heat_rate::FT
    maximum_photochemical_rate::FT
    closed_quenching_rate::FT
    function PSIRateParameters(
        fluorescence_rate::FT,
        constitutive_heat_rate::FT,
        maximum_photochemical_rate::FT,
        closed_quenching_rate::FT,
    ) where {FT<:AbstractFloat}
        kf, kd, kp, kx = fluorescence_rate, constitutive_heat_rate,
            maximum_photochemical_rate, closed_quenching_rate
        s_open = kf + kd + kp
        s_closed = kf + kd + kx
        phi_open = kp / s_open
        phi_closed = kx / s_closed
        all(isfinite, (kf, kd, kp, kx, s_open, s_closed)) &&
            kf > zero(FT) && kd >= zero(FT) && kp > zero(FT) && kx > zero(FT) &&
            zero(FT) < phi_open < one(FT) && zero(FT) < phi_closed < one(FT) ||
            throw(ArgumentError("Invalid or unrepresentable PSI rate parameters"))
        new{FT}(kf, kd, kp, kx)
    end
end
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

"""
    CytochromeFluorescence(; psii, psi, yield_tolerance)

Mechanistic steady-state PSII + PSI fluorescence model coupled through
Cytochrome b₆f redox poise and energy-dependent feedback quenching `K_N`
(Johnson & Berry 2021, Eqs. 7–30). `yield_tolerance` is a roundoff-sized
allowance in dimensionless yield units, matching `RateCompetitionFluorescence`.
"""
struct CytochromeFluorescence{FT<:AbstractFloat} <: AbstractFluorescence
    psii::PSIIRateParameters{FT}
    psi::PSIRateParameters{FT}
    yield_tolerance::FT
    function CytochromeFluorescence(
        psii::PSIIRateParameters{FT},
        psi::PSIRateParameters{FT},
        yield_tolerance::FT,
    ) where {FT<:AbstractFloat}
        phi2_max = psii_maximum_photochemical_yield(psii, one(FT))
        phi1_max = psi_maximum_photochemical_yield(psi)
        headroom = min(one(FT) - phi2_max, one(FT) - phi1_max)
        isfinite(yield_tolerance) && zero(FT) <= yield_tolerance < headroom ||
            throw(ArgumentError("Invalid Cytochrome fluorescence yield tolerance"))
        new{FT}(psii, psi, yield_tolerance)
    end
end
CytochromeFluorescence(; psii, psi, yield_tolerance) =
    CytochromeFluorescence(psii, psi, yield_tolerance)

"""Open-center PSI maximum photochemical quantum yield `K_P1 / (K_F1 + K_D1 + K_P1)`."""
@inline psi_maximum_photochemical_yield(p::PSIRateParameters) =
    p.maximum_photochemical_rate /
    (p.fluorescence_rate + p.constitutive_heat_rate + p.maximum_photochemical_rate)
@inline psi_maximum_photochemical_yield(m::CytochromeFluorescence) =
    psi_maximum_photochemical_yield(m.psi)

"""
    psii_maximum_photochemical_yield(psii_or_model[, psii_open_fraction = 1])

PSII photochemical yield at zero dynamic feedback quenching (`K_N = 0`) and open
reaction-center fraction `psii_open_fraction ∈ [0, 1]`, accounting for connected
exciton sharing (`K_U2`) and any dark-sustained quenching (`K_N,sus`).
"""
@inline function psii_maximum_photochemical_yield(
    p::PSIIRateParameters{FT},
    psii_open_fraction::FT = one(FT),
) where {FT}
    d0 = p.fluorescence_rate + p.constitutive_heat_rate + p.sustained_quenching_rate
    s_base = d0 + p.exciton_transfer_rate
    d0_hat = d0 / s_base
    ku_hat = p.exciton_transfer_rate / s_base
    kp_hat = p.maximum_photochemical_rate / s_base
    m0_hat = d0_hat + kp_hat * (d0_hat + psii_open_fraction * ku_hat)
    return (psii_open_fraction * kp_hat) / m0_hat
end
@inline psii_maximum_photochemical_yield(
    m::CytochromeFluorescence{FT},
    psii_open_fraction::FT = one(FT),
) where {FT} = psii_maximum_photochemical_yield(m.psii, psii_open_fraction)

"""
    CoupledPhotosystemInputs{FT}(psii_excitation, psi_excitation,
        linear_electron_flow, psi_electron_flow, psii_open_fraction)

Scalar inputs to `CytochromeFluorescence`:
- `psii_excitation`, `psi_excitation`: absorbed photon fluxes delivered to PSII
  and PSI [μmol photons m⁻² one-sided leaf s⁻¹]
- `linear_electron_flow`: realized PSII linear electron flux `J_P680` [μmol e⁻ m⁻² leaf s⁻¹]
- `psi_electron_flow`: realized PSI electron flux `J_P700 = J_P680 + J_CEF` [μmol e⁻ m⁻² leaf s⁻¹]
- `psii_open_fraction`: open PSII reaction-center fraction `q₂ ∈ [0, 1]`
"""
Base.@kwdef struct CoupledPhotosystemInputs{FT<:AbstractFloat}
    psii_excitation::FT
    psi_excitation::FT
    linear_electron_flow::FT
    psi_electron_flow::FT
    psii_open_fraction::FT
end

"""
Concrete scalar result of `CytochromeFluorescence`.
- `psii_photon_flux`, `psi_photon_flux`, `total_photon_flux` are internal leaf
  fluorescence emission rates [μmol photons m⁻² one-sided leaf s⁻¹] **before**
  within-leaf optical escape/reabsorption and spectral shaping.
- Per-photosystem yields (`psii_*_yield` and `psi_*_yield`) each partition their
  photosystem's excitation events (`∑ Φ = 1`), not absorbed joules.
- `npq = (K_N + K_N,sus) / (K_F2 + K_D2)` is total internal PSII Stern–Volmer
  non-photochemical quenching `(Fm - Fm') / Fm'` referenced to the **unquenched
  pre-stress dark maximum `Fm = K_F2 / (K_F2 + K_D2)`** (whereas dark-adapted PAM
  measurements on a photoinhibited leaf with unrelaxed `K_N,sus > 0` measure the
  dynamic component `K_N / (K_F2 + K_D2 + K_N,sus)`).
- `FluorescenceDark` (`Q₂ = Q₁ = J_P680 = J_P700 = 0`) reports zero photon fluxes
  and `NaN` flux-derived yields/rates. Lit steady-state `K_N` inversion requires
  `q₂ > 0` and `0 < Φ_P2 ≤ Φ_P2,max(q₂) < 1`; zero-turnover illuminated states
  (`Vcmax = 0` where `K_N → ∞`, or `Vqmax = 0` where `q₂ = Φ_P2 = 0` is indeterminate)
  and PSII-starved far-red `PrescribedDualExcitation` (`Φ_P2 > Φ_P2,max(q₂)`) return
  `FluorescenceIncompatibleRates` with `NaN` fluxes without throwing on GPU devices.
"""
struct CytochromeFluorescenceResult{FT}
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

"""
    fluorescence(model::CytochromeFluorescence, inputs::CoupledPhotosystemInputs)

Evaluate coupled PSII + PSI fluorescence. Invalid inputs return `FluorescenceInvalidInput`;
domains outside steady-state `K_N ≥ 0` inversion (including lit `q₂ = 0`, lit `Φ_P2 = 0`,
or far-red shade with `Φ_P2 > Φ_P2,max(q₂) + yield_tolerance`) return
`FluorescenceIncompatibleRates` without throwing.
"""
@inline function fluorescence(
    model::CytochromeFluorescence{FT},
    x::CoupledPhotosystemInputs{FT},
) where {FT}
    q2_exc, q1_exc = x.psii_excitation, x.psi_excitation
    j2, j1, q2_open = x.linear_electron_flow, x.psi_electron_flow, x.psii_open_fraction
    all(isfinite, (q2_exc, q1_exc, j2, j1, q2_open)) &&
        q2_exc >= zero(FT) && q1_exc >= zero(FT) &&
        j2 >= zero(FT) && j1 >= zero(FT) &&
        zero(FT) <= q2_open <= one(FT) ||
        return unavailable_cytochrome_fluorescence(FT, FluorescenceInvalidInput)

    if iszero(q2_exc) && iszero(q1_exc)
        return unavailable_cytochrome_fluorescence(
            FT,
            iszero(j2) && iszero(j1) ? FluorescenceDark : FluorescenceInvalidInput,
        )
    elseif iszero(q2_exc)
        return unavailable_cytochrome_fluorescence(
            FT,
            iszero(j2) ? FluorescenceIncompatibleRates : FluorescenceInvalidInput,
        )
    elseif iszero(q1_exc)
        return unavailable_cytochrome_fluorescence(
            FT,
            iszero(j1) ? FluorescenceIncompatibleRates : FluorescenceInvalidInput,
        )
    end

    # 1. PSI excitation partitioning (oxidized P700⁺ quenches closed traps at K_X1).
    phi_p1 = j1 / q1_exc
    phi_p1_max = psi_maximum_photochemical_yield(model.psi)
    isfinite(phi_p1) && phi_p1 < one(FT) && phi_p1 - phi_p1_max <= model.yield_tolerance ||
        return unavailable_cytochrome_fluorescence(FT, FluorescenceIncompatibleRates)
    boundary1 = phi_p1 > phi_p1_max
    kf1 = model.psi.fluorescence_rate
    kd1 = model.psi.constitutive_heat_rate
    kp1 = model.psi.maximum_photochemical_rate
    kx1 = model.psi.closed_quenching_rate
    s_open1 = kf1 + kd1 + kp1
    s_closed1 = kf1 + kd1 + kx1
    if boundary1
        q1_open = one(FT)
        rem1 = one(FT) - phi_p1
        base1 = kf1 + kd1
        phi_f1 = rem1 * (kf1 / base1)
        phi_d1 = rem1 * (kd1 / base1)
        phi_n1 = zero(FT)
    else
        q1_open = phi_p1 / phi_p1_max
        closed1 = one(FT) - q1_open
        phi_f1 = q1_open * (kf1 / s_open1) + closed1 * (kf1 / s_closed1)
        phi_d1 = q1_open * (kd1 / s_open1) + closed1 * (kd1 / s_closed1)
        phi_n1 = closed1 * (kx1 / s_closed1)
    end

    # 2. PSII connected-units feedback quenching K_N >= 0 and excitation partitioning.
    phi_p2 = j2 / q2_exc
    q2_open > zero(FT) && isfinite(phi_p2) && zero(FT) < phi_p2 < one(FT) ||
        return unavailable_cytochrome_fluorescence(FT, FluorescenceIncompatibleRates)

    kf2 = model.psii.fluorescence_rate
    kd2 = model.psii.constitutive_heat_rate
    kp2 = model.psii.maximum_photochemical_rate
    ku2 = model.psii.exciton_transfer_rate
    kns = model.psii.sustained_quenching_rate

    d0 = kf2 + kd2 + kns
    s_base = d0 + ku2
    d0_hat = d0 / s_base
    ku_hat = ku2 / s_base
    kp_hat = kp2 / s_base
    m0_hat = d0_hat + kp_hat * (d0_hat + q2_open * ku_hat)
    phi_p2_max = (q2_open * kp_hat) / m0_hat

    phi_p2 - phi_p2_max <= model.yield_tolerance ||
        return unavailable_cytochrome_fluorescence(FT, FluorescenceIncompatibleRates)
    boundary2 = phi_p2 > phi_p2_max

    kn = if boundary2
        zero(FT)
    else
        # Solve phi_p2 * kn_hat^2 + b_hat * kn_hat - p0_hat = 0 directly for kn_hat >= 0,
        # avoiding both S_1 - S_base cancellation near kn = 0 and rate^2 overflow in Float32.
        p0_hat = max(zero(FT), q2_open * kp_hat - phi_p2 * m0_hat)
        b_hat = phi_p2 * (one(FT) + d0_hat + kp_hat) - q2_open * kp_hat
        two = one(FT) + one(FT)
        radical = hypot(b_hat, two * sqrt(phi_p2) * sqrt(p0_hat))
        kn_hat = if iszero(p0_hat)
            zero(FT)
        elseif b_hat >= zero(FT)
            (two * p0_hat) / (b_hat + radical)
        else
            (radical - b_hat) / (two * phi_p2)
        end
        s_base * kn_hat
    end

    kdiss2 = d0 + kn
    isfinite(kn) && kn >= zero(FT) && isfinite(kdiss2) && kdiss2 > zero(FT) ||
        return unavailable_cytochrome_fluorescence(FT, FluorescenceIncompatibleRates)

    # Because K_F2, K_D2, K_N,sus and K_N act with equal rates in open and closed
    # PSII units, non-photochemical excitons branch strictly in proportion to their
    # rates even when K_U2 > 0 transfers excitons across units.
    rem2 = one(FT) - phi_p2
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
        return unavailable_cytochrome_fluorescence(FT, FluorescenceIncompatibleRates)

    return CytochromeFluorescenceResult(
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

Extract coupled PSII/PSI excitations, realized linear (`J_P680`) and PSI (`J_P700`)
electron flows, and open PSII fraction (`q₂`) from signed Cytochrome turnover at
biochemical-site `site.co2` (`Cc`). Refuses ambient `GasExchangeEnvironment` so
callers explicitly verify gas-exchange convergence first.
"""
@inline function photosystem_flux_inputs(
    p::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation},
    site::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L,FT}
    nan = FT(NaN)
    isfinite(site.temperature) && site.temperature > zero(FT) &&
        isfinite(site.absorbed_par) && site.absorbed_par >= zero(FT) &&
        isfinite(site.co2) && site.co2 >= zero(FT) &&
        isfinite(site.oxygen) && site.oxygen >= zero(FT) &&
        valid_excitation(excitation, site.absorbed_par) ||
        return CoupledPhotosystemInputs(nan, nan, nan, nan, nan)
    turnover = biochemical_turnover(p, site, excitation)
    return CoupledPhotosystemInputs(
        psii_excitation(p.electron_transport, site.absorbed_par, excitation),
        psi_excitation(p.electron_transport, site.absorbed_par, excitation),
        turnover.linear_electron_flow,
        turnover.psi_electron_flow,
        turnover.psii_open_fraction,
    )
end

"""
    minimum_psii_partition(fluorescence_model::CytochromeFluorescence,
        stoichiometry::CytochromeStoichiometry;
        linear_atp_per_electron, cyclic_atp_per_electron)

Minimum broadband PSII photon fraction `β_min` ensuring `Φ_P2 ≤ Φ_P2,max(q₂)`
(`K_N ≥ 0`) across all light intensities `Q_abs > 0`, all `Vqmax > 0`, and all
`Cc ∈ [0, ∞)` under `SerialLimitation()` (Johnson & Berry 2021, Eq. 20):
`β_min = Φ_P1,max / (Φ_P1,max + η_min * Φ_P2,max(1))`.
With Johnson & Berry (2021) Table 1 rates (`Φ_P1,max = 14.5/15.1`, `Φ_P2,max(1) = 4.5/5.1`,
`η_min = 1`), `β_min ≈ 0.52114` (explaining why naive `β = 0.50` produces negative `K_N`).
"""
function minimum_psii_partition(
    fl::CytochromeFluorescence{FT},
    s::CytochromeStoichiometry{FT};
    linear_atp_per_electron::FT,
    cyclic_atp_per_electron::FT,
) where {FT}
    isfinite(linear_atp_per_electron) && linear_atp_per_electron >= zero(FT) &&
        isfinite(cyclic_atp_per_electron) && cyclic_atp_per_electron > zero(FT) ||
        throw(ArgumentError("Invalid ATP/electron coupling ratios"))
    all(x -> isfinite(x) && x > zero(FT), (s.electrons_per_carboxylation,
        s.electrons_per_photorespiratory_co2, s.atp_per_carboxylation,
        s.atp_per_photorespiratory_co2)) ||
        throw(ArgumentError("Invalid Cytochrome stoichiometry"))
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        linear_atp_per_electron, cyclic_atp_per_electron, s)
    isfinite(cef_c1) && isfinite(cef_c2) && cef_c1 >= zero(FT) && cef_c2 >= zero(FT) ||
        throw(ArgumentError("Implied minimum PSI/PSII electron ratio η_min < 1 implies negative cyclic flow"))
    eta_min = one(FT) + min(cef_c1 / s.electrons_per_carboxylation,
                            cef_c2 / s.electrons_per_photorespiratory_co2)
    phi1_max = psi_maximum_photochemical_yield(fl)
    phi2_max = psii_maximum_photochemical_yield(fl, one(FT))
    return phi1_max / (phi1_max + eta_min * phi2_max)
end

"""
    validate_fluorescence_coupling(fl::CytochromeFluorescence,
        et::CytochromeElectronParameters, stoichiometry::CytochromeStoichiometry;
        linear_atp_per_electron, cyclic_atp_per_electron)
    validate_fluorescence_coupling(fl::CytochromeFluorescence,
        photosynthesis::CytochromePhotosynthesis, temperature, oxygen)

Host setup check verifying that:
1. `et.limitation isa SerialLimitation` (required so `q₂ = Vqmax / (Vqmax + Φ_P1,max Q₁)`
   guarantees `Φ_P2(Q) ≤ Φ_P2,max(q₂(Q))` across all `Q_abs > 0`);
2. `et.psi_maximum_yield` agrees with `psi_maximum_photochemical_yield(fl)` within `fl.yield_tolerance`;
3. ATP/NADPH stoichiometry implies nonnegative cyclic electron flow (`η_min ≥ 1`);
4. broadband `et.psii_fraction` satisfies `β ≥ β_min` so low-light PSII yield never
   exceeds `Φ_P2,max(1)` (preventing unphysical negative `K_N`).
"""
function validate_fluorescence_coupling(
    fl::CytochromeFluorescence{FT},
    et::CytochromeElectronParameters{FT},
    s::CytochromeStoichiometry{FT};
    linear_atp_per_electron::FT,
    cyclic_atp_per_electron::FT,
) where {FT}
    et.limitation isa SerialLimitation ||
        throw(ArgumentError("Coupled Cytochrome fluorescence requires SerialLimitation() in CytochromeElectronParameters"))
    expected_psi = psi_maximum_photochemical_yield(fl)
    abs(et.psi_maximum_yield - expected_psi) <= fl.yield_tolerance ||
        throw(ArgumentError("Cytochrome electron transport and PSI fluorescence maximum yields disagree"))
    beta_min = minimum_psii_partition(fl, s;
        linear_atp_per_electron, cyclic_atp_per_electron)
    # Check low-light PSII yield directly in yield units against fl.yield_tolerance.
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        linear_atp_per_electron, cyclic_atp_per_electron, s)
    eta_min = one(FT) + min(cef_c1 / s.electrons_per_carboxylation,
                            cef_c2 / s.electrons_per_photorespiratory_co2)
    phi2_max = psii_maximum_photochemical_yield(fl, one(FT))
    et.psii_fraction > zero(FT) ||
        throw(ArgumentError("Broadband PSII fraction must exceed β_min = $(beta_min)"))
    low_light_phi2 = ((one(FT) - et.psii_fraction) * et.psi_maximum_yield) /
                     (et.psii_fraction * eta_min)
    low_light_phi2 - phi2_max <= fl.yield_tolerance ||
        throw(ArgumentError("Broadband PSII fraction $(et.psii_fraction) < β_min ($(beta_min)) implies super-maximal PSII yield"))
    return nothing
end

function validate_fluorescence_coupling(
    fl::CytochromeFluorescence{FT},
    photo::CytochromePhotosynthesis,
    temperature::FT,
    oxygen::FT,
) where {FT}
    r = biochemical_rates(photo.temperature, temperature, oxygen)
    return validate_fluorescence_coupling(
        fl,
        photo.electron_transport,
        photo.stoichiometry;
        linear_atp_per_electron = r.linear_atp_per_electron,
        cyclic_atp_per_electron = r.cyclic_atp_per_electron,
    )
end
