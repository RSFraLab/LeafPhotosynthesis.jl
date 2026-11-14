"""
    CytochromeCapacity(; site_density, turnover_rate)
    CytochromeCapacity(site_density, turnover_rate)

Parameterize Cytochrome b₆f maximum plastoquinol-oxidation capacity from active
site density `site_density` [μmol sites m⁻² leaf] and catalytic turnover law
`turnover_rate` [mol e⁻ mol⁻¹ sites s⁻¹] (Johnson & Berry 2021):
``V_{q\\max}(T) = B \\cdot k_q(T)`` [μmol e⁻ m⁻² leaf s⁻¹].
"""
struct CytochromeCapacity{FT<:Real,K}
    site_density::FT
    turnover_rate::K
    function CytochromeCapacity(site_density::FT, turnover_rate::K) where {FT<:Real,K}
        require_continuous_scalar(FT, CytochromeCapacity, (site_density, turnover_rate))
        sp = primal_scalar(site_density)
        isfinite(sp) && sp >= zero(sp) ||
            throw(ArgumentError("Cytochrome b6f site_density must be finite and nonnegative"))
        return new{FT,K}(site_density, turnover_rate)
    end
end
CytochromeCapacity(; site_density, turnover_rate) =
    CytochromeCapacity(site_density, turnover_rate)

@inline temperature_response(c::CytochromeCapacity, temperature) =
    c.site_density * temperature_response(c.turnover_rate, temperature)

"""
    CytochromeTemperatureResponse(; vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron)

Bundle temperature and oxygen response laws for Cytochrome b₆f C3 biochemistry
(Johnson & Berry 2021).

`vqmax` may be a [`CytochromeCapacity`](@ref) or any temperature-response law
returning μmol e⁻ m⁻² leaf s⁻¹; `gamma_star` may be a temperature law or an
[`AbstractCompensationPoint`](@ref); `linear_atp_per_electron` (`η_L`) and
`cyclic_atp_per_electron` (`η_C`) use [`ConstantResponse`](@ref) or temperature laws.
"""
Base.@kwdef struct CytochromeTemperatureResponse{V,Q,R,C,O,G,NL,NC}
    vcmax::V
    vqmax::Q
    respiration::R
    kc::C
    ko::O
    gamma_star::G
    linear_atp_per_electron::NL
    cyclic_atp_per_electron::NC
end

"""
    CytochromeRates{FT}(; vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron)
    CytochromeRates(vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron)

Store evaluated Cytochrome b₆f kinetic capacities and ATP coupling efficiencies
at leaf temperature [K] and biochemical-site O₂ [Pa]:
- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹.
- `vqmax`: μmol electrons m⁻² leaf s⁻¹.
- `kc`, `ko`, `gamma_star`: Pa.
- `linear_atp_per_electron` (`η_L`), `cyclic_atp_per_electron` (`η_C`): mol ATP (mol e⁻)⁻¹.
"""
Base.@kwdef struct CytochromeRates{FT<:Real}
    vcmax::FT
    vqmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    linear_atp_per_electron::FT
    cyclic_atp_per_electron::FT
    function CytochromeRates{FT}(
        vcmax,
        vqmax,
        respiration,
        kc,
        ko,
        gamma_star,
        linear_atp_per_electron,
        cyclic_atp_per_electron,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeRates{FT},
            (
                vcmax,
                vqmax,
                respiration,
                kc,
                ko,
                gamma_star,
                linear_atp_per_electron,
                cyclic_atp_per_electron,
            ),
        )
        return new{FT}(
            vcmax,
            vqmax,
            respiration,
            kc,
            ko,
            gamma_star,
            linear_atp_per_electron,
            cyclic_atp_per_electron,
        )
    end
    function CytochromeRates(
        vcmax::FT,
        vqmax::FT,
        respiration::FT,
        kc::FT,
        ko::FT,
        gamma_star::FT,
        linear_atp_per_electron::FT,
        cyclic_atp_per_electron::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeRates,
            (
                vcmax,
                vqmax,
                respiration,
                kc,
                ko,
                gamma_star,
                linear_atp_per_electron,
                cyclic_atp_per_electron,
            ),
        )
        return new{FT}(
            vcmax,
            vqmax,
            respiration,
            kc,
            ko,
            gamma_star,
            linear_atp_per_electron,
            cyclic_atp_per_electron,
        )
    end
end
CytochromeRates(
    vc::Real, vq::Real, rd::Real, kc::Real, ko::Real, gs::Real, nl::Real, nc::Real,
) = CytochromeRates(
    promote_compatible(CytochromeRates, vc, vq, rd, kc, ko, gs, nl, nc)...,
)

@inline biochemical_rates(p::CytochromeTemperatureResponse, temperature) = CytochromeRates(
    temperature_response(p.vcmax, temperature),
    temperature_response(p.vqmax, temperature),
    temperature_response(p.respiration, temperature),
    temperature_response(p.kc, temperature),
    temperature_response(p.ko, temperature),
    temperature_response(p.gamma_star, temperature),
    temperature_response(p.linear_atp_per_electron, temperature),
    temperature_response(p.cyclic_atp_per_electron, temperature),
)

@inline biochemical_rates(p::CytochromeTemperatureResponse, temperature, oxygen) = CytochromeRates(
    temperature_response(p.vcmax, temperature),
    temperature_response(p.vqmax, temperature),
    temperature_response(p.respiration, temperature),
    temperature_response(p.kc, temperature),
    temperature_response(p.ko, temperature),
    compensation_point(p.gamma_star, temperature, oxygen),
    temperature_response(p.linear_atp_per_electron, temperature),
    temperature_response(p.cyclic_atp_per_electron, temperature),
)

"""
    CytochromeStoichiometry{FT}(; electrons_per_carboxylation,
        electrons_per_photorespiratory_co2, atp_per_carboxylation,
        atp_per_photorespiratory_co2)
    CytochromeStoichiometry(electrons_per_carboxylation,
        electrons_per_photorespiratory_co2, atp_per_carboxylation,
        atp_per_photorespiratory_co2)

Specify C3 Calvin–Benson and photorespiratory linear-electron and ATP demands
per CO₂ event (classically `4, 8, 3, 7` in Johnson & Berry 2021, Eq. 24).
"""
Base.@kwdef struct CytochromeStoichiometry{FT<:Real}
    electrons_per_carboxylation::FT
    electrons_per_photorespiratory_co2::FT
    atp_per_carboxylation::FT
    atp_per_photorespiratory_co2::FT
    function CytochromeStoichiometry{FT}(
        electrons_per_carboxylation,
        electrons_per_photorespiratory_co2,
        atp_per_carboxylation,
        atp_per_photorespiratory_co2,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeStoichiometry{FT},
            (
                electrons_per_carboxylation,
                electrons_per_photorespiratory_co2,
                atp_per_carboxylation,
                atp_per_photorespiratory_co2,
            ),
        )
        return new{FT}(
            electrons_per_carboxylation,
            electrons_per_photorespiratory_co2,
            atp_per_carboxylation,
            atp_per_photorespiratory_co2,
        )
    end
    function CytochromeStoichiometry(
        electrons_per_carboxylation::FT,
        electrons_per_photorespiratory_co2::FT,
        atp_per_carboxylation::FT,
        atp_per_photorespiratory_co2::FT,
    ) where {FT<:Real}
        require_continuous_scalar(
            FT,
            CytochromeStoichiometry,
            (
                electrons_per_carboxylation,
                electrons_per_photorespiratory_co2,
                atp_per_carboxylation,
                atp_per_photorespiratory_co2,
            ),
        )
        return new{FT}(
            electrons_per_carboxylation,
            electrons_per_photorespiratory_co2,
            atp_per_carboxylation,
            atp_per_photorespiratory_co2,
        )
    end
end
CytochromeStoichiometry(ec::Real, eg::Real, ac::Real, ag::Real) =
    CytochromeStoichiometry(promote_compatible(CytochromeStoichiometry, ec, eg, ac, ag)...)

"""
    cytochrome_cyclic_coefficients(linear_atp_per_electron, cyclic_atp_per_electron,
        stoichiometry::CytochromeStoichiometry) -> (cef_c1, cef_c2)

Compute cyclic PSI electron requirements per carboxylation and per
photorespiratory CO₂ release directly from the ATP deficit prior to adding
linear electron demand so `a_c == η_L * e_c` evaluates to exact zero in `Float32`.
"""
@inline function cytochrome_cyclic_coefficients(
    linear_atp_per_electron,
    cyclic_atp_per_electron,
    s::CytochromeStoichiometry,
)
    cef_c1 = (s.atp_per_carboxylation -
              linear_atp_per_electron * s.electrons_per_carboxylation) /
             cyclic_atp_per_electron
    cef_c2 = (s.atp_per_photorespiratory_co2 -
              linear_atp_per_electron * s.electrons_per_photorespiratory_co2) /
             cyclic_atp_per_electron
    return cef_c1, cef_c2
end

"""
    cytochrome_carbon_coefficients(rates::CytochromeRates,
        stoichiometry::CytochromeStoichiometry) -> (c1, c2)

Compute effective PSI electron costs `(c1, c2)` multiplying `Cc` and `Γ*` in
``J_{P700} = (c_1 C_c + c_2 Γ^*) V_c / C_c``.
"""
@inline function cytochrome_carbon_coefficients(r::CytochromeRates, s::CytochromeStoichiometry)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s,
    )
    return s.electrons_per_carboxylation + cef_c1,
           s.electrons_per_photorespiratory_co2 + cef_c2
end

"""
    psi_to_psii_electron_ratio(rates::CytochromeRates, co2,
        stoichiometry::CytochromeStoichiometry)

Compute the steady-state PSI-to-PSII electron flux ratio
``η = J_{P700} / J_{P680} = 1 + J_{\\mathrm{CEF}} / J_{P680}`` required to
balance Calvin–Benson and photorespiratory ATP/NADPH stoichiometry at
biochemical-site `co2` [Pa].
"""
@inline function psi_to_psii_electron_ratio(r::CytochromeRates, co2, s::CytochromeStoichiometry)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    lef_denom = s.electrons_per_carboxylation * co2 +
                s.electrons_per_photorespiratory_co2 * r.gamma_star
    S = typeof(c1 + c2 + lef_denom)
    iszero(primal_scalar(lef_denom)) && return S(c1 / s.electrons_per_carboxylation)
    return S((c1 * co2 + c2 * r.gamma_star) / lef_denom)
end

"""
    CytochromeElectronParameters(; psii_fraction, psi_maximum_yield, limitation)
    CytochromeElectronParameters(psii_fraction, psi_maximum_yield, limitation)

Specify broadband PSII photon partition `psii_fraction` (`β ∈ (0, 1)`),
open-center PSI maximum photochemical yield `psi_maximum_yield` (`Φ_P1,max ∈ (0, 1)`),
and rate co-limitation model `limitation` ([`SerialLimitation`](@ref) for
Johnson & Berry 2021, Eq. 20).
"""
struct CytochromeElectronParameters{FT<:Real,L}
    psii_fraction::FT
    psi_maximum_yield::FT
    limitation::L
    function CytochromeElectronParameters(
        psii_fraction::FT,
        psi_maximum_yield::FT,
        limitation::L,
    ) where {FT<:Real,L}
        require_continuous_scalar(
            FT, CytochromeElectronParameters, (psii_fraction, psi_maximum_yield, limitation),
        )
        isfinite(psii_fraction) && zero(FT) < psii_fraction < one(FT) ||
            throw(ArgumentError("invalid PSII fraction: coupled Z-scheme transport requires 0 < psii_fraction < 1"))
        isfinite(psi_maximum_yield) && zero(FT) < psi_maximum_yield < one(FT) ||
            throw(ArgumentError("invalid PSI maximum photochemical yield"))
        return new{FT,L}(psii_fraction, psi_maximum_yield, limitation)
    end
end
CytochromeElectronParameters(psii_fraction::Real, psi_maximum_yield::Real, limitation) =
    CytochromeElectronParameters(
        promote_compatible(CytochromeElectronParameters, psii_fraction, psi_maximum_yield)...,
        limitation,
    )
CytochromeElectronParameters(; psii_fraction, psi_maximum_yield, limitation) =
    CytochromeElectronParameters(psii_fraction, psi_maximum_yield, limitation)

"""
    cytochrome_electron_supply(parameters, total_absorbed_par, vqmax[, excitation])

Evaluate potential PSI electron supply `psi_electron_supply` [μmol e⁻ m⁻² leaf s⁻¹]
co-limited by Cytochrome b₆f capacity `vqmax`, together with the steady-state
open PSII reaction-center fraction `psii_open_fraction` (`q₂`, Johnson & Berry
2021, Eq. 22).

Single-photosystem illumination (`iszero(Q₂) != iszero(Q₁)`) lies outside the
coupled Z-scheme steady-state domain and returns `NaN` fields.
"""
@inline function cytochrome_electron_supply(
    p::CytochromeElectronParameters,
    absorbed_par,
    vqmax,
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    # Dispatch psi_excitation/psii_excitation unconditionally so unsupported
    # excitation types (e.g. PrescribedPSIIExcitation) fail via MethodError.
    q1 = psi_excitation(p, absorbed_par, excitation)
    q2 = psii_excitation(p, absorbed_par, excitation)
    OutFT = typeof(p.psi_maximum_yield + q1 + q2 + zero(vqmax))
    vqp = primal_scalar(vqmax)
    q1p, q2p = primal_scalar(q1), primal_scalar(q2)
    valid = valid_excitation(excitation, absorbed_par) && !isnan(vqp) && vqp >= zero(vqp) &&
        iszero(q2p) == iszero(q1p)
    psi_light = p.psi_maximum_yield * q1
    psi_supply = OutFT(colimit(psi_light, vqmax, p.limitation))
    psi_light_p = primal_scalar(psi_light)
    q2_open = if isinf(vqp)
        one(OutFT)
    elseif exact_zero_scalar(vqmax) || (iszero(vqp) && !(p.limitation isa SerialLimitation && !iszero(psi_light_p)))
        iszero(psi_light_p) && iszero(q2p) ? one(OutFT) : zero(OutFT)
    elseif p.limitation isa SerialLimitation
        if isinf(psi_light_p)
            zero(OutFT)
        else
            scale = max(vqmax, psi_light)
            x, y = vqmax / scale, psi_light / scale
            OutFT(x / (x + y))
        end
    else
        OutFT(clamp(one(psi_supply) - psi_supply / vqmax, zero(psi_supply), one(psi_supply)))
    end
    return checked_result(OutFT, valid, (; psi_electron_supply = psi_supply, psii_open_fraction = q2_open))
end

"""
    CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
        limitation, compensation)

Compose steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021) from immutable component structs.
"""
Base.@kwdef struct CytochromePhotosynthesis{T,E,S,L,C} <: AbstractC3Photosynthesis{L,C}
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    compensation::C
end

@inline function electron_demand(r::CytochromeRates, co2, s::CytochromeStoichiometry)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    return c1 * co2 + c2 * r.gamma_star
end

@inline function capacity_and_slope(
    r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation,
)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    dc = co2 + r.kc * (one(oxygen) + oxygen / r.ko)
    dj = c1 * co2 + c2 * r.gamma_star
    ac = r.vcmax / dc
    aj = psi_electron_supply / dj
    return colimit_slope(ac, aj, -ac / dc, -c1 * aj / dj, limitation)
end

@inline function biochemical_turnover(
    r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation, ::SignedAssimilation,
)
    capacity = assimilation_capacity(r, co2, oxygen, psi_electron_supply, s, limitation)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s,
    )
    carboxylation = capacity * co2
    photorespiratory_co2 = capacity * r.gamma_star
    linear_flow = capacity * (s.electrons_per_carboxylation * co2 +
                              s.electrons_per_photorespiratory_co2 * r.gamma_star)
    psi_flow = capacity * electron_demand(r, co2, s)
    # Direct ATP-deficit coefficients avoid Float32 cancellation in c1 - e_c and J_P700 - J_P680.
    cyclic_flow = capacity * (cef_c1 * co2 + cef_c2 * r.gamma_star)
    S = typeof(carboxylation + photorespiratory_co2 + linear_flow + psi_flow + cyclic_flow + r.respiration)
    return typed_result(S, (;
        carboxylation,
        photorespiratory_co2,
        r.respiration,
        linear_electron_flow = linear_flow,
        psi_electron_flow = psi_flow,
        cyclic_electron_flow = cyclic_flow,
    ))
end

@inline function valid_cytochrome_rates_and_stoichiometry(
    r::CytochromeRates,
    s::CytochromeStoichiometry;
    require_positive_gamma::Bool = true,
)
    vcp, vqp, rdp, kcp, kop, gsp, nlp, ncp = map(
        primal_scalar,
        (
            r.vcmax,
            r.vqmax,
            r.respiration,
            r.kc,
            r.ko,
            r.gamma_star,
            r.linear_atp_per_electron,
            r.cyclic_atp_per_electron,
        ),
    )
    zr = zero(vcp)
    zs = zero(primal_scalar(s.electrons_per_carboxylation))
    all(isfinite, (vcp, vqp, rdp, kcp, kop, gsp, nlp, ncp)) || return false
    gamma_ok = require_positive_gamma ? (gsp > zr) : (gsp >= zr)
    min(vcp, vqp, rdp, nlp) >= zr && min(kcp, kop, ncp) > zr && gamma_ok || return false
    all(
        x -> isfinite(primal_scalar(x)) && primal_scalar(x) > zs,
        (
            s.electrons_per_carboxylation,
            s.electrons_per_photorespiratory_co2,
            s.atp_per_carboxylation,
            s.atp_per_photorespiratory_co2,
        ),
    ) || return false
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s,
    )
    c1 = s.electrons_per_carboxylation + cef_c1
    c2 = s.electrons_per_photorespiratory_co2 + cef_c2
    c1p, c2p, cef1p, cef2p = map(primal_scalar, (c1, c2, cef_c1, cef_c2))
    z = zero(c1p)
    return isfinite(c1p) && isfinite(c2p) && cef1p >= z && cef2p >= z && c1p > z && c2p > z
end

@inline function valid_zscheme_supply(psi_supply, q2, eta)
    ps, q2p, ep = map(primal_scalar, (psi_supply, q2, eta))
    return isfinite(ps) && ps >= zero(ps) &&
        isfinite(q2p) && q2p >= zero(q2p) &&
        isfinite(ep) && ps <= ep * q2p
end

@inline function light_supply(p::CytochromePhotosynthesis, env, excitation)
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(
        p.electron_transport, env.absorbed_par, rates.vqmax, excitation,
    )
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    return rates, supply.psi_electron_supply,
        (; supply.psi_electron_supply, supply.psii_open_fraction, psii_excitation = q2)
end

# PSI electron cost per carboxylation is c₁ = e_c + cef_c1, so the light-limited branch
# is J_P700 (C - Γ*) / (c₁ C + c₂ Γ*) = (J_P700 / c₁) (C - Γ*) / (C + (c₂ / c₁) Γ*).
@inline function light_limited_kinetics(p::CytochromePhotosynthesis, r, psi_electron_supply)
    c1, c2 = cytochrome_carbon_coefficients(r, p.stoichiometry)
    return psi_electron_supply / c1, (c2 / c1) * r.gamma_star
end

# Linear electron transport J_P680 = J_P700 / η(C) and the Z-scheme domain at `co2`.
@inline function electron_state(p::CytochromePhotosynthesis, r, supply, co2)
    eta = psi_to_psii_electron_ratio(r, co2, p.stoichiometry)
    valid = !(iszero(primal_scalar(co2)) && iszero(primal_scalar(r.gamma_star))) &&
        valid_zscheme_supply(supply.psi_electron_supply, supply.psii_excitation, eta)
    return valid, supply.psi_electron_supply / eta,
        (; supply.psi_electron_supply, supply.psii_open_fraction, psi_to_psii_ratio = eta)
end

@inline valid_biochemistry(p::CytochromePhotosynthesis, r) =
    valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false)

@inline valid_equilibrium_biochemistry(p::CytochromePhotosynthesis, r) =
    valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry)

# Cyclic PSI flow evolves no O₂; it is reported beside the realized linear flow.
@inline turnover_electron_flows(::CytochromePhotosynthesis, turnover) =
    (; turnover.psi_electron_flow, turnover.cyclic_electron_flow)

