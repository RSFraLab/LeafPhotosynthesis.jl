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
    if !valid_excitation(excitation, absorbed_par) || isnan(vqp) || vqp < zero(vqp) ||
       (iszero(q2p) != iszero(q1p))
        n = OutFT(NaN)
        return (; psi_electron_supply = n, psii_open_fraction = n)
    end
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
    return (; psi_electron_supply = psi_supply, psii_open_fraction = q2_open)
end

"""
    CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
        limitation, compensation)

Compose steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021) from immutable component structs.
"""
Base.@kwdef struct CytochromePhotosynthesis{T,E,S,L,C} <: AbstractPhotosynthesis
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    compensation::C
end

@inline function assimilation_capacity(
    r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation,
)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    rubisco = r.vcmax / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    cytochrome = psi_electron_supply / (c1 * co2 + c2 * r.gamma_star)
    return colimit(rubisco, cytochrome, limitation)
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

@inline function carbon_assimilation(
    r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation, compensation,
)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    excess = co2_excess(co2, r.gamma_star, compensation)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    cytochrome = psi_electron_supply * excess / (c1 * co2 + c2 * r.gamma_star)
    assimilation = assimilation_capacity(r, co2, oxygen, psi_electron_supply, s, limitation) * excess
    net = assimilation - r.respiration
    S = typeof(rubisco + cytochrome + net)
    return (;
        rubisco_limited = S(rubisco),
        electron_limited = S(cytochrome),
        assimilation = S(assimilation),
        respiration = S(r.respiration),
        net_assimilation = S(net),
    )
end

@inline function biochemical_turnover(
    r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation, ::SignedAssimilation,
)
    capacity = assimilation_capacity(r, co2, oxygen, psi_electron_supply, s, limitation)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s,
    )
    c1 = s.electrons_per_carboxylation + cef_c1
    c2 = s.electrons_per_photorespiratory_co2 + cef_c2
    carboxylation = capacity * co2
    photorespiratory_co2 = capacity * r.gamma_star
    linear_flow = capacity * (s.electrons_per_carboxylation * co2 +
                              s.electrons_per_photorespiratory_co2 * r.gamma_star)
    psi_flow = capacity * (c1 * co2 + c2 * r.gamma_star)
    # Direct ATP-deficit coefficients avoid Float32 cancellation in c1 - e_c and J_P700 - J_P680.
    cyclic_flow = capacity * (cef_c1 * co2 + cef_c2 * r.gamma_star)
    S = typeof(carboxylation + photorespiratory_co2 + linear_flow + psi_flow + cyclic_flow + r.respiration)
    return (;
        carboxylation = S(carboxylation),
        photorespiratory_co2 = S(photorespiratory_co2),
        respiration = S(r.respiration),
        linear_electron_flow = S(linear_flow),
        psi_electron_flow = S(psi_flow),
        cyclic_electron_flow = S(cyclic_flow),
    )
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

@inline function photosynthesize(
    p::CytochromePhotosynthesis,
    env::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {FT<:Real}
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(
        p.electron_transport, env.absorbed_par, rates.vqmax, excitation,
    )
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    OutFT = typeof(
        rates.vcmax + p.stoichiometry.electrons_per_carboxylation +
        env.co2 + env.oxygen + supply.psi_electron_supply + q2,
    )
    n = OutFT(NaN)
    co2_p, o2_p = primal_scalar(env.co2), primal_scalar(env.oxygen)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(co2_p) && co2_p >= zero(co2_p) && isfinite(o2_p) && o2_p >= zero(o2_p))
        return (; invalid_carbon(OutFT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    if !valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) ||
       (iszero(co2_p) && iszero(primal_scalar(rates.gamma_star)))
        return (; invalid_carbon(OutFT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    eta = psi_to_psii_electron_ratio(rates, env.co2, p.stoichiometry)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; invalid_carbon(OutFT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    carbon = carbon_assimilation(
        rates, env.co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation,
    )
    return (;
        rubisco_limited = OutFT(carbon.rubisco_limited),
        electron_limited = OutFT(carbon.electron_limited),
        assimilation = OutFT(carbon.assimilation),
        respiration = OutFT(carbon.respiration),
        net_assimilation = OutFT(carbon.net_assimilation),
        electron_transport = OutFT(supply.psi_electron_supply / eta),
        psi_electron_supply = OutFT(supply.psi_electron_supply),
        psii_open_fraction = OutFT(supply.psii_open_fraction),
        psi_to_psii_ratio = OutFT(eta),
    )
end

@inline function biochemical_turnover(
    p::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L,FT<:Real}
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(
        p.electron_transport, env.absorbed_par, r.vqmax, excitation,
    )
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    OutFT = typeof(
        r.vcmax + p.stoichiometry.electrons_per_carboxylation +
        env.co2 + env.oxygen + supply.psi_electron_supply + q2,
    )
    n = OutFT(NaN)
    co2_p, o2_p = primal_scalar(env.co2), primal_scalar(env.oxygen)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(co2_p) && co2_p >= zero(co2_p) && isfinite(o2_p) && o2_p >= zero(o2_p))
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    if !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false) ||
       (iszero(co2_p) && iszero(primal_scalar(r.gamma_star)))
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    eta = psi_to_psii_electron_ratio(r, env.co2, p.stoichiometry)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    turnover = biochemical_turnover(
        r, env.co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation,
    )
    return (;
        carboxylation = OutFT(turnover.carboxylation),
        photorespiratory_co2 = OutFT(turnover.photorespiratory_co2),
        respiration = OutFT(turnover.respiration),
        linear_electron_flow = OutFT(turnover.linear_electron_flow),
        psi_electron_flow = OutFT(turnover.psi_electron_flow),
        cyclic_electron_flow = OutFT(turnover.cyclic_electron_flow),
        electron_supply = OutFT(supply.psi_electron_supply / eta),
        psi_electron_supply = OutFT(supply.psi_electron_supply),
        psii_open_fraction = OutFT(supply.psii_open_fraction),
    )
end

@inline function oxygen_exchange(
    p::CoreOxygenExchange,
    model::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,L}
    # Gross O₂ evolution comes strictly from PSII water oxidation (realized linear_electron_flow);
    # cyclic PSI electron flow evolves zero O₂.
    turnover = biochemical_turnover(model, env, excitation)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    inputs = OxygenFluxInputs(turnover.linear_electron_flow, oxygenation, turnover.respiration)
    ox = oxygen_exchange(p, inputs)
    OutFT = typeof(ox.net_release + turnover.carboxylation)
    return (;
        gross_production = OutFT(ox.gross_production),
        photorespiratory_uptake = OutFT(ox.photorespiratory_uptake),
        respiratory_uptake = OutFT(ox.respiratory_uptake),
        net_release = OutFT(ox.net_release),
        carboxylation = OutFT(turnover.carboxylation),
        rubisco_oxygenation = OutFT(oxygenation),
        linear_electron_flow = OutFT(turnover.linear_electron_flow),
        psi_electron_flow = OutFT(turnover.psi_electron_flow),
        cyclic_electron_flow = OutFT(turnover.cyclic_electron_flow),
    )
end

@inline function gas_exchange(
    p::CytochromePhotosynthesis{T,E,S,MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment{FT},
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation(),
) where {T,E,S,FT<:Real}
    conductance = inv(inv(d.boundary_co2) + d.diffusivity_ratio / state.conductance + inv(d.mesophyll_co2))
    g = conductance * d.mole_to_micromole / env.pressure
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    OutFT = typeof(
        r.vcmax + p.stoichiometry.electrons_per_carboxylation +
        env.co2 + g + conductance + supply.psi_electron_supply + q2,
    )
    n = OutFT(NaN)
    t_p, co2_p, o2_p, p_p = map(
        primal_scalar,
        (env.temperature, env.co2, env.oxygen, env.pressure),
    )
    gs_p, g_p = primal_scalar(state.conductance), primal_scalar(g)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(t_p) && t_p > zero(t_p) &&
         isfinite(co2_p) && co2_p >= zero(co2_p) && isfinite(o2_p) && o2_p >= zero(o2_p) &&
         isfinite(p_p) && p_p > zero(p_p) &&
         isfinite(gs_p) && gs_p > zero(gs_p)) ||
       !(isfinite(g_p) && g_p > zero(g_p)) ||
       !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false)
        return (; carbon = invalid_carbon(OutFT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    c1, c2 = cytochrome_carbon_coefficients(r, p.stoichiometry)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    cytochrome_co2 = diffusion_co2_root(
        supply.psi_electron_supply / c1,
        (c2 / c1) * r.gamma_star, r.gamma_star, r.respiration, env.co2, g,
    )
    above_compensation = g * (env.co2 - r.gamma_star) + r.respiration >= zero(g)
    co2 = OutFT(above_compensation ? max(rubisco_co2, cytochrome_co2) : min(rubisco_co2, cytochrome_co2))
    eta = psi_to_psii_electron_ratio(r, co2, p.stoichiometry)
    if !isfinite(primal_scalar(co2)) || (iszero(primal_scalar(co2)) && iszero(primal_scalar(r.gamma_star))) ||
       !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; carbon = invalid_carbon(OutFT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    carbon_raw = carbon_assimilation(
        r, co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation,
    )
    carbon = (;
        rubisco_limited = OutFT(carbon_raw.rubisco_limited),
        electron_limited = OutFT(carbon_raw.electron_limited),
        assimilation = OutFT(carbon_raw.assimilation),
        respiration = OutFT(carbon_raw.respiration),
        net_assimilation = OutFT(carbon_raw.net_assimilation),
    )
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    surface_co2 = OutFT(env.co2 - pressure_flux / d.boundary_co2)
    intercellular_co2 = OutFT(co2 + pressure_flux / d.mesophyll_co2)
    return (;
        carbon,
        electron_transport = OutFT(supply.psi_electron_supply / eta),
        psi_electron_supply = OutFT(supply.psi_electron_supply),
        psii_open_fraction = OutFT(supply.psii_open_fraction),
        surface_co2,
        intercellular_co2,
        biochemical_co2 = co2,
        conductance = OutFT(conductance),
        diffusion_residual = OutFT(g * (env.co2 - co2) - carbon.net_assimilation),
    )
end

@inline function cytochrome_equilibrium_result(
    ::Type{OutFT}, co2, p, r, psi_supply, q2, env, d, g0, gain,
    width, iterations, status,
) where {OutFT}
    co2_s = OutFT(co2)
    eta = psi_to_psii_electron_ratio(r, co2_s, p.stoichiometry)
    valid_zscheme_supply(psi_supply, q2, eta) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    e = equilibrium_residual(co2_s, p, r, psi_supply, env, d, g0, gain)
    carbon_raw = carbon_assimilation(
        r, co2_s, env.oxygen, psi_supply,
        p.stoichiometry, p.limitation, p.compensation,
    )
    carbon = (;
        rubisco_limited = OutFT(carbon_raw.rubisco_limited),
        electron_limited = OutFT(carbon_raw.electron_limited),
        assimilation = OutFT(carbon_raw.assimilation),
        respiration = OutFT(carbon_raw.respiration),
        net_assimilation = OutFT(carbon_raw.net_assimilation),
    )
    conductance = OutFT(inv(e.resistance))
    flux_residual = OutFT(e.residual * conductance * d.mole_to_micromole / env.pressure)
    intercellular = OutFT(co2_s + carbon.net_assimilation * env.pressure / (d.mole_to_micromole * d.mesophyll_co2))
    return EquilibriumSolution{OutFT,typeof(carbon)}(
        carbon, OutFT(psi_supply / eta), co2_s, intercellular, OutFT(e.surface), OutFT(e.gs),
        conductance, OutFT(e.residual), flux_residual, OutFT(width), iterations, status, status == EquilibriumConverged,
    )
end

@inline function equilibrium_gas_exchange(
    p::CytochromePhotosynthesis,
    stomata::SteadyStateStomata,
    env::GasExchangeEnvironment{FT},
    d::SeriesDiffusion,
    solver::EquilibriumSolver,
    excitation::AbstractExcitation = BroadbandExcitation();
    initial_co2 = nothing,
) where {FT<:Real}
    g0, gain = equilibrium_coefficients(stomata.target, env)
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    psi_supply = supply.psi_electron_supply
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    OutFT = typeof(
        r.vcmax + p.stoichiometry.electrons_per_carboxylation +
        env.co2 + d.diffusivity_ratio + g0 + gain + psi_supply + q2 + zero(solver.co2_atol),
    )
    valid_equilibrium_inputs(env, d, stomata.target) || return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    valid_excitation(excitation, env.absorbed_par) || return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    isfinite(gain) || return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    isfinite(psi_supply) && primal_scalar(psi_supply) >= zero(primal_scalar(psi_supply)) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidInput)
    resistance0 = inv(d.boundary_co2) + d.diffusivity_ratio / g0 + inv(d.mesophyll_co2)
    if iszero(primal_scalar(psi_supply)) || iszero(primal_scalar(r.vcmax))
        x = OutFT(env.co2 + r.respiration * env.pressure / d.mole_to_micromole * resistance0)
        e = equilibrium_residual(x, p, r, psi_supply, env, d, g0, gain)
        status = equilibrium_resolved(e, x, env, d, solver) ? EquilibriumConverged : EquilibriumPrecisionLimit
        x_pol = polish_equilibrium_tangent(
            x, zero(OutFT), x + max(one(OutFT), x), p, r, psi_supply, env, d, g0, gain, status,
        )
        return cytochrome_equilibrium_result(
            OutFT, x_pol, p, r, psi_supply, q2, env, d, g0, gain, zero(OutFT), 0, status,
        )
    end
    lo = zero(OutFT)
    elo = equilibrium_residual(lo, p, r, psi_supply, env, d, g0, gain)
    hi = OutFT(env.co2 + max(-elo.assimilation, zero(elo.assimilation)) * env.pressure / d.mole_to_micromole * resistance0)
    ehi = equilibrium_residual(hi, p, r, psi_supply, env, d, g0, gain)
    if equilibrium_resolved(elo, lo, env, d, solver)
        x_lo = polish_equilibrium_tangent(lo, lo, hi, p, r, psi_supply, env, d, g0, gain, EquilibriumConverged)
        return cytochrome_equilibrium_result(
            OutFT, x_lo, p, r, psi_supply, q2, env, d, g0, gain, hi - lo, 0, EquilibriumConverged,
        )
    end
    if equilibrium_resolved(ehi, hi, env, d, solver)
        x_hi = polish_equilibrium_tangent(hi, lo, hi, p, r, psi_supply, env, d, g0, gain, EquilibriumConverged)
        return cytochrome_equilibrium_result(
            OutFT, x_hi, p, r, psi_supply, q2, env, d, g0, gain, hi - lo, 0, EquilibriumConverged,
        )
    end
    isfinite(hi) && isfinite(elo.residual) && isfinite(ehi.residual) &&
        elo.residual >= zero(OutFT) && ehi.residual <= zero(OutFT) ||
        return invalid_equilibrium(OutFT, EquilibriumInvalidBracket)
    lo_orig, hi_orig = lo, hi
    x = if isnothing(initial_co2)
        midpoint(lo, hi)
    elseif initial_co2 isa Real && compatible_scalar(typeof(initial_co2), OutFT)
        x0 = primal_scalar(initial_co2)
        isfinite(x0) && lo < x0 < hi ? OutFT(x0) : midpoint(lo, hi)
    else
        midpoint(lo, hi)
    end
    best, best_error = lo, OutFT(abs(elo.residual))
    if ehi.physical && abs(ehi.residual) < best_error
        best, best_error = hi, OutFT(abs(ehi.residual))
    end
    for iteration in 1:solver.max_iterations
        e = equilibrium_residual(x, p, r, psi_supply, env, d, g0, gain)
        isfinite(e.residual) || return invalid_equilibrium(OutFT, EquilibriumInvalidBracket)
        if e.physical && abs(e.residual) < best_error
            best, best_error = x, OutFT(abs(e.residual))
        end
        if equilibrium_resolved(e, x, env, d, solver)
            x_pol = polish_equilibrium_tangent(
                x, lo_orig, hi_orig, p, r, psi_supply, env, d, g0, gain, EquilibriumConverged,
            )
            return cytochrome_equilibrium_result(
                OutFT, x_pol, p, r, psi_supply, q2, env, d, g0, gain,
                hi - lo, iteration, EquilibriumConverged,
            )
        end
        if e.residual > zero(OutFT)
            lo = x
        else
            hi = x
        end
        candidate = equilibrium_step(solver.method, x, e, lo, hi, iteration)
        if !(lo < candidate < hi)
            return cytochrome_equilibrium_result(
                OutFT, best, p, r, psi_supply, q2, env, d, g0, gain,
                hi - lo, iteration, EquilibriumPrecisionLimit,
            )
        end
        x = candidate
    end
    return cytochrome_equilibrium_result(
        OutFT, best, p, r, psi_supply, q2, env, d, g0, gain,
        hi - lo, solver.max_iterations, EquilibriumIterationLimit,
    )
end

@inline function assimilation_sensitivity(
    p::CytochromePhotosynthesis,
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    exchange = gas_exchange(p, env, state, d),
    excitation::AbstractExcitation = BroadbandExcitation(),
)
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    OutFT = typeof(
        r.vcmax + p.stoichiometry.electrons_per_carboxylation +
        env.co2 + state.conductance + d.diffusivity_ratio +
        exchange.biochemical_co2 + supply.psi_electron_supply + q2,
    )
    n = OutFT(NaN)
    if !valid_excitation(excitation, env.absorbed_par) || !isfinite(exchange.biochemical_co2) ||
       !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    eta = psi_to_psii_electron_ratio(r, exchange.biochemical_co2, p.stoichiometry)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    slope = carbon_assimilation_derivative(
        r, exchange.biochemical_co2, env.oxygen,
        supply.psi_electron_supply, p.stoichiometry, p.limitation, p.compensation,
    )
    conductance = exchange.conductance
    conversion = d.mole_to_micromole / env.pressure
    supply_slope = conductance * conversion
    dg_dgs = (conductance / state.conductance)^2 * d.diffusivity_ratio * conversion
    dco2_dgs = OutFT(dg_dgs * (env.co2 - exchange.biochemical_co2) / (slope + supply_slope))
    derivative = OutFT(slope * dco2_dgs)
    return (; dco2_dgs, dnet_assimilation_dgs = derivative, dassimilation_dgs = derivative)
end

@inline assimilation_sensitivity(
    p::CytochromePhotosynthesis,
    env::GasExchangeEnvironment,
    state::StomatalState,
    d::SeriesDiffusion,
    excitation::AbstractExcitation,
) = assimilation_sensitivity(
    p, env, state, d,
    gas_exchange(p, env, state, d, excitation), excitation,
)
