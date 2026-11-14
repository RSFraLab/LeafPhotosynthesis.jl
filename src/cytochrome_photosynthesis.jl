"""
    CytochromeCapacity(; site_density, turnover_rate)
    CytochromeCapacity(site_density, turnover_rate)

Parameterize Cytochrome b₆f maximum plastoquinol-oxidation capacity from active
site density `site_density` [μmol sites m⁻² leaf] and catalytic turnover law
`turnover_rate` [mol e⁻ mol⁻¹ sites s⁻¹] (Johnson & Berry 2021):
``V_{q\\max}(T) = B \\cdot k_q(T)`` [μmol e⁻ m⁻² leaf s⁻¹].
"""
Base.@kwdef struct CytochromeCapacity{FT<:Real,K<:AbstractTemperatureResponse} <: AbstractTemperatureResponse
    site_density::FT
    turnover_rate::K
    function CytochromeCapacity{FT,K}(
        fields::Vararg{Any,2},
    ) where {FT<:Real,K<:AbstractTemperatureResponse}
        c = new{FT,K}(continuous_fields(FT, CytochromeCapacity{FT,K}, fields)...)
        is_nonnegative(c.site_density) ||
            throw(ArgumentError("Cytochrome b6f site_density must be finite and nonnegative"))
        return c
    end
end
CytochromeCapacity(site_density::Real, turnover_rate::AbstractTemperatureResponse) =
    CytochromeCapacity{scalar_type(CytochromeCapacity, site_density),typeof(turnover_rate)}(
        site_density, turnover_rate,
    )

@inline function temperature_response(c::CytochromeCapacity, temperature)
    require_compatible(temperature_response, c, temperature)
    return c.site_density * temperature_response(c.turnover_rate, temperature)
end

"""
    CytochromeTemperatureResponse(; vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron)

Bundle temperature and oxygen response laws for Cytochrome b₆f C3 biochemistry
(Johnson & Berry 2021).

Every field is an [`AbstractTemperatureResponse`](@ref): `vqmax` a [`CytochromeCapacity`](@ref)
or any law returning μmol e⁻ m⁻² leaf s⁻¹, `linear_atp_per_electron` (`η_L`) and
`cyclic_atp_per_electron` (`η_C`) a [`ConstantResponse`](@ref) or a temperature law;
`gamma_star` may also be an [`AbstractCompensationPoint`](@ref).
"""
Base.@kwdef struct CytochromeTemperatureResponse{
    V<:AbstractTemperatureResponse,
    Q<:AbstractTemperatureResponse,
    R<:AbstractTemperatureResponse,
    C<:AbstractTemperatureResponse,
    O<:AbstractTemperatureResponse,
    G<:GammaStarLaw,
    NL<:AbstractTemperatureResponse,
    NC<:AbstractTemperatureResponse,
} <: AbstractLeafScalar
    vcmax::V
    vqmax::Q
    respiration::R
    kc::C
    ko::O
    gamma_star::G
    linear_atp_per_electron::NL
    cyclic_atp_per_electron::NC
    CytochromeTemperatureResponse{V,Q,R,C,O,G,NL,NC}(fields::Vararg{Any,8}) where {V,Q,R,C,O,G,NL,NC} =
        new{V,Q,R,C,O,G,NL,NC}(component_fields(CytochromeTemperatureResponse{V,Q,R,C,O,G,NL,NC}, fields)...)
end
CytochromeTemperatureResponse(
    vcmax::AbstractTemperatureResponse,
    vqmax::AbstractTemperatureResponse,
    respiration::AbstractTemperatureResponse,
    kc::AbstractTemperatureResponse,
    ko::AbstractTemperatureResponse,
    gamma_star::GammaStarLaw,
    linear_atp_per_electron::AbstractTemperatureResponse,
    cyclic_atp_per_electron::AbstractTemperatureResponse,
) = CytochromeTemperatureResponse{
    map(typeof, (vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron))...,
}(vcmax, vqmax, respiration, kc, ko, gamma_star, linear_atp_per_electron, cyclic_atp_per_electron)

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
Base.@kwdef struct CytochromeRates{FT<:Real} <: AbstractLeafScalar
    vcmax::FT
    vqmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    linear_atp_per_electron::FT
    cyclic_atp_per_electron::FT
    CytochromeRates{FT}(fields::Vararg{Any,8}) where {FT<:Real} =
        new{FT}(continuous_fields(FT, CytochromeRates{FT}, fields)...)
end
CytochromeRates(fields::Vararg{Real,8}) =
    CytochromeRates{scalar_type(CytochromeRates, fields...)}(fields...)

@inline function biochemical_rates(p::CytochromeTemperatureResponse, temperature, oxygen = nothing)
    require_compatible(biochemical_rates, p, temperature, oxygen)
    return CytochromeRates(
        temperature_response(p.vcmax, temperature),
        temperature_response(p.vqmax, temperature),
        temperature_response(p.respiration, temperature),
        temperature_response(p.kc, temperature),
        temperature_response(p.ko, temperature),
        gamma_star_response(p.gamma_star, temperature, oxygen),
        temperature_response(p.linear_atp_per_electron, temperature),
        temperature_response(p.cyclic_atp_per_electron, temperature),
    )
end

"""
    CytochromeStoichiometry{FT}(; electrons_per_carboxylation,
        electrons_per_photorespiratory_co2, atp_per_carboxylation,
        atp_per_photorespiratory_co2)
    CytochromeStoichiometry(electrons_per_carboxylation,
        electrons_per_photorespiratory_co2, atp_per_carboxylation,
        atp_per_photorespiratory_co2)

Specify C3 Calvin–Benson and photorespiratory linear-electron and ATP demands
per CO₂ event (classically `4, 8, 3, 7` in Johnson & Berry 2021, Eq. 24).
Throws `ArgumentError` unless all four are finite and positive.
"""
Base.@kwdef struct CytochromeStoichiometry{FT<:Real} <: AbstractLeafScalar
    electrons_per_carboxylation::FT
    electrons_per_photorespiratory_co2::FT
    atp_per_carboxylation::FT
    atp_per_photorespiratory_co2::FT
    function CytochromeStoichiometry{FT}(fields::Vararg{Any,4}) where {FT<:Real}
        s = new{FT}(continuous_fields(FT, CytochromeStoichiometry{FT}, fields)...)
        all(is_positive, (s.electrons_per_carboxylation, s.electrons_per_photorespiratory_co2,
            s.atp_per_carboxylation, s.atp_per_photorespiratory_co2)) ||
            throw(ArgumentError("CytochromeStoichiometry requires finite positive coefficients"))
        return s
    end
end
CytochromeStoichiometry(fields::Vararg{Real,4}) =
    CytochromeStoichiometry{scalar_type(CytochromeStoichiometry, fields...)}(fields...)

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
Base.@kwdef struct CytochromeElectronParameters{FT<:Real,L<:AbstractColimitation} <: AbstractLeafScalar
    psii_fraction::FT
    psi_maximum_yield::FT
    limitation::L
    function CytochromeElectronParameters{FT,L}(
        fields::Vararg{Any,3},
    ) where {FT<:Real,L<:AbstractColimitation}
        p = new{FT,L}(continuous_fields(FT, CytochromeElectronParameters{FT,L}, fields)...)
        is_open_fraction(p.psii_fraction) ||
            throw(ArgumentError("invalid PSII fraction: coupled Z-scheme transport requires 0 < psii_fraction < 1"))
        is_open_fraction(p.psi_maximum_yield) ||
            throw(ArgumentError("invalid PSI maximum photochemical yield"))
        return p
    end
end
CytochromeElectronParameters(
    psii_fraction::Real,
    psi_maximum_yield::Real,
    limitation::AbstractColimitation,
) =
    CytochromeElectronParameters{
        scalar_type(CytochromeElectronParameters, psii_fraction, psi_maximum_yield),typeof(limitation),
    }(psii_fraction, psi_maximum_yield, limitation)

"""
    cytochrome_electron_supply(parameters, total_absorbed_par, vqmax,
        excitation=BroadbandExcitation())

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
    require_compatible(cytochrome_electron_supply, p, absorbed_par, vqmax, excitation)
    # Dispatch psi_excitation/psii_excitation unconditionally so unsupported
    # excitation types (e.g. PrescribedPSIIExcitation) fail via MethodError.
    q1 = psi_excitation(p, absorbed_par, excitation)
    q2 = psii_excitation(p, absorbed_par, excitation)
    vqp = primal_scalar(vqmax)
    q1p, q2p = primal_scalar(q1), primal_scalar(q2)
    valid = valid_excitation(excitation, absorbed_par) && !isnan(vqp) && vqp >= zero(vqp) &&
        iszero(q2p) == iszero(q1p)
    psi_light = p.psi_maximum_yield * q1
    colimited = colimit(psi_light, vqmax, p.limitation)
    OutFT = promote_type(typeof(psi_light + q2 + zero(vqmax)), typeof(colimited))
    psi_supply = OutFT(colimited)
    q2_open = open_psii_fraction(OutFT, p.limitation, vqmax, psi_light, psi_supply, q2p)
    return checked_result(valid, (; psi_electron_supply = psi_supply, psii_open_fraction = q2_open))
end

# Open PSII fraction q₂ (Johnson & Berry 2021, Eq. 22) at Cytochrome capacity `vqmax`, PSI
# light supply `psi_light`, co-limited PSI supply `psi_supply` and primal PSII excitation
# `q2p`. Without capacity the centres are open only in the dark.
@inline function open_psii_fraction(::Type{OutFT}, ::SerialLimitation, vqmax, psi_light, psi_supply, q2p) where {OutFT}
    vqp, psi_light_p = primal_scalar(vqmax), primal_scalar(psi_light)
    isinf(vqp) && return one(OutFT)
    if exact_zero_scalar(vqmax) || (iszero(vqp) && iszero(psi_light_p))
        return iszero(psi_light_p) && iszero(q2p) ? one(OutFT) : zero(OutFT)
    end
    isinf(psi_light_p) && return zero(OutFT)
    scale = max(vqmax, psi_light)
    x, y = vqmax / scale, psi_light / scale
    return OutFT(x / (x + y))
end
@inline function open_psii_fraction(::Type{OutFT}, ::AbstractColimitation, vqmax, psi_light, psi_supply, q2p) where {OutFT}
    vqp, psi_light_p = primal_scalar(vqmax), primal_scalar(psi_light)
    isinf(vqp) && return one(OutFT)
    iszero(vqp) && return iszero(psi_light_p) && iszero(q2p) ? one(OutFT) : zero(OutFT)
    return OutFT(primal_clamp(one(psi_supply) - psi_supply / vqmax, zero(psi_supply), one(psi_supply)))
end

"""
    CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
        limitation, assimilation_sign)
    CytochromePhotosynthesis(temperature, electron_transport, stoichiometry,
        limitation, assimilation_sign)

Compose steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021) from immutable components: `temperature`
([`CytochromeTemperatureResponse`](@ref)), `electron_transport`
([`CytochromeElectronParameters`](@ref)), `stoichiometry` ([`CytochromeStoichiometry`](@ref)),
the carbon `limitation` ([`AbstractColimitation`](@ref)) and `assimilation_sign`
([`AbstractAssimilationSign`](@ref)). Components of different binary precisions are a
`MethodError`.
"""
Base.@kwdef struct CytochromePhotosynthesis{
    T<:AbstractLeafScalar,
    E<:AbstractLeafScalar,
    S<:AbstractLeafScalar,
    L<:AbstractColimitation,
    C<:AbstractAssimilationSign,
} <: AbstractC3Photosynthesis{L,C}
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    assimilation_sign::C
    CytochromePhotosynthesis{T,E,S,L,C}(fields::Vararg{Any,5}) where {T,E,S,L,C} =
        new{T,E,S,L,C}(component_fields(CytochromePhotosynthesis{T,E,S,L,C}, fields)...)
end
CytochromePhotosynthesis(
    temperature::AbstractLeafScalar,
    electron_transport::AbstractLeafScalar,
    stoichiometry::AbstractLeafScalar,
    limitation::AbstractColimitation,
    assimilation_sign::AbstractAssimilationSign,
) = CytochromePhotosynthesis{map(typeof, (temperature, electron_transport, stoichiometry,
    limitation, assimilation_sign))...}(temperature, electron_transport, stoichiometry, limitation,
    assimilation_sign)

# Cytochrome models whose PSI supply is co-limited in series with Cytochrome b₆f, the only
# electron co-limitation for which the state-transition and dark-partition closed forms hold.
const SerialCytochromePhotosynthesis =
    CytochromePhotosynthesis{<:Any,<:CytochromeElectronParameters{<:Real,SerialLimitation}}

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

@inline valid_zscheme_supply(psi_supply, q2, eta) =
    is_nonnegative(psi_supply) && is_nonnegative(q2) && is_finite(eta) &&
    primal_scalar(psi_supply) <= primal_scalar(eta) * primal_scalar(q2)

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

# Rate domain: finite rates, nonnegative capacities, respiration, Γ* and η_L, positive
# Michaelis constants and η_C, and an ATP balance that needs no negative cyclic flow
# (cef ≥ 0, c₁, c₂ > 0). The stoichiometry is checked by its constructor.
@inline function valid_biochemistry(p::CytochromePhotosynthesis, r)
    s = p.stoichiometry
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s,
    )
    c1 = s.electrons_per_carboxylation + cef_c1
    c2 = s.electrons_per_photorespiratory_co2 + cef_c2
    return all(is_nonnegative, (r.vcmax, r.vqmax, r.respiration, r.gamma_star,
            r.linear_atp_per_electron)) &&
        all(is_positive, (r.kc, r.ko, r.cyclic_atp_per_electron)) &&
        is_nonnegative(cef_c1) && is_nonnegative(cef_c2) && is_positive(c1) && is_positive(c2)
end

# Cyclic PSI flow evolves no O₂; it is reported beside the realized linear flow.
@inline turnover_electron_flows(::CytochromePhotosynthesis, turnover) =
    (; turnover.psi_electron_flow, turnover.cyclic_electron_flow)

