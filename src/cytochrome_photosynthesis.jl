"""
    CytochromeCapacity(; site_density, turnover_rate)

Cytochrome b₆f maximum plastoquinol-oxidation capacity parameterised from
active-site density `site_density` [μmol sites m⁻² leaf] and catalytic turnover
response `turnover_rate` [mol e⁻ mol⁻¹ sites s⁻¹], as in Johnson & Berry (2021):
`Vqmax(T) = site_density * kq(T)` [μmol e⁻ m⁻² leaf s⁻¹].
"""
struct CytochromeCapacity{FT<:AbstractFloat,K}
    site_density::FT
    turnover_rate::K
    function CytochromeCapacity(site_density::FT, turnover_rate::K) where {FT<:AbstractFloat,K}
        isfinite(site_density) && site_density >= zero(FT) ||
            throw(ArgumentError("Cytochrome b6f site_density must be finite and nonnegative"))
        new{FT,K}(site_density, turnover_rate)
    end
end
CytochromeCapacity(; site_density, turnover_rate) =
    CytochromeCapacity(site_density, turnover_rate)

@inline temperature_response(c::CytochromeCapacity, temperature) =
    c.site_density * temperature_response(c.turnover_rate, temperature)

"""
    CytochromeTemperatureResponse(; vcmax, vqmax, respiration, kc, ko, gamma_star,
        linear_atp_per_electron, cyclic_atp_per_electron)

Temperature/oxygen laws for Cytochrome b₆f C3 biochemistry (Johnson & Berry 2021).
`vqmax` may be a `CytochromeCapacity` or any temperature-response law returning
μmol e⁻ m⁻² leaf s⁻¹. `gamma_star` may be a temperature law or an
`AbstractCompensationPoint`. `linear_atp_per_electron` (η_L) and
`cyclic_atp_per_electron` (η_C) use `ConstantResponse` or temperature laws.
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

Evaluated Cytochrome b₆f biochemical rates and coupling efficiencies at leaf
temperature [K] and biochemical-site O₂ [Pa].
- `vcmax`, `respiration`: μmol CO₂ m⁻² leaf s⁻¹
- `vqmax`: μmol electrons m⁻² leaf s⁻¹
- `kc`, `ko`, `gamma_star`: Pa
- `linear_atp_per_electron` (η_L), `cyclic_atp_per_electron` (η_C): mol ATP (mol e⁻)⁻¹
"""
Base.@kwdef struct CytochromeRates{FT}
    vcmax::FT
    vqmax::FT
    respiration::FT
    kc::FT
    ko::FT
    gamma_star::FT
    linear_atp_per_electron::FT
    cyclic_atp_per_electron::FT
end

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

Explicit C3 Calvin–Benson and photorespiratory electron and ATP demands per CO₂
event (classically `4, 8, 3, 7` in Johnson & Berry 2021, Eq. 24):
- linear electron demand = `(e_c * Cc + e_Γ * Γ*) * Vc / Cc`
- ATP demand = `(a_c * Cc + a_Γ * Γ*) * Vc / Cc`

Balancing linear ATP synthesis (`η_L * J_P680`) and PSI cyclic photophosphorylation
(`η_C * (J_P700 - J_P680)`) yields `J_P700 = (c1 * Cc + c2 * Γ*) * Vc / Cc`,
preserving an exact hyperbolic `Cc` dependence with zero hidden constants.
"""
Base.@kwdef struct CytochromeStoichiometry{FT}
    electrons_per_carboxylation::FT
    electrons_per_photorespiratory_co2::FT
    atp_per_carboxylation::FT
    atp_per_photorespiratory_co2::FT
end

"""
    cytochrome_cyclic_coefficients(linear_atp_per_electron, cyclic_atp_per_electron,
        stoichiometry::CytochromeStoichiometry)

Cyclic PSI electron requirements `(cef_c1, cef_c2)` per carboxylation and per
photorespiratory CO₂ release, computed directly from the ATP deficit before adding
linear electron demand so `a_c == η_L * e_c` evaluates to exact `0` in `Float32`.
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
    cytochrome_carbon_coefficients(rates::CytochromeRates, stoichiometry::CytochromeStoichiometry)

Effective PSI electron costs `(c1, c2)` multiplying `Cc` and `Γ*` in
`η(Cc, Γ*) * (e_c * Cc + e_Γ * Γ*) = c1 * Cc + c2 * Γ*`. For classical C3
stoichiometry (`4, 8, 3, 7`) with `η_L = 0.75` and `η_C = 1.0`, `(c1, c2) = (4, 9)`
because each photorespiratory CO₂ release requires one cyclic PSI turnover to supply
its seventh ATP.
"""
@inline function cytochrome_carbon_coefficients(r::CytochromeRates, s::CytochromeStoichiometry)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s)
    return s.electrons_per_carboxylation + cef_c1,
           s.electrons_per_photorespiratory_co2 + cef_c2
end

"""
    psi_to_psii_electron_ratio(rates::CytochromeRates, co2, stoichiometry::CytochromeStoichiometry)

Steady-state ratio `η = J_P700 / J_P680 = 1 + J_CEF / J_P680` required to balance
Calvin–Benson and photorespiratory ATP/NADPH stoichiometry at biochemical-site
`co2` [Pa]. Smooth and finite at `co2 = Γ* = 0`.
"""
@inline function psi_to_psii_electron_ratio(r::CytochromeRates, co2, s::CytochromeStoichiometry)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    lef_denom = s.electrons_per_carboxylation * co2 +
                s.electrons_per_photorespiratory_co2 * r.gamma_star
    iszero(lef_denom) && return c1 / s.electrons_per_carboxylation
    return (c1 * co2 + c2 * r.gamma_star) / lef_denom
end

"""
    CytochromeElectronParameters(; psii_fraction, psi_maximum_yield, limitation)

Broadband PSII photon partition `psii_fraction` (`β ∈ (0, 1)`; coupled Z-scheme
electron transport requires non-zero excitation of both PSII and PSI), open-center PSI
maximum photochemical quantum yield `psi_maximum_yield` (`Φ_P1,max ∈ (0, 1)`),
and rate-limitation law `limitation` (`SerialLimitation()` for Johnson & Berry 2021, Eq. 20).
"""
struct CytochromeElectronParameters{FT,L}
    psii_fraction::FT
    psi_maximum_yield::FT
    limitation::L
end

function CytochromeElectronParameters(; psii_fraction::FT, psi_maximum_yield::FT,
    limitation) where {FT}
    isfinite(psii_fraction) && zero(FT) < psii_fraction < one(FT) ||
        throw(ArgumentError("invalid PSII fraction: coupled Z-scheme transport requires 0 < psii_fraction < 1"))
    isfinite(psi_maximum_yield) && zero(FT) < psi_maximum_yield < one(FT) ||
        throw(ArgumentError("invalid PSI maximum photochemical yield"))
    return CytochromeElectronParameters(psii_fraction, psi_maximum_yield, limitation)
end

"""
    cytochrome_electron_supply(parameters, total_absorbed_par, vqmax[, excitation])

Evaluate potential PSI electron supply `psi_electron_supply` [μmol e⁻ m⁻² leaf s⁻¹]
co-limited by Cytochrome b₆f capacity `vqmax`, together with the steady-state
open PSII reaction-center fraction `psii_open_fraction` (`q₂ = 1 - J_P700,avail / Vqmax`,
Johnson & Berry 2021, Eq. 22). For `SerialLimitation`, `q₂` is evaluated via
scale-normalized `Vqmax / (Vqmax + Φ_P1,max Q₁)` to avoid both high-light subtractive
cancellation and `Float32` overflow. Single-photosystem illumination
(`iszero(Q₂) != iszero(Q₁)`) is outside the coupled Z-scheme steady-state domain
and returns `NaN`.
"""
@inline function cytochrome_electron_supply(p::CytochromeElectronParameters, absorbed_par,
    vqmax, excitation::AbstractExcitation = BroadbandExcitation())
    if !valid_excitation(excitation, absorbed_par) || isnan(vqmax) || vqmax < zero(vqmax)
        n = oftype(absorbed_par + zero(vqmax), NaN)
        return (; psi_electron_supply = n, psii_open_fraction = n)
    end
    q1 = psi_excitation(p, absorbed_par, excitation)
    q2 = psii_excitation(p, absorbed_par, excitation)
    if iszero(q2) != iszero(q1)
        n = oftype(absorbed_par + zero(vqmax), NaN)
        return (; psi_electron_supply = n, psii_open_fraction = n)
    end
    psi_light = p.psi_maximum_yield * q1
    psi_supply = colimit(psi_light, vqmax, p.limitation)
    q2_open = if isinf(vqmax)
        one(psi_supply)
    elseif iszero(vqmax)
        iszero(psi_light) && iszero(q2) ? one(psi_supply) : zero(psi_supply)
    elseif p.limitation isa SerialLimitation
        if isinf(psi_light)
            zero(psi_supply)
        else
            scale = max(vqmax, psi_light)
            x, y = vqmax / scale, psi_light / scale
            x / (x + y)
        end
    else
        clamp(one(psi_supply) - psi_supply / vqmax, zero(psi_supply), one(psi_supply))
    end
    return (; psi_electron_supply = psi_supply, psii_open_fraction = q2_open)
end

"""
    CytochromePhotosynthesis(; temperature, electron_transport, stoichiometry,
        limitation, compensation)

Composable steady-state C3 photosynthesis governed by Cytochrome b₆f plastoquinol
oxidation, PSI photochemistry, and coupled linear/cyclic ATP–NADPH stoichiometry
(Johnson & Berry 2021). Contains immutable parameters only, with no mutable
scratch arrays or embedded physical constants.
"""
Base.@kwdef struct CytochromePhotosynthesis{T,E,S,L,C} <: AbstractPhotosynthesis
    temperature::T
    electron_transport::E
    stoichiometry::S
    limitation::L
    compensation::C
end

@inline function assimilation_capacity(r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    rubisco = r.vcmax / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    cytochrome = psi_electron_supply / (c1 * co2 + c2 * r.gamma_star)
    return colimit(rubisco, cytochrome, limitation)
end

@inline function capacity_and_slope(r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    dc = co2 + r.kc * (one(oxygen) + oxygen / r.ko)
    dj = c1 * co2 + c2 * r.gamma_star
    ac = r.vcmax / dc
    aj = psi_electron_supply / dj
    return colimit_slope(ac, aj, -ac / dc, -c1 * aj / dj, limitation)
end

"""
    carbon_assimilation(rates::CytochromeRates, co2, oxygen, psi_electron_supply,
        stoichiometry::CytochromeStoichiometry, limitation, compensation)

Cytochrome-coupled C3 carbon assimilation [μmol CO₂ m⁻² leaf s⁻¹] at prescribed
biochemical-site `co2` and `oxygen` [Pa], given available `psi_electron_supply`
[μmol e⁻ m⁻² leaf s⁻¹]. `electron_limited` denotes the Cytochrome b₆f / PSI
electron-transport-limited pre-respiration rate.
"""
@inline function carbon_assimilation(r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation, compensation)
    c1, c2 = cytochrome_carbon_coefficients(r, s)
    excess = co2_excess(co2, r.gamma_star, compensation)
    rubisco = r.vcmax * excess / (co2 + r.kc * (one(oxygen) + oxygen / r.ko))
    cytochrome = psi_electron_supply * excess / (c1 * co2 + c2 * r.gamma_star)
    assimilation = assimilation_capacity(r, co2, oxygen, psi_electron_supply, s, limitation) * excess
    return (; rubisco_limited = rubisco, electron_limited = cytochrome,
        assimilation, respiration = r.respiration,
        net_assimilation = assimilation - r.respiration)
end

@inline function biochemical_turnover(r::CytochromeRates, co2, oxygen, psi_electron_supply,
    s::CytochromeStoichiometry, limitation, ::SignedAssimilation)
    capacity = assimilation_capacity(r, co2, oxygen, psi_electron_supply, s, limitation)
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s)
    c1 = s.electrons_per_carboxylation + cef_c1
    c2 = s.electrons_per_photorespiratory_co2 + cef_c2
    linear_flow = capacity * (s.electrons_per_carboxylation * co2 +
                              s.electrons_per_photorespiratory_co2 * r.gamma_star)
    psi_flow = capacity * (c1 * co2 + c2 * r.gamma_star)
    # Direct ATP-deficit coefficients avoid Float32 cancellation in c1 - e_c and J_P700 - J_P680.
    cyclic_flow = capacity * (cef_c1 * co2 + cef_c2 * r.gamma_star)
    return (; carboxylation = capacity * co2,
        photorespiratory_co2 = capacity * r.gamma_star,
        respiration = r.respiration,
        linear_electron_flow = linear_flow,
        psi_electron_flow = psi_flow,
        cyclic_electron_flow = cyclic_flow)
end

@inline function valid_cytochrome_rates_and_stoichiometry(r::CytochromeRates,
    s::CytochromeStoichiometry; require_positive_gamma::Bool = true)
    zr = zero(r.vcmax)
    zs = zero(s.electrons_per_carboxylation)
    all(isfinite, (r.vcmax, r.vqmax, r.respiration, r.kc, r.ko, r.gamma_star,
        r.linear_atp_per_electron, r.cyclic_atp_per_electron)) || return false
    gamma_ok = require_positive_gamma ? (r.gamma_star > zr) : (r.gamma_star >= zr)
    min(r.vcmax, r.vqmax, r.respiration, r.linear_atp_per_electron) >= zr &&
        min(r.kc, r.ko, r.cyclic_atp_per_electron) > zr && gamma_ok || return false
    all(x -> isfinite(x) && x > zs, (s.electrons_per_carboxylation,
        s.electrons_per_photorespiratory_co2, s.atp_per_carboxylation,
        s.atp_per_photorespiratory_co2)) || return false
    cef_c1, cef_c2 = cytochrome_cyclic_coefficients(
        r.linear_atp_per_electron, r.cyclic_atp_per_electron, s)
    c1 = s.electrons_per_carboxylation + cef_c1
    c2 = s.electrons_per_photorespiratory_co2 + cef_c2
    z = zero(c1)
    return isfinite(c1) && isfinite(c2) && cef_c1 >= z && cef_c2 >= z && c1 > z && c2 > z
end

@inline valid_zscheme_supply(psi_supply, q2, eta) =
    isfinite(psi_supply) && psi_supply >= zero(psi_supply) &&
    isfinite(q2) && q2 >= zero(q2) &&
    isfinite(eta) && psi_supply <= eta * q2

@inline function photosynthesize(p::CytochromePhotosynthesis, env::LeafEnvironment{FT},
    excitation::AbstractExcitation = BroadbandExcitation()) where {FT}
    n = FT(NaN)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(env.co2) && env.co2 >= zero(FT) && isfinite(env.oxygen) && env.oxygen >= zero(FT))
        return (; invalid_carbon(FT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    rates = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    if !valid_cytochrome_rates_and_stoichiometry(rates, p.stoichiometry; require_positive_gamma = false) ||
       (iszero(env.co2) && iszero(rates.gamma_star))
        return (; invalid_carbon(FT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par,
        rates.vqmax, excitation)
    eta = psi_to_psii_electron_ratio(rates, env.co2, p.stoichiometry)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; invalid_carbon(FT)..., electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n, psi_to_psii_ratio = n)
    end
    carbon = carbon_assimilation(rates, env.co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation)
    return (; carbon...,
        electron_transport = supply.psi_electron_supply / eta,
        psi_electron_supply = supply.psi_electron_supply,
        psii_open_fraction = supply.psii_open_fraction,
        psi_to_psii_ratio = eta)
end

@inline function biochemical_turnover(p::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation},
    env::LeafEnvironment{FT}, excitation::AbstractExcitation = BroadbandExcitation()) where {T,E,S,L,FT}
    n = FT(NaN)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(env.co2) && env.co2 >= zero(FT) && isfinite(env.oxygen) && env.oxygen >= zero(FT))
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    if !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false) ||
       (iszero(env.co2) && iszero(r.gamma_star))
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par,
        r.vqmax, excitation)
    eta = psi_to_psii_electron_ratio(r, env.co2, p.stoichiometry)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; carboxylation = n, photorespiratory_co2 = n, respiration = n,
            linear_electron_flow = n, psi_electron_flow = n, cyclic_electron_flow = n,
            electron_supply = n, psi_electron_supply = n, psii_open_fraction = n)
    end
    turnover = biochemical_turnover(r, env.co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation)
    return (; turnover...,
        electron_supply = supply.psi_electron_supply / eta,
        psi_electron_supply = supply.psi_electron_supply,
        psii_open_fraction = supply.psii_open_fraction)
end

@inline function oxygen_exchange(p::CoreOxygenExchange,
    model::CytochromePhotosynthesis{T,E,S,L,SignedAssimilation}, env::LeafEnvironment,
    excitation::AbstractExcitation = BroadbandExcitation()) where {T,E,S,L}
    # Gross O₂ evolution comes strictly from PSII water oxidation (realized linear_electron_flow);
    # cyclic PSI electron flow evolves zero O₂.
    turnover = biochemical_turnover(model, env, excitation)
    oxygenation = turnover.photorespiratory_co2 / p.co2_per_oxygenation
    inputs = OxygenFluxInputs(turnover.linear_electron_flow, oxygenation, turnover.respiration)
    return (; oxygen_exchange(p, inputs)...,
        carboxylation = turnover.carboxylation, rubisco_oxygenation = oxygenation,
        linear_electron_flow = turnover.linear_electron_flow,
        psi_electron_flow = turnover.psi_electron_flow,
        cyclic_electron_flow = turnover.cyclic_electron_flow)
end

@inline function gas_exchange(p::CytochromePhotosynthesis{T,E,S,MinimumLimitation,SignedAssimilation},
    env::GasExchangeEnvironment{FT}, state::StomatalState, d::SeriesDiffusion,
    excitation::AbstractExcitation = BroadbandExcitation()) where {T,E,S,FT}
    n = FT(NaN)
    if !valid_excitation(excitation, env.absorbed_par) ||
       !(isfinite(env.temperature) && env.temperature > zero(FT) &&
         isfinite(env.co2) && env.co2 >= zero(FT) && isfinite(env.oxygen) && env.oxygen >= zero(FT) &&
         isfinite(env.pressure) && env.pressure > zero(FT) &&
         isfinite(state.conductance) && state.conductance > zero(FT))
        return (; carbon = invalid_carbon(FT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    conductance = inv(inv(d.boundary_co2) + d.diffusivity_ratio / state.conductance + inv(d.mesophyll_co2))
    g = conductance * d.mole_to_micromole / env.pressure
    if !(isfinite(g) && g > zero(g))
        return (; carbon = invalid_carbon(FT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    if !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false)
        return (; carbon = invalid_carbon(FT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    c1, c2 = cytochrome_carbon_coefficients(r, p.stoichiometry)
    km = r.kc * (one(env.oxygen) + env.oxygen / r.ko)
    rubisco_co2 = diffusion_co2_root(r.vcmax, km, r.gamma_star, r.respiration, env.co2, g)
    cytochrome_co2 = diffusion_co2_root(supply.psi_electron_supply / c1,
        (c2 / c1) * r.gamma_star, r.gamma_star, r.respiration, env.co2, g)
    above_compensation = g * (env.co2 - r.gamma_star) + r.respiration >= zero(g)
    co2 = above_compensation ? max(rubisco_co2, cytochrome_co2) : min(rubisco_co2, cytochrome_co2)
    eta = psi_to_psii_electron_ratio(r, co2, p.stoichiometry)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    if !isfinite(co2) || (iszero(co2) && iszero(r.gamma_star)) ||
       !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; carbon = invalid_carbon(FT), electron_transport = n,
            psi_electron_supply = n, psii_open_fraction = n,
            surface_co2 = n, intercellular_co2 = n, biochemical_co2 = n,
            conductance = n, diffusion_residual = n)
    end
    carbon = carbon_assimilation(r, co2, env.oxygen, supply.psi_electron_supply,
        p.stoichiometry, p.limitation, p.compensation)
    pressure_flux = carbon.net_assimilation * env.pressure / d.mole_to_micromole
    surface_co2 = env.co2 - pressure_flux / d.boundary_co2
    intercellular_co2 = co2 + pressure_flux / d.mesophyll_co2
    return (; carbon,
        electron_transport = supply.psi_electron_supply / eta,
        psi_electron_supply = supply.psi_electron_supply,
        psii_open_fraction = supply.psii_open_fraction,
        surface_co2, intercellular_co2, biochemical_co2 = co2, conductance,
        diffusion_residual = g * (env.co2 - co2) - carbon.net_assimilation)
end

@inline function cytochrome_equilibrium_result(co2::FT, p, r, psi_supply, q2, env, d, g0, gain,
    width, iterations, status) where {FT}
    eta = psi_to_psii_electron_ratio(r, co2, p.stoichiometry)
    valid_zscheme_supply(psi_supply, q2, eta) ||
        return invalid_equilibrium(FT, EquilibriumInvalidInput)
    e = equilibrium_residual(co2, p, r, psi_supply, env, d, g0, gain)
    carbon = carbon_assimilation(r, co2, env.oxygen, psi_supply,
        p.stoichiometry, p.limitation, p.compensation)
    conductance = inv(e.resistance)
    flux_residual = e.residual * conductance * d.mole_to_micromole / env.pressure
    intercellular = co2 + carbon.net_assimilation * env.pressure / (d.mole_to_micromole * d.mesophyll_co2)
    return EquilibriumSolution(carbon, psi_supply / eta, co2, intercellular, e.surface, e.gs,
        conductance, e.residual, flux_residual, width, iterations, status, status == EquilibriumConverged)
end

@inline function equilibrium_gas_exchange(p::CytochromePhotosynthesis, stomata::SteadyStateStomata,
    env::GasExchangeEnvironment{FT}, d::SeriesDiffusion, solver::EquilibriumSolver,
    excitation::AbstractExcitation = BroadbandExcitation();
    initial_co2 = nothing) where {FT}
    valid_equilibrium_inputs(env, d, stomata.target) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    valid_excitation(excitation, env.absorbed_par) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    g0, gain = equilibrium_coefficients(stomata.target, env)
    isfinite(gain) || return invalid_equilibrium(FT, EquilibriumInvalidInput)
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry) ||
        return invalid_equilibrium(FT, EquilibriumInvalidInput)
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    psi_supply = supply.psi_electron_supply
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    isfinite(psi_supply) && psi_supply >= zero(FT) ||
        return invalid_equilibrium(FT, EquilibriumInvalidInput)
    resistance0 = inv(d.boundary_co2) + d.diffusivity_ratio / g0 + inv(d.mesophyll_co2)
    if iszero(psi_supply) || iszero(r.vcmax)
        x = env.co2 + r.respiration * env.pressure / d.mole_to_micromole * resistance0
        e = equilibrium_residual(x, p, r, psi_supply, env, d, g0, gain)
        status = equilibrium_resolved(e, x, env, d, solver) ? EquilibriumConverged : EquilibriumPrecisionLimit
        return cytochrome_equilibrium_result(x, p, r, psi_supply, q2, env, d, g0, gain, zero(FT), 0, status)
    end
    lo = zero(FT)
    elo = equilibrium_residual(lo, p, r, psi_supply, env, d, g0, gain)
    hi = env.co2 + max(-elo.assimilation, zero(FT)) * env.pressure / d.mole_to_micromole * resistance0
    ehi = equilibrium_residual(hi, p, r, psi_supply, env, d, g0, gain)
    for (x, e) in ((lo, elo), (hi, ehi))
        equilibrium_resolved(e, x, env, d, solver) &&
            return cytochrome_equilibrium_result(x, p, r, psi_supply, q2, env, d, g0, gain,
                hi - lo, 0, EquilibriumConverged)
    end
    isfinite(hi) && isfinite(elo.residual) && isfinite(ehi.residual) &&
        elo.residual >= zero(FT) && ehi.residual <= zero(FT) ||
        return invalid_equilibrium(FT, EquilibriumInvalidBracket)
    x = isnothing(initial_co2) ? midpoint(lo, hi) :
        (isfinite(initial_co2) && lo < initial_co2 < hi ? FT(initial_co2) : midpoint(lo, hi))
    best, best_error = lo, abs(elo.residual)
    if ehi.physical && abs(ehi.residual) < best_error
        best, best_error = hi, abs(ehi.residual)
    end
    for iteration in 1:solver.max_iterations
        e = equilibrium_residual(x, p, r, psi_supply, env, d, g0, gain)
        isfinite(e.residual) || return invalid_equilibrium(FT, EquilibriumInvalidBracket)
        if e.physical && abs(e.residual) < best_error
            best, best_error = x, abs(e.residual)
        end
        equilibrium_resolved(e, x, env, d, solver) &&
            return cytochrome_equilibrium_result(x, p, r, psi_supply, q2, env, d, g0, gain,
                hi - lo, iteration, EquilibriumConverged)
        if e.residual > zero(FT)
            lo = x
        else
            hi = x
        end
        candidate = equilibrium_step(solver.method, x, e, lo, hi, iteration)
        if !(lo < candidate < hi)
            return cytochrome_equilibrium_result(best, p, r, psi_supply, q2, env, d, g0, gain,
                hi - lo, iteration, EquilibriumPrecisionLimit)
        end
        x = candidate
    end
    return cytochrome_equilibrium_result(best, p, r, psi_supply, q2, env, d, g0, gain,
        hi - lo, solver.max_iterations, EquilibriumIterationLimit)
end

@inline function assimilation_sensitivity(p::CytochromePhotosynthesis, env::GasExchangeEnvironment,
    state::StomatalState, d::SeriesDiffusion, exchange = gas_exchange(p, env, state, d),
    excitation::AbstractExcitation = BroadbandExcitation())
    n = oftype(env.absorbed_par, NaN)
    if !valid_excitation(excitation, env.absorbed_par) || !isfinite(exchange.biochemical_co2)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    r = biochemical_rates(p.temperature, env.temperature, env.oxygen)
    if !valid_cytochrome_rates_and_stoichiometry(r, p.stoichiometry; require_positive_gamma = false)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    supply = cytochrome_electron_supply(p.electron_transport, env.absorbed_par, r.vqmax, excitation)
    eta = psi_to_psii_electron_ratio(r, exchange.biochemical_co2, p.stoichiometry)
    q2 = psii_excitation(p.electron_transport, env.absorbed_par, excitation)
    if !valid_zscheme_supply(supply.psi_electron_supply, q2, eta)
        return (; dco2_dgs = n, dnet_assimilation_dgs = n, dassimilation_dgs = n)
    end
    slope = carbon_assimilation_derivative(r, exchange.biochemical_co2, env.oxygen,
        supply.psi_electron_supply, p.stoichiometry, p.limitation, p.compensation)
    conductance = exchange.conductance
    conversion = d.mole_to_micromole / env.pressure
    supply_slope = conductance * conversion
    dg_dgs = (conductance / state.conductance)^2 * d.diffusivity_ratio * conversion
    dco2_dgs = dg_dgs * (env.co2 - exchange.biochemical_co2) / (slope + supply_slope)
    derivative = slope * dco2_dgs
    return (; dco2_dgs, dnet_assimilation_dgs = derivative, dassimilation_dgs = derivative)
end

@inline assimilation_sensitivity(p::CytochromePhotosynthesis, env::GasExchangeEnvironment,
    state::StomatalState, d::SeriesDiffusion, excitation::AbstractExcitation) =
    assimilation_sensitivity(p, env, state, d,
        gas_exchange(p, env, state, d, excitation), excitation)
