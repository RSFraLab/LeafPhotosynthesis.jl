# Minimal standalone dual number to verify temperature AD propagation through ConstantResponse.
struct CytoTestDual{T<:AbstractFloat} <: Real
    v::T
    d::T
end
CytoTestDual(x::T) where {T<:AbstractFloat} = CytoTestDual{T}(x, zero(T))
(::Type{CytoTestDual{T}})(x::Real) where {T<:AbstractFloat} = CytoTestDual{T}(T(x), zero(T))
(::Type{CytoTestDual{T}})(x::CytoTestDual{T}) where {T<:AbstractFloat} = x
Base.promote_rule(::Type{CytoTestDual{T}}, ::Type{S}) where {T<:AbstractFloat,S<:Real} =
    CytoTestDual{promote_type(T, S)}
Base.zero(::Type{CytoTestDual{T}}) where {T<:AbstractFloat} = CytoTestDual(zero(T), zero(T))
Base.zero(x::CytoTestDual) = zero(typeof(x))
Base.one(::Type{CytoTestDual{T}}) where {T<:AbstractFloat} = CytoTestDual(one(T), zero(T))
Base.one(x::CytoTestDual) = one(typeof(x))
Base.isnan(x::CytoTestDual) = isnan(x.v) || isnan(x.d)
Base.isfinite(x::CytoTestDual) = isfinite(x.v) && isfinite(x.d)
Base.:+(a::CytoTestDual{T}, b::CytoTestDual{T}) where {T} = CytoTestDual(a.v + b.v, a.d + b.d)
Base.:-(a::CytoTestDual{T}, b::CytoTestDual{T}) where {T} = CytoTestDual(a.v - b.v, a.d - b.d)
Base.:-(a::CytoTestDual) = CytoTestDual(-a.v, -a.d)
Base.:*(a::CytoTestDual{T}, b::CytoTestDual{T}) where {T} = CytoTestDual(a.v * b.v, a.d * b.v + a.v * b.d)
Base.:/(a::CytoTestDual{T}, b::CytoTestDual{T}) where {T} = CytoTestDual(a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v))
Base.inv(a::CytoTestDual) = CytoTestDual(inv(a.v), -a.d / (a.v * a.v))
Base.abs(a::CytoTestDual) = a.v >= zero(a.v) ? a : -a
Base.exp(a::CytoTestDual) = (ev = exp(a.v); CytoTestDual(ev, ev * a.d))
Base.log1p(a::CytoTestDual) = CytoTestDual(log1p(a.v), a.d / (one(a.v) + a.v))
Base.:<(a::CytoTestDual, b::CytoTestDual) = a.v < b.v
Base.:<=(a::CytoTestDual, b::CytoTestDual) = a.v <= b.v

# Independent 256-bit BigFloat oracle evaluating Johnson & Berry (2021) directly
# from Eqs. 7–30 (solving for S_1 and summing the connected-units geometric series),
# without sharing the production c1/c2 or normalized kn_hat shortcuts.
function cytochrome_oracle(photo::CytochromePhotosynthesis, fl::CytochromeFluorescence,
    site::LeafEnvironment, excitation::AbstractExcitation = BroadbandExcitation())
    T = BigFloat(site.temperature)
    par = BigFloat(site.absorbed_par)
    co2 = BigFloat(site.co2)
    o2 = BigFloat(site.oxygen)

    eval_tr(p::ConstantResponse) = BigFloat(p.value)
    eval_tr(p::Arrhenius) = BigFloat(p.value_ref) * exp(BigFloat(p.activation_energy) /
        BigFloat(p.gas_constant) * (inv(BigFloat(p.temperature_ref)) - inv(T)))
    function eval_tr(p::PeakedArrhenius)
        a = p.arrhenius
        r = BigFloat(a.gas_constant)
        tref = BigFloat(a.temperature_ref)
        ir = (BigFloat(p.entropy) - BigFloat(p.deactivation_energy) / tref) / r
        it = (BigFloat(p.entropy) - BigFloat(p.deactivation_energy) / T) / r
        sp(z) = max(z, big(0)) + log1p(exp(-abs(z)))
        return eval_tr(a) * exp(sp(ir) - sp(it))
    end
    eval_tr(p::CytochromeCapacity) = BigFloat(p.site_density) * eval_tr(p.turnover_rate)
    eval_tr(p::PrescribedCompensationPoint) = eval_tr(p.response)

    tr = photo.temperature
    vcmax = eval_tr(tr.vcmax)
    vqmax = eval_tr(tr.vqmax)
    rd = eval_tr(tr.respiration)
    kc = eval_tr(tr.kc)
    ko = eval_tr(tr.ko)
    gamma = tr.gamma_star isa OxygenScaledCompensationPoint ?
        eval_tr(tr.gamma_star.response) * (o2 / BigFloat(tr.gamma_star.reference_oxygen)) :
        eval_tr(tr.gamma_star)
    nl = eval_tr(tr.linear_atp_per_electron)
    nc = eval_tr(tr.cyclic_atp_per_electron)

    et = photo.electron_transport
    q2_exc = excitation isa BroadbandExcitation ? BigFloat(et.psii_fraction) * par :
        BigFloat(excitation.psii_photons)
    q1_exc = excitation isa BroadbandExcitation ? (big(1) - BigFloat(et.psii_fraction)) * par :
        BigFloat(excitation.psi_photons)

    phi1_max_et = BigFloat(et.psi_maximum_yield)
    j1_light = phi1_max_et * q1_exc
    j1_avail = iszero(vqmax + j1_light) ? big(0) : (vqmax * j1_light) / (vqmax + j1_light)
    q2_open = iszero(vqmax) ? (iszero(j1_light) && iszero(q2_exc) ? big(1) : big(0)) :
        big(1) - j1_avail / vqmax

    s = photo.stoichiometry
    ec = BigFloat(s.electrons_per_carboxylation)
    eg = BigFloat(s.electrons_per_photorespiratory_co2)
    ac_st = BigFloat(s.atp_per_carboxylation)
    ag_st = BigFloat(s.atp_per_photorespiratory_co2)

    lef_cost = ec * co2 + eg * gamma
    eta = iszero(lef_cost) ? (big(1) - nl / nc + (ac_st / ec) / nc) :
        (big(1) - nl / nc + (ac_st * co2 + ag_st * gamma) / (lef_cost * nc))
    j2_avail = j1_avail / eta

    cap_c = vcmax / (co2 + kc * (big(1) + o2 / ko))
    cap_j = iszero(lef_cost) ? j1_avail / ((big(1) - nl / nc) * ec + ac_st / nc) / big(0) :
        j2_avail / lef_cost
    if iszero(lef_cost)
        c1_bf = (big(1) - nl / nc) * ec + ac_st / nc
        c2_bf = (big(1) - nl / nc) * eg + ag_st / nc
        cap_j = j1_avail / (c1_bf * co2 + c2_bf * gamma)
    end
    cap = min(cap_c, cap_j)
    excess = co2 - gamma
    anet = cap * excess - rd
    j2_real = cap * lef_cost
    j1_real = eta * j2_real
    j_cef = j1_real - j2_real

    # Direct Johnson & Berry (2021) Eqs. 17–30 fluorescence solution via S_1 and phi_U2:
    kf2 = BigFloat(fl.psii.fluorescence_rate)
    kd2 = BigFloat(fl.psii.constitutive_heat_rate)
    kp2 = BigFloat(fl.psii.maximum_photochemical_rate)
    ku2 = BigFloat(fl.psii.exciton_transfer_rate)
    kns = BigFloat(fl.psii.sustained_quenching_rate)
    kf1 = BigFloat(fl.psi.fluorescence_rate)
    kd1 = BigFloat(fl.psi.constitutive_heat_rate)
    kp1 = BigFloat(fl.psi.maximum_photochemical_rate)
    kx1 = BigFloat(fl.psi.closed_quenching_rate)

    if iszero(q2_exc) && iszero(q1_exc)
        return (; anet, j1_avail, j2_avail, q2_open, eta, j2_real, j1_real, j_cef, dark = true)
    end

    phi_p1 = j1_real / q1_exc
    s_open1 = kf1 + kd1 + kp1
    s_closed1 = kf1 + kd1 + kx1
    phi1_max = kp1 / s_open1
    q1_open = min(big(1), phi_p1 / phi1_max)
    phi_f1 = q1_open * (kf1 / s_open1) + (big(1) - q1_open) * (kf1 / s_closed1)
    phi_d1 = q1_open * (kd1 / s_open1) + (big(1) - q1_open) * (kd1 / s_closed1)
    phi_n1 = (big(1) - q1_open) * (kx1 / s_closed1)

    phi_p2 = j2_real / q2_exc
    a_q = phi_p2
    b_q = -(ku2 * phi_p2 + kp2 * (q2_open - phi_p2))
    c_q = -(phi_p2 * (big(1) - q2_open) * ku2 * kp2)
    s1_raw = (-b_q + sqrt(max(big(0), b_q^2 - big(4) * a_q * c_q))) / (big(2) * a_q)
    s_base = kf2 + kd2 + ku2 + kns
    kn = max(big(0), s1_raw - s_base)
    s1 = s_base + kn
    s2 = s1 + kp2
    phi_u2 = q2_open * (ku2 / s2) + (big(1) - q2_open) * (ku2 / s1)
    inv_one_minus_u2 = inv(big(1) - phi_u2)
    phi_f2 = (q2_open * (kf2 / s2) + (big(1) - q2_open) * (kf2 / s1)) * inv_one_minus_u2
    phi_d2 = (q2_open * (kd2 / s2) + (big(1) - q2_open) * (kd2 / s1)) * inv_one_minus_u2
    phi_n2 = (q2_open * ((kns + kn) / s2) + (big(1) - q2_open) * ((kns + kn) / s1)) * inv_one_minus_u2
    npq = (kn + kns) / (kf2 + kd2)
    f2_flux = q2_exc * phi_f2
    f1_flux = q1_exc * phi_f1
    f_tot = f2_flux + f1_flux
    phi_f_tot = f_tot / (q1_exc + q2_exc)

    return (; anet, j1_avail, j2_avail, q2_open, eta, j2_real, j1_real, j_cef, dark = false,
        phi_p2, phi_f2, phi_d2, phi_n2, phi_p1, phi_f1, phi_d1, phi_n1,
        phi_f_tot, q1_open, kn, npq, f2_flux, f1_flux, f_tot)
end

function cytochrome_allocations(photo, fl, site, gas_env, state, diff, solver, leaf)
    photosynthesize(photo, site)
    biochemical_turnover(photo, site)
    inputs = photosystem_flux_inputs(photo, site)
    fluorescence(fl, inputs)
    gas_exchange(photo, gas_env, state, diff)
    equilibrium_gas_exchange(leaf, gas_env, diff, solver)
    return @allocated(photosynthesize(photo, site)) +
           @allocated(biochemical_turnover(photo, site)) +
           @allocated(photosystem_flux_inputs(photo, site)) +
           @allocated(fluorescence(fl, inputs)) +
           @allocated(gas_exchange(photo, gas_env, state, diff)) +
           @allocated(equilibrium_gas_exchange(leaf, gas_env, diff, solver))
end

@testset "Cytochrome b6f photosynthesis and coupled PSII/PSI fluorescence" begin
    for FT in (Float32, Float64)
        case = setup_cytochrome_stress(FT)
        photo = case.cyto_leaf.photosynthesis
        x, v = case.settings["experiment"], case.settings["validation"]
        key = FT === Float32 ? "float32" : "float64"
        rtol, atol = FT(v[key * "_relative_tolerance"]), FT(v[key * "_flux_absolute_tolerance"])
        T, O2 = x["temperature_K"], x["oxygen_Pa"]

        @test isbitstype(typeof(photo))
        @test isbitstype(typeof(case.jb_connected))
        @test isbitstype(typeof(case.cyto_leaf))

        # 1. ConstantResponse & SerialLimitation primitives.
        cr = ConstantResponse(FT(0.75))
        @test (@inferred temperature_response(cr, T)) === FT(0.75)
        @test isnan(temperature_response(cr, FT(NaN)))
        @test_throws ArgumentError ConstantResponse(FT(NaN))
        @test_throws ArgumentError ConstantResponse(FT(Inf))
        @test_throws MethodError ConstantResponse(1)
        @test_throws MethodError CytochromeCapacity(1, photo.temperature.vqmax.turnover_rate)
        @test_throws ArgumentError CytochromeCapacity(FT(-0.1), photo.temperature.vqmax.turnover_rate)
        @test_throws ArgumentError CytochromeElectronParameters(;
            psii_fraction = zero(FT),
            psi_maximum_yield = photo.electron_transport.psi_maximum_yield,
            limitation = SerialLimitation(),
        )
        @test_throws ArgumentError CytochromeElectronParameters(;
            psii_fraction = one(FT),
            psi_maximum_yield = photo.electron_transport.psi_maximum_yield,
            limitation = SerialLimitation(),
        )

        # Temperature AD propagation through ConstantResponse + CytochromeRates.
        t_dual = CytoTestDual(T, one(FT))
        r_dual = @inferred biochemical_rates(photo.temperature, t_dual, CytoTestDual(O2))
        @test r_dual isa CytochromeRates{CytoTestDual{FT}}
        @test iszero(r_dual.linear_atp_per_electron.d) && iszero(r_dual.cyclic_atp_per_electron.d)
        @test isnan(temperature_response(cr, CytoTestDual(FT(NaN), zero(FT))).v)

        # Float32 exact cyclic deficit (0 ULP cancellation when a_c == eta_L * e_c at non-power-of-two eta_C).
        if FT === Float32
            r_077 = CytochromeRates(FT(60), FT(195), FT(0.9), FT(40), FT(27000), FT(0), FT(0.75), FT(0.77))
            c1_077, c2_077 = cytochrome_carbon_coefficients(r_077, photo.stoichiometry)
            @test c1_077 === FT(4)
            @test c2_077 > FT(8)
            @test validate_fluorescence_coupling(case.jb_connected, photo.electron_transport,
                photo.stoichiometry; linear_atp_per_electron = FT(0.75), cyclic_atp_per_electron = FT(0.77)) === nothing
            turn_077 = biochemical_turnover(r_077, FT(1000), O2, FT(120), photo.stoichiometry,
                MinimumLimitation(), SignedAssimilation())
            @test turn_077.cyclic_electron_flow === zero(FT)
        end

        sl = SerialLimitation()
        @test (@inferred colimit(zero(FT), zero(FT), sl)) === zero(FT)
        @test colimit(FT(3), FT(6), sl) ≈ FT(2)
        @test colimit(FT(1e-30), FT(1e30), sl) ≈ FT(1e-30)
        @test colimit(FT(5), FT(Inf), sl) === FT(5)
        @test colimit(FT(Inf), FT(5), sl) === FT(5)
        @test isnan(colimit(-one(FT), one(FT), sl))

        # 2. Setup coupling validation: SerialLimitation, negative-CEF rejection, and beta_min check.
        beta_min = minimum_psii_partition(case.jb_connected, photo.stoichiometry;
            linear_atp_per_electron = FT(0.75), cyclic_atp_per_electron = one(FT))
        @test isapprox(beta_min, FT(0.5211416490486258); rtol)
        @test validate_fluorescence_coupling(case.jb_connected, photo, T, O2) === nothing

        for non_serial_lim in (MinimumLimitation(), QuadraticLimitation(FT(0.9)))
            non_serial_et = CytochromeElectronParameters(;
                psii_fraction = photo.electron_transport.psii_fraction,
                psi_maximum_yield = photo.electron_transport.psi_maximum_yield,
                limitation = non_serial_lim,
            )
            @test_throws ArgumentError validate_fluorescence_coupling(
                case.jb_connected, non_serial_et, photo.stoichiometry;
                linear_atp_per_electron = FT(0.75), cyclic_atp_per_electron = one(FT))
        end

        bad_beta_et = CytochromeElectronParameters(;
            psii_fraction = FT(0.50),
            psi_maximum_yield = photo.electron_transport.psi_maximum_yield,
            limitation = SerialLimitation(),
        )
        @test_throws ArgumentError validate_fluorescence_coupling(
            case.jb_connected, bad_beta_et, photo.stoichiometry;
            linear_atp_per_electron = FT(0.75), cyclic_atp_per_electron = one(FT))
        # Even if constructed without setup validation, runtime rejects super-maximal PSII yield (never negative KN).
        bad_photo = CytochromePhotosynthesis(;
            photo.temperature, electron_transport = bad_beta_et,
            photo.stoichiometry, photo.limitation, photo.compensation)
        low_light_site = LeafEnvironment(T, FT(5), FT(35), O2)
        bad_res = fluorescence(case.jb_connected, photosystem_flux_inputs(bad_photo, low_light_site))
        @test bad_res.status === FluorescenceIncompatibleRates
        @test isnan(bad_res.psii_photon_flux) && isnan(bad_res.psii_quenching_rate)

        # Negative cyclic electron flow (eta < 1, e.g. linear_atp_per_electron = 0.85) rejected everywhere.
        bad_cef_tr = CytochromeTemperatureResponse(;
            vcmax = photo.temperature.vcmax,
            vqmax = photo.temperature.vqmax,
            respiration = photo.temperature.respiration,
            kc = photo.temperature.kc,
            ko = photo.temperature.ko,
            gamma_star = photo.temperature.gamma_star,
            linear_atp_per_electron = ConstantResponse(FT(0.85)),
            cyclic_atp_per_electron = ConstantResponse(one(FT)),
        )
        bad_cef_photo = CytochromePhotosynthesis(;
            temperature = bad_cef_tr, electron_transport = photo.electron_transport,
            stoichiometry = photo.stoichiometry, limitation = photo.limitation,
            compensation = photo.compensation)
        probe_site = LeafEnvironment(T, FT(800), FT(28), O2)
        probe_gas = GasExchangeEnvironment(T, FT(800), FT(41), O2, FT(101325), FT(0.55), FT(1500))
        probe_state = StomatalState(FT(0.18))
        @test_throws ArgumentError minimum_psii_partition(case.jb_connected, photo.stoichiometry;
            linear_atp_per_electron = FT(0.85), cyclic_atp_per_electron = one(FT))
        @test isnan((@inferred photosynthesize(bad_cef_photo, probe_site)).net_assimilation)
        @test isnan((@inferred biochemical_turnover(bad_cef_photo, probe_site)).cyclic_electron_flow)
        @test isnan((@inferred oxygen_exchange(case.cyto_leaf.oxygen, bad_cef_photo, probe_site)).gross_production)
        @test isnan((@inferred gas_exchange(bad_cef_photo, probe_gas, probe_state, case.diffusion)).biochemical_co2)
        @test isnan((@inferred assimilation_sensitivity(bad_cef_photo, probe_gas, probe_state, case.diffusion)).dco2_dgs)
        @test (@inferred equilibrium_gas_exchange(bad_cef_photo, case.cyto_leaf.stomata,
            probe_gas, case.diffusion, case.solver)).status === EquilibriumInvalidInput

        # 3. Oracle check across light, sub-compensation, exact compensation, connected & puddle, sustained KN, and K_X1 != K_P1.
        sustained_fl = CytochromeFluorescence(;
            psii = PSIIRateParameters(;
                fluorescence_rate = FT(0.05), constitutive_heat_rate = FT(0.55),
                maximum_photochemical_rate = FT(4.5), exciton_transfer_rate = FT(1.5),
                sustained_quenching_rate = FT(0.4),
            ),
            psi = PSIRateParameters(;
                fluorescence_rate = FT(0.05), constitutive_heat_rate = FT(0.55),
                maximum_photochemical_rate = FT(14.5), closed_quenching_rate = FT(10.0),
            ),
            yield_tolerance = case.jb_connected.yield_tolerance,
        )
        rates = biochemical_rates(photo.temperature, T, O2)
        for light in FT.(x["test_light"]), factor in FT.(x["test_co2_gamma_factors"])
            site = LeafEnvironment(T, light, rates.gamma_star * factor, O2)
            c_out = @inferred photosynthesize(photo, site)
            t_out = @inferred biochemical_turnover(photo, site)
            o_out = @inferred oxygen_exchange(case.cyto_leaf, site)
            inputs = @inferred photosystem_flux_inputs(photo, site)

            # Gross O2 comes strictly from water-splitting linear electron flow, not cyclic PSI flow.
            @test isapprox(o_out.gross_production, t_out.linear_electron_flow / FT(4); rtol, atol)
            @test isapprox(t_out.linear_electron_flow + t_out.cyclic_electron_flow, t_out.psi_electron_flow; rtol, atol)

            for fl_model in (case.jb_connected, case.jb_puddle, sustained_fl)
                # Note: sustained_fl has lower phi2_max(1), so low light at high Cc can exceed its quenched maximum.
                oracle = cytochrome_oracle(photo, fl_model, site)
                res = @inferred fluorescence(fl_model, inputs)
                @test isbitstype(typeof(res))
                if oracle.dark
                    @test res.status === FluorescenceDark
                    @test res.psii_photon_flux === zero(FT)
                    @test res.psi_photon_flux === zero(FT)
                    @test res.total_photon_flux === zero(FT)
                    @test all(isnan, (res.psii_fluorescence_yield, res.psi_fluorescence_yield,
                        res.total_fluorescence_yield, res.psii_quenching_rate, res.npq))
                elseif fl_model === sustained_fl && !fluorescence_valid(res)
                    @test res.status === FluorescenceIncompatibleRates
                else
                    @test res.status === FluorescenceLit
                    @test isapprox(c_out.net_assimilation, oracle.anet; rtol, atol)
                    @test isapprox(c_out.psi_electron_supply, oracle.j1_avail; rtol, atol)
                    @test isapprox(c_out.electron_transport, oracle.j2_avail; rtol, atol)
                    @test isapprox(c_out.psii_open_fraction, oracle.q2_open; rtol, atol)
                    @test isapprox(c_out.psi_to_psii_ratio, oracle.eta; rtol, atol)
                    @test isapprox(t_out.linear_electron_flow, oracle.j2_real; rtol, atol)
                    @test isapprox(t_out.psi_electron_flow, oracle.j1_real; rtol, atol)
                    @test isapprox(t_out.cyclic_electron_flow, oracle.j_cef; rtol, atol)

                    @test isapprox(res.psii_photochemical_yield, oracle.phi_p2; rtol, atol)
                    @test isapprox(res.psii_fluorescence_yield, oracle.phi_f2; rtol, atol)
                    @test isapprox(res.psii_constitutive_heat_yield, oracle.phi_d2; rtol, atol)
                    @test isapprox(res.psii_regulated_heat_yield, oracle.phi_n2; rtol, atol)
                    @test isapprox(res.psi_photochemical_yield, oracle.phi_p1; rtol, atol)
                    @test isapprox(res.psi_fluorescence_yield, oracle.phi_f1; rtol, atol)
                    @test isapprox(res.psi_constitutive_heat_yield, oracle.phi_d1; rtol, atol)
                    @test isapprox(res.psi_regulated_heat_yield, oracle.phi_n1; rtol, atol)
                    @test isapprox(res.total_fluorescence_yield, oracle.phi_f_tot; rtol, atol)
                    @test isapprox(res.psi_open_fraction, oracle.q1_open; rtol, atol)
                    @test isapprox(res.psii_quenching_rate, oracle.kn; rtol, atol)
                    @test isapprox(res.npq, oracle.npq; rtol, atol)
                    @test isapprox(res.psii_photon_flux, oracle.f2_flux; rtol, atol)
                    @test isapprox(res.psi_photon_flux, oracle.f1_flux; rtol, atol)
                    @test isapprox(res.total_photon_flux, oracle.f_tot; rtol, atol)

                    psii_yields = (res.psii_photochemical_yield, res.psii_fluorescence_yield,
                        res.psii_constitutive_heat_yield, res.psii_regulated_heat_yield)
                    psi_yields = (res.psi_photochemical_yield, res.psi_fluorescence_yield,
                        res.psi_constitutive_heat_yield, res.psi_regulated_heat_yield)
                    @test all(y -> zero(FT) <= y <= one(FT), psii_yields)
                    @test all(y -> zero(FT) <= y <= one(FT), psi_yields)
                    @test isapprox(sum(psii_yields), one(FT); rtol)
                    @test isapprox(sum(psi_yields), one(FT); rtol)
                    factor == one(FT) && @test t_out.linear_electron_flow > zero(FT)
                end
            end
        end

        # 4. Zero heap allocations across all scalar entry points.
        gas_env = GasExchangeEnvironment(T, FT(800), FT(41), O2, FT(101325), FT(0.55), FT(1500))
        state = StomatalState(FT(0.18))
        lit_site = LeafEnvironment(T, FT(800), FT(28), O2)
        @test cytochrome_allocations(photo, case.jb_connected, lit_site, gas_env,
            state, case.diffusion, case.solver, case.cyto_leaf) == 0

        # 5. Analytical gas_exchange, Bisection vs SafeguardedNewton equilibrium, and BigFloat sensitivity FD check.
        ex = @inferred gas_exchange(case.cyto_leaf, gas_env, state, case.diffusion)
        @test abs(ex.diffusion_residual) <= atol
        eq_newton = @inferred equilibrium_gas_exchange(case.cyto_leaf, gas_env, case.diffusion, case.solver)
        bisect_solver = EquilibriumSolver(;
            method = Bisection(),
            case.solver.co2_atol, case.solver.co2_rtol,
            case.solver.flux_atol, case.solver.flux_rtol,
            max_iterations = 120,
        )
        eq_bisect = @inferred equilibrium_gas_exchange(case.cyto_leaf, gas_env, case.diffusion, bisect_solver)
        @test eq_newton.converged && eq_bisect.converged
        tol_c = case.solver.co2_atol + case.solver.co2_rtol * max(abs(gas_env.co2), abs(eq_newton.biochemical_co2))
        tol_a = case.solver.flux_atol + case.solver.flux_rtol * abs(eq_newton.carbon.net_assimilation)
        @test isapprox(eq_newton.biochemical_co2, eq_bisect.biochemical_co2; atol = 8 * tol_c, rtol = 64 * eps(FT))
        @test isapprox(eq_newton.carbon.net_assimilation, eq_bisect.carbon.net_assimilation; atol = 8 * tol_a, rtol = 64 * eps(FT))

        # Smooth carbon co-limitation (SerialLimitation & QuadraticLimitation) in equilibrium.
        for lim in (SerialLimitation(), QuadraticLimitation(FT(0.95)))
            smooth_photo = CytochromePhotosynthesis(;
                photo.temperature, photo.electron_transport,
                photo.stoichiometry, limitation = lim, photo.compensation)
            s_newt = equilibrium_gas_exchange(smooth_photo, case.cyto_leaf.stomata, gas_env, case.diffusion, case.solver)
            s_bis = equilibrium_gas_exchange(smooth_photo, case.cyto_leaf.stomata, gas_env, case.diffusion, bisect_solver)
            @test s_newt.converged && s_bis.converged
            sc_tol = case.solver.co2_atol + case.solver.co2_rtol * max(abs(gas_env.co2), abs(s_newt.biochemical_co2))
            @test isapprox(s_newt.biochemical_co2, s_bis.biochemical_co2; atol = 8 * sc_tol, rtol = 64 * eps(FT))
        end

        # Analytic assimilation_sensitivity vs 256-bit BigFloat central finite differences (including eta(Cc) derivative).
        sens = @inferred assimilation_sensitivity(photo, gas_env, state, case.diffusion, ex)
        if FT === Float64
            bf_case = setup_cytochrome_stress(BigFloat)
            bf_env = GasExchangeEnvironment(BigFloat(T), BigFloat(800), BigFloat(41),
                BigFloat(O2), BigFloat(101325), BigFloat(0.55), BigFloat(1500))
            gs0 = BigFloat(0.18)
            h = big"1e-18"
            ex_hi = gas_exchange(bf_case.cyto_leaf.photosynthesis, bf_env, StomatalState(gs0 + h), bf_case.diffusion)
            ex_lo = gas_exchange(bf_case.cyto_leaf.photosynthesis, bf_env, StomatalState(gs0 - h), bf_case.diffusion)
            fd_dco2 = (ex_hi.biochemical_co2 - ex_lo.biochemical_co2) / (2 * h)
            fd_da = (ex_hi.carbon.net_assimilation - ex_lo.carbon.net_assimilation) / (2 * h)
            @test isapprox(sens.dco2_dgs, Float64(fd_dco2); rtol = 1e-11)
            @test isapprox(sens.dnet_assimilation_dgs, Float64(fd_da); rtol = 1e-11)
        end

        # Prognostic leaf_tendencies integration.
        prog_leaf = LeafModel(;
            photosynthesis = photo,
            stomata = PrognosticStomata(; target = case.cyto_leaf.stomata.target,
                opening_timescale = FT(600), closing_timescale = FT(300)),
            oxygen = case.cyto_leaf.oxygen,
        )
        tend = @inferred leaf_tendencies(prog_leaf, gas_env, state, case.diffusion)
        @test isfinite(tend.conductance_tendency) && isfinite(tend.oxygen.net_release)

        # PrescribedCompensationPoint & OxygenScaledCompensationPoint oracle checks across O2 levels.
        for cp in (
            PrescribedCompensationPoint(photo.temperature.gamma_star),
            OxygenScaledCompensationPoint(photo.temperature.gamma_star, O2),
        )
            cp_tr = CytochromeTemperatureResponse(;
                vcmax = photo.temperature.vcmax, vqmax = photo.temperature.vqmax,
                respiration = photo.temperature.respiration, kc = photo.temperature.kc,
                ko = photo.temperature.ko, gamma_star = cp,
                linear_atp_per_electron = photo.temperature.linear_atp_per_electron,
                cyclic_atp_per_electron = photo.temperature.cyclic_atp_per_electron,
            )
            cp_photo = CytochromePhotosynthesis(;
                temperature = cp_tr, electron_transport = photo.electron_transport,
                stoichiometry = photo.stoichiometry, limitation = photo.limitation,
                compensation = photo.compensation)
            for o2_val in (zero(FT), O2 / FT(2), O2, FT(35000))
                cp_site = LeafEnvironment(T, FT(600), FT(28), o2_val)
                cp_oracle = cytochrome_oracle(cp_photo, case.jb_connected, cp_site)
                cp_c = @inferred photosynthesize(cp_photo, cp_site)
                cp_fl = @inferred fluorescence(case.jb_connected, photosystem_flux_inputs(cp_photo, cp_site))
                @test isapprox(cp_c.net_assimilation, FT(cp_oracle.anet); rtol, atol)
                @test fluorescence_valid(cp_fl)
                @test isapprox(cp_fl.psii_fluorescence_yield, FT(cp_oracle.phi_f2); rtol)
            end
        end

        # 6. PrescribedDualExcitation, Z-scheme photon conservation, & rate-scaling invariance.
        dual = PrescribedDualExcitation(
            photo.electron_transport.psii_fraction * lit_site.absorbed_par,
            (one(FT) - photo.electron_transport.psii_fraction) * lit_site.absorbed_par,
        )
        @test valid_excitation(dual, lit_site.absorbed_par)
        @test photosynthesize(photo, lit_site, dual) == photosynthesize(photo, lit_site)
        @test fluorescence(case.jb_connected, photosystem_flux_inputs(photo, lit_site, dual)) ==
              fluorescence(case.jb_connected, photosystem_flux_inputs(photo, lit_site))
        # Sub-unity sum (non-photosynthetic pigment absorption) is valid.
        sub_dual = PrescribedDualExcitation(FT(300), FT(280))
        @test valid_excitation(sub_dual, FT(800))
        @test fluorescence_valid(fluorescence(case.jb_connected, photosystem_flux_inputs(photo, lit_site, sub_dual)))
        # Also supported on FvCB ElectronTransportParameters.
        @test isfinite(photosynthesize(case.fvcb_leaf, lit_site, sub_dual).net_assimilation)
        # Super-PAR sum, wrong FT, non-finite, single-photosystem, and PSII-starved (Phi_P2 > 1) rejected across all entry points.
        @test !valid_excitation(PrescribedDualExcitation(FT(500), FT(400)), FT(800))
        other_ft = FT === Float32 ? Float64 : Float32
        @test !valid_excitation(PrescribedDualExcitation(other_ft(100), other_ft(100)), FT(800))
        for bad_exc in (
            PrescribedDualExcitation(FT(500), FT(400)),
            PrescribedDualExcitation(-one(FT), FT(100)),
            PrescribedDualExcitation(FT(NaN), FT(100)),
            PrescribedDualExcitation(FT(Inf), FT(100)),
            PrescribedDualExcitation(other_ft(100), other_ft(100)),
            PrescribedDualExcitation(zero(FT), FT(50)),   # zero PSII photons, lit PSI
            PrescribedDualExcitation(FT(50), zero(FT)),   # lit PSII, zero PSI photons
            PrescribedDualExcitation(one(FT), FT(50)),    # PSII-starved far-red (J_P680,avail > Q2)
        )
            @test isnan((@inferred photosynthesize(photo, lit_site, bad_exc)).net_assimilation)
            @test isnan((@inferred biochemical_turnover(photo, lit_site, bad_exc)).linear_electron_flow)
            @test isnan((@inferred oxygen_exchange(case.cyto_leaf, lit_site, bad_exc)).gross_production)
            @test isnan((@inferred gas_exchange(photo, gas_env, state, case.diffusion, bad_exc)).biochemical_co2)
            @test isnan((@inferred assimilation_sensitivity(photo, gas_env, state, case.diffusion, bad_exc)).dco2_dgs)
            @test (@inferred equilibrium_gas_exchange(photo, case.cyto_leaf.stomata,
                gas_env, case.diffusion, case.solver, bad_exc)).status === EquilibriumInvalidInput
            bad_fl = @inferred fluorescence(case.jb_connected, photosystem_flux_inputs(photo, lit_site, bad_exc))
            @test !fluorescence_valid(bad_fl)
        end
        @test_throws MethodError photosynthesize(photo, lit_site, PrescribedPSIIExcitation(FT(300)))
        @test_throws MethodError photosystem_flux_inputs(photo, gas_env)

        for scale in (FT(7), FT(1e12))
            p2, p1 = case.jb_connected.psii, case.jb_connected.psi
            scaled_fl = CytochromeFluorescence(;
                psii = PSIIRateParameters(p2.fluorescence_rate * scale, p2.constitutive_heat_rate * scale,
                    p2.maximum_photochemical_rate * scale, p2.exciton_transfer_rate * scale,
                    p2.sustained_quenching_rate * scale),
                psi = PSIRateParameters(p1.fluorescence_rate * scale, p1.constitutive_heat_rate * scale,
                    p1.maximum_photochemical_rate * scale, p1.closed_quenching_rate * scale),
                yield_tolerance = case.jb_connected.yield_tolerance,
            )
            inp = photosystem_flux_inputs(photo, lit_site)
            r_base = fluorescence(case.jb_connected, inp)
            r_scaled = fluorescence(scaled_fl, inp)
            for f in (:psii_photon_flux, :psi_photon_flux, :total_photon_flux,
                      :psii_fluorescence_yield, :psi_fluorescence_yield, :total_fluorescence_yield, :npq)
                @test isapprox(getproperty(r_base, f), getproperty(r_scaled, f); rtol)
            end
            @test isapprox(r_base.psii_quenching_rate * scale, r_scaled.psii_quenching_rate; rtol)
        end

        # 7. Boundary roundoff & invalid/incompatible CoupledPhotosystemInputs.
        phi1_max = psi_maximum_photochemical_yield(case.jb_connected)
        phi2_max_half = psii_maximum_photochemical_yield(case.jb_connected, FT(0.5))
        b_inp = CoupledPhotosystemInputs(one(FT), one(FT), nextfloat(phi2_max_half), nextfloat(phi1_max), FT(0.5))
        b_res = fluorescence(case.jb_connected, b_inp)
        @test fluorescence_valid(b_res) && b_res.boundary_roundoff && iszero(b_res.psii_quenching_rate)
        over_inp = CoupledPhotosystemInputs(one(FT), one(FT),
            nextfloat(phi2_max_half + case.jb_connected.yield_tolerance), phi1_max / FT(2), FT(0.5))
        @test fluorescence(case.jb_connected, over_inp).status === FluorescenceIncompatibleRates
        for bad_inp in (
            CoupledPhotosystemInputs(-one(FT), one(FT), zero(FT), zero(FT), FT(0.5)),
            CoupledPhotosystemInputs(zero(FT), zero(FT), one(FT), zero(FT), one(FT)),
            CoupledPhotosystemInputs(one(FT), one(FT), FT(0.2), FT(0.2), FT(1.2)),
            CoupledPhotosystemInputs(FT(NaN), one(FT), FT(0.2), FT(0.2), FT(0.5)),
        )
            @test fluorescence(case.jb_connected, bad_inp).status === FluorescenceInvalidInput
        end
    end
end
