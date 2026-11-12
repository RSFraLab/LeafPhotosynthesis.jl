function additions_allocations(rule, obs, fl, photo, leaf, site, res, a2, a1)
    dark_balanced_psii_partition(fl, photo, site)
    state_transition_partition(rule, fl, photo, site)
    state_transition_excitation(rule, fl, leaf, site)
    dark_pam_levels(obs, fl, a2, a1)
    pam_fluorescence(obs, fl, res, a2, a1)
    n1 = @allocated dark_balanced_psii_partition(fl, photo, site)
    n2 = @allocated state_transition_partition(rule, fl, photo, site)
    n3 = @allocated state_transition_excitation(rule, fl, leaf, site)
    n4 = @allocated dark_pam_levels(obs, fl, a2, a1)
    n5 = @allocated pam_fluorescence(obs, fl, res, a2, a1)
    return n1 + n2 + n3 + n4 + n5
end

# Independent 256-bit BigFloat oracle for Johnson & Berry (2021) Eqs. 38–42 PAM levels.
function pam_bigfloat_oracle(obs::PAMFluorescenceObserver, fl::CytochromeFluorescence,
    res::CytochromeFluorescenceResult, a2, a1, a2_d = a2, a1_d = a1)
    kf2, kd2, kp2, kns = BigFloat(fl.psii.fluorescence_rate), BigFloat(fl.psii.constitutive_heat_rate),
        BigFloat(fl.psii.maximum_photochemical_rate), BigFloat(fl.psii.sustained_quenching_rate)
    kf1, kd1, kp1, kx1 = BigFloat(fl.psi.fluorescence_rate), BigFloat(fl.psi.constitutive_heat_rate),
        BigFloat(fl.psi.maximum_photochemical_rate), BigFloat(fl.psi.closed_quenching_rate)
    e2, e1 = BigFloat(obs.psii_sensitivity), BigFloat(obs.psi_sensitivity)
    kn = BigFloat(res.psii_quenching_rate)
    ba2, ba1, ba2d, ba1d = BigFloat(a2), BigFloat(a1), BigFloat(a2_d), BigFloat(a1_d)

    fm_d = ba2d * (kf2 / (kf2 + kd2 + kns)) * e2 + ba1d * (kf1 / (kf1 + kd1 + kx1)) * e1
    fo_d = ba2d * (kf2 / (kf2 + kd2 + kns + kp2)) * e2 + ba1d * (kf1 / (kf1 + kd1 + kp1)) * e1
    fmp = ba2 * (kf2 / (kf2 + kd2 + kns + kn)) * e2 + ba1 * (kf1 / (kf1 + kd1 + kx1)) * e1
    fop = ba2 * (kf2 / (kf2 + kd2 + kns + kn + kp2)) * e2 + ba1 * (kf1 / (kf1 + kd1 + kp1)) * e1
    fs = ba2 * BigFloat(res.psii_fluorescence_yield) * e2 + ba1 * BigFloat(res.psi_fluorescence_yield) * e1

    return (;
        fs, fm_dark = fm_d, fo_dark = fo_d, fm_prime = fmp, fo_prime = fop,
        apparent_photochemical_yield = big(1) - fs / fmp,
        apparent_regulated_heat_yield = fs / fmp - fs / fm_d,
        apparent_residual_yield = fs / fm_d,
        apparent_npq = fm_d / fmp - big(1),
        qp = (fmp - fs) / (fmp - fop),
        ql = ((fmp - fs) / (fmp - fop)) * (fop / fs),
        fv_fm = (fm_d - fo_d) / fm_d,
    )
end

@testset "State transitions (qT) and PAM observation operator" begin
    @test isempty(Test.detect_ambiguities(LeafPhotosynthesis))
    @test_throws UndefKeywordError StateTransitionRule()
    @test_throws UndefKeywordError PAMFluorescenceObserver()

    for FT in (Float32, Float64)
        @testset "$FT" begin
            case = setup_cytochrome_stress(FT)
            photo = case.cyto_leaf.photosynthesis
            fl_conn = case.jb_connected
            fl_pud = case.jb_puddle
            v = case.settings["validation"]
            key = FT === Float32 ? "float32" : "float64"
            rtol = FT(v[key * "_relative_tolerance"])
            atol = FT(v[key * "_flux_absolute_tolerance"])
            T_K = case.settings["experiment"]["temperature_K"]
            O2_Pa = case.settings["experiment"]["oxygen_Pa"]

            # 1. Constructor & domain guards.
            rule = @inferred StateTransitionRule{FT}()
            @test rule.minimum_psii_fraction === zero(FT) && rule.maximum_psii_fraction === one(FT)
            bounded_rule = @inferred StateTransitionRule(; minimum_psii_fraction = FT(0.35), maximum_psii_fraction = FT(0.65))
            @test bounded_rule == StateTransitionRule(FT(0.35), FT(0.65))
            for (lo, hi) in ((-FT(0.1), one(FT)), (zero(FT), FT(1.1)), (FT(0.5), FT(0.5)),
                             (FT(0.6), FT(0.4)), (FT(NaN), one(FT)), (zero(FT), FT(Inf)))
                @test_throws ArgumentError StateTransitionRule(lo, hi)
            end

            obs_ideal = @inferred PAMFluorescenceObserver(; psii_sensitivity = one(FT), psi_sensitivity = zero(FT))
            obs_pam = @inferred PAMFluorescenceObserver(one(FT), FT(2))
            for (e2, e1) in ((zero(FT), one(FT)), (-one(FT), one(FT)), (one(FT), -FT(0.1)),
                             (FT(NaN), one(FT)), (one(FT), FT(Inf)))
                @test_throws ArgumentError PAMFluorescenceObserver(e2, e1)
            end

            # 2. Exact Johnson & Berry (2021) Table 2 Case 10 benchmarks (0%, ~17%, ~22% excursions) & K_N = 0 below saturation.
            S_jb = (FT(3.6) / FT(260e-6)) * (FT(179e-3) / FT(3.6 * 0.27))
            jb_temp = CytochromeTemperatureResponse(;
                vcmax = ConstantResponse(FT(100)),
                vqmax = CytochromeCapacity(FT(350 / 300), ConstantResponse(FT(300))),
                respiration = ConstantResponse(one(FT)),
                kc = ConstantResponse(FT(26)),
                ko = ConstantResponse(FT(17900)),
                gamma_star = OxygenScaledCompensationPoint(ConstantResponse(FT(0.5 * 1e5) / S_jb), FT(1e5)),
                linear_atp_per_electron = ConstantResponse(FT(0.75)),
                cyclic_atp_per_electron = ConstantResponse(one(FT)),
            )
            jb_photo = CytochromePhotosynthesis(;
                temperature = jb_temp,
                electron_transport = CytochromeElectronParameters(;
                    psii_fraction = FT(0.525),
                    psi_maximum_yield = psi_maximum_photochemical_yield(fl_conn),
                    limitation = SerialLimitation(),
                ),
                stoichiometry = photo.stoichiometry,
                limitation = MinimumLimitation(),
                compensation = SignedAssimilation(),
            )
            fl_ku10 = CytochromeFluorescence(;
                psii = PSIIRateParameters(FT(0.05), FT(0.55), FT(4.5), FT(10), zero(FT)),
                psi = fl_conn.psi,
                yield_tolerance = fl_conn.yield_tolerance,
            )

            dark_site = LeafEnvironment(FT(298.15), zero(FT), FT(20), FT(20900))
            sat_site = LeafEnvironment(FT(298.15), FT(2000), FT(20), FT(20900))
            beta0_pud = @inferred dark_balanced_psii_partition(fl_pud, jb_photo, dark_site)
            beta0_ku2 = @inferred dark_balanced_psii_partition(fl_conn, jb_photo, dark_site)
            beta0_ku10 = @inferred dark_balanced_psii_partition(fl_ku10, jb_photo, dark_site)
            @test isapprox(beta0_pud, beta0_ku2; rtol) && isapprox(beta0_ku2, beta0_ku10; rtol)
            # Independent BigFloat formula check for β₀ = Φ_P1,max / (Φ_P1,max + η * Φ_P2,max(1))
            S_bf = (big"3.6" / big"260e-6") * (big"179e-3" / (big"3.6" * big"0.27"))
            g_bf = big"0.5" * big"20900" / S_bf
            eta_bf = (big(4) * big(20) + big(9) * g_bf) / (big(4) * big(20) + big(8) * g_bf)
            phi1_bf = big"14.5" / big"15.1"
            phi2_bf = big"4.5" / big"5.1"
            beta0_bf = phi1_bf / (phi1_bf + eta_bf * phi2_bf)
            @test isapprox(beta0_ku2, FT(beta0_bf); rtol)

            # At Q_abs = 0, state_transition_partition equals dark_balanced_psii_partition.
            @test isapprox(state_transition_partition(rule, fl_conn, jb_photo, dark_site), beta0_ku2; rtol)
            dark_exc = @inferred state_transition_excitation(rule, fl_conn, jb_photo, dark_site)
            @test iszero(dark_exc.psii_photons) && iszero(dark_exc.psi_photons)
            @test valid_excitation(dark_exc, zero(FT))
            @test fluorescence(fl_conn, photosystem_flux_inputs(jb_photo, dark_site, dark_exc)).status === FluorescenceDark

            # Puddle (K_U2 = 0) has 0% state transition; K_U2 = 2 gives 17.1%; K_U2 = 10 gives 21.8% (Johnson & Berry 2021, p. 126).
            bsat_pud = state_transition_partition(rule, fl_pud, jb_photo, sat_site)
            bsat_ku2 = state_transition_partition(rule, fl_conn, jb_photo, sat_site)
            bsat_ku10 = state_transition_partition(rule, fl_ku10, jb_photo, sat_site)
            km_bf = big(20) + big(26) * (big(1) + big(20900) / big(17900))
            jc_bf = (big(100) / km_bf) * (big(4) * big(20) + big(9) * g_bf)
            q2sat_bf = big(1) - jc_bf / big(350)
            phi2_sat_bf(ku) = (q2sat_bf * big"4.5" * (big"0.6" + ku)) /
                (big"0.6" * (big"5.1" + ku) + q2sat_bf * big"4.5" * ku)
            bsat_ku2_bf = (q2sat_bf * phi1_bf) / (q2sat_bf * phi1_bf + eta_bf * phi2_sat_bf(big(2)))
            bsat_ku10_bf = (q2sat_bf * phi1_bf) / (q2sat_bf * phi1_bf + eta_bf * phi2_sat_bf(big(10)))
            @test isapprox((beta0_pud - bsat_pud) / beta0_pud, zero(FT); atol = 16 * eps(FT))
            @test isapprox(bsat_ku2, FT(bsat_ku2_bf); rtol)
            @test isapprox(bsat_ku10, FT(bsat_ku10_bf); rtol)
            @test FT(0.169) < (beta0_ku2 - bsat_ku2) / beta0_ku2 < FT(0.172)
            @test FT(0.217) < (beta0_ku10 - bsat_ku10) / beta0_ku10 < FT(0.219)

            # Across light sweep: K_N == 0 below crossover, frozen β == β_sat and K_N > 0 above crossover.
            prev_beta = beta0_ku10
            for q_abs in FT.((1e-4, 5, 50, 200, 400, 500, 650, 1000, 2000))
                s = LeafEnvironment(FT(298.15), q_abs, FT(20), FT(20900))
                b = @inferred state_transition_partition(rule, fl_ku10, jb_photo, s)
                @test b <= prev_beta + 8 * eps(FT)
                prev_beta = b
                exc = @inferred state_transition_excitation(rule, fl_ku10, jb_photo, s)
                @test exc.psii_photons + exc.psi_photons === q_abs
                @test valid_excitation(exc, q_abs)
                c_out = photosynthesize(jb_photo, s, exc)
                f_out = fluorescence(fl_ku10, photosystem_flux_inputs(jb_photo, s, exc))
                @test fluorescence_valid(f_out)
                if c_out.electron_limited < c_out.rubisco_limited - FT(1e-4)
                    # Light-limited regime: dynamic K_N vanishes to numerical precision, and BigFloat balance holds.
                    @test f_out.psii_quenching_rate <= 64 * eps(FT) * FT(15)
                    bf_orc = cytochrome_oracle(jb_photo, fl_ku10, s, exc)
                    @test isapprox(FT(bf_orc.kn), zero(FT); atol = 64 * eps(FT) * FT(15))
                elseif c_out.electron_limited > c_out.rubisco_limited + FT(1e-3)
                    # Rubisco-limited regime: β is frozen at β_sat and K_N > 0 absorbs excess light.
                    @test b == bsat_ku10
                    @test f_out.psii_quenching_rate > FT(0.05)
                end
            end

            # 3. Sterbenz-exact photon conservation across 60+ IEEE-754 binades up to floatmax(FT).
            exponents = FT === Float32 ? (-40:120) : (-300:300)
            for k in exponents
                q_bin = ldexp(FT(1.37591), k)
                isfinite(q_bin) && q_bin > zero(FT) || continue
                for test_rule in (rule, StateTransitionRule(FT(0.2), FT(0.42))) # exercises both β >= 0.5 and β < 0.5 branches!
                    s_bin = LeafEnvironment(T_K, q_bin, FT(28), O2_Pa)
                    exc_bin = state_transition_excitation(test_rule, fl_conn, photo, s_bin)
                    @test exc_bin.psii_photons + exc_bin.psi_photons === q_bin
                    @test valid_excitation(exc_bin, q_bin)
                end
            end
            s_max = LeafEnvironment(T_K, floatmax(FT), FT(28), O2_Pa)
            exc_max = state_transition_excitation(rule, fl_conn, photo, s_max)
            @test exc_max.psii_photons + exc_max.psi_photons === floatmax(FT)
            @test valid_excitation(exc_max, floatmax(FT))

            # 4. LHCII physiological clamping, never-saturating Rubisco (J_c >= Vqmax), Vqmax = 0, Vcmax = 0, and smooth colimitations.
            s_mid = LeafEnvironment(FT(298.15), FT(300), FT(20), FT(20900))
            floor_rule = StateTransitionRule(FT(0.49), FT(0.60))
            @test state_transition_partition(floor_rule, fl_ku10, jb_photo, s_mid) == FT(0.49)
            @test fluorescence_valid(fluorescence(fl_ku10,
                photosystem_flux_inputs(jb_photo, s_mid, state_transition_excitation(floor_rule, fl_ku10, jb_photo, s_mid))))
            s_low = LeafEnvironment(FT(298.15), FT(10), FT(20), FT(20900))
            # Mild PSII starvation (Φ_P2,max < Φ_P2 < 1): valid Z-scheme turnover, incompatible K_N >= 0 inversion.
            mild_starve_rule = StateTransitionRule(FT(0.49), FT(0.50))
            @test fluorescence(fl_ku10,
                photosystem_flux_inputs(jb_photo, s_low, state_transition_excitation(mild_starve_rule, fl_ku10, jb_photo, s_low))).status ===
                FluorescenceIncompatibleRates
            # Severe PSII starvation (Φ_P2 > 1): violates Z-scheme electron-per-photon conservation -> NaN inputs -> InvalidInput.
            severe_starve_rule = StateTransitionRule(FT(0.25), FT(0.45))
            @test fluorescence(fl_ku10,
                photosystem_flux_inputs(jb_photo, s_low, state_transition_excitation(severe_starve_rule, fl_ku10, jb_photo, s_low))).status ===
                FluorescenceInvalidInput

            # Huge Vcmax (J_P700,c > Vqmax): leaf never Rubisco-saturates, K_N remains ~0 even at high light.
            high_vc = setup_cytochrome_stress(FT; vcmax_factor = FT(5))
            s_hi = LeafEnvironment(T_K, FT(1800), FT(28), O2_Pa)
            exc_hi = state_transition_excitation(rule, high_vc.jb_connected, high_vc.cyto_leaf, s_hi)
            fl_hi = fluorescence(high_vc.jb_connected, photosystem_flux_inputs(high_vc.cyto_leaf.photosynthesis, s_hi, exc_hi))
            @test fluorescence_valid(fl_hi) && fl_hi.psii_quenching_rate <= 64 * eps(FT) * FT(15)

            # Concave colimitations (QuadraticLimitation & SerialLimitation) remain K_N >= 0 compatible.
            for lim in (QuadraticLimitation(FT(0.92)), SerialLimitation())
                colim_photo = CytochromePhotosynthesis(;
                    jb_photo.temperature, jb_photo.electron_transport,
                    jb_photo.stoichiometry, limitation = lim, jb_photo.compensation)
                for q_val in FT.((25, 550, 1500))
                    sc = LeafEnvironment(FT(298.15), q_val, FT(20), FT(20900))
                    ec = state_transition_excitation(rule, fl_conn, colim_photo, sc)
                    @test fluorescence_valid(fluorescence(fl_conn, photosystem_flux_inputs(colim_photo, sc, ec)))
                end
            end

            # Zero Vcmax freezes β at dark-balanced β₀; zero Vqmax at Q=0 also equals β₀.
            zero_vc = setup_cytochrome_stress(FT; vcmax_factor = zero(FT))
            @test isapprox(state_transition_partition(rule, zero_vc.jb_connected, zero_vc.cyto_leaf, s_hi),
                dark_balanced_psii_partition(zero_vc.jb_connected, zero_vc.cyto_leaf.photosynthesis, s_hi); rtol)
            zero_vq = setup_cytochrome_stress(FT; transport_factor = zero(FT))
            @test isapprox(state_transition_partition(rule, zero_vq.jb_connected, zero_vq.cyto_leaf,
                LeafEnvironment(T_K, zero(FT), FT(28), O2_Pa)),
                dark_balanced_psii_partition(zero_vq.jb_connected, zero_vq.cyto_leaf.photosynthesis, s_hi); rtol)

            # Invalid / mismatched inputs return NaN without throwing.
            other_ft = FT === Float32 ? Float64 : Float32
            bad_et_photo = CytochromePhotosynthesis(;
                photo.temperature,
                electron_transport = CytochromeElectronParameters(FT(0.525), FT(0.80), SerialLimitation()),
                photo.stoichiometry, photo.limitation, photo.compensation)
            min_et_photo = CytochromePhotosynthesis(;
                photo.temperature,
                electron_transport = CytochromeElectronParameters(FT(0.525), photo.electron_transport.psi_maximum_yield, MinimumLimitation()),
                photo.stoichiometry, photo.limitation, photo.compensation)
            @test isnan(state_transition_partition(rule, fl_conn, bad_et_photo, s_hi))
            @test isnan(state_transition_partition(rule, fl_conn, min_et_photo, s_hi))
            @test isnan(state_transition_partition(rule, fl_conn, photo, LeafEnvironment(other_ft(298.15), other_ft(500), other_ft(28), other_ft(21000))))
            for bad_s in (
                LeafEnvironment(zero(FT), FT(500), FT(28), O2_Pa),
                LeafEnvironment(T_K, -one(FT), FT(28), O2_Pa),
                LeafEnvironment(T_K, FT(500), -one(FT), O2_Pa),
                LeafEnvironment(T_K, FT(500), FT(28), -one(FT)),
                LeafEnvironment(T_K, FT(Inf), FT(28), O2_Pa),
                LeafEnvironment(T_K, FT(NaN), FT(28), O2_Pa),
            )
                @test isnan(dark_balanced_psii_partition(fl_conn, photo, bad_s))
                @test isnan(state_transition_partition(rule, fl_conn, photo, bad_s))
                bad_e = state_transition_excitation(rule, fl_conn, photo, bad_s)
                @test isnan(bad_e.psii_photons) && isnan(bad_e.psi_photons)
            end

            # 5. AD derivative verification (CytoTestDual) vs. 256-bit BigFloat finite differences.
            if FT === Float64
                bf_case = setup_cytochrome_stress(BigFloat)
                bf_rule = StateTransitionRule{BigFloat}()
                h = big"1e-18"
                # Light-limited (PAR = 180) and Rubisco-limited (PAR = 900) points:
                for q0 in (180.0, 900.0)
                    # d/d(absorbed_par)
                    d_site_q = LeafEnvironment(CytoTestDual(T_K), CytoTestDual(q0, 1.0), CytoTestDual(28.0), CytoTestDual(O2_Pa))
                    d_beta_q = @inferred state_transition_partition(rule, fl_conn, photo, d_site_q)
                    d_exc_q = @inferred state_transition_excitation(rule, fl_conn, photo, d_site_q)
                    d_an_q = (@inferred photosynthesize(photo, d_site_q, d_exc_q)).net_assimilation
                    s_hi_bf = LeafEnvironment(BigFloat(T_K), BigFloat(q0) + h, BigFloat(28), BigFloat(O2_Pa))
                    s_lo_bf = LeafEnvironment(BigFloat(T_K), BigFloat(q0) - h, BigFloat(28), BigFloat(O2_Pa))
                    fd_beta_q = (state_transition_partition(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, s_hi_bf) -
                                 state_transition_partition(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, s_lo_bf)) / (2 * h)
                    fd_an_q = (photosynthesize(bf_case.cyto_leaf, s_hi_bf, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, s_hi_bf)).net_assimilation -
                               photosynthesize(bf_case.cyto_leaf, s_lo_bf, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, s_lo_bf)).net_assimilation) / (2 * h)
                    @test isapprox(d_beta_q.d, Float64(fd_beta_q); rtol = 1e-11, atol = 1e-14)
                    @test isapprox(d_an_q.d, Float64(fd_an_q); rtol = 1e-11)

                    # d/d(temperature) & d/d(co2)
                    d_site_t = LeafEnvironment(CytoTestDual(T_K, 1.0), CytoTestDual(q0), CytoTestDual(28.0), CytoTestDual(O2_Pa))
                    d_an_t = (@inferred photosynthesize(photo, d_site_t, state_transition_excitation(rule, fl_conn, photo, d_site_t))).net_assimilation
                    st_hi = LeafEnvironment(BigFloat(T_K) + h, BigFloat(q0), BigFloat(28), BigFloat(O2_Pa))
                    st_lo = LeafEnvironment(BigFloat(T_K) - h, BigFloat(q0), BigFloat(28), BigFloat(O2_Pa))
                    fd_an_t = (photosynthesize(bf_case.cyto_leaf, st_hi, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, st_hi)).net_assimilation -
                               photosynthesize(bf_case.cyto_leaf, st_lo, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, st_lo)).net_assimilation) / (2 * h)
                    @test isapprox(d_an_t.d, Float64(fd_an_t); rtol = 1e-10)

                    d_site_c = LeafEnvironment(CytoTestDual(T_K), CytoTestDual(q0), CytoTestDual(28.0, 1.0), CytoTestDual(O2_Pa))
                    d_an_c = (@inferred photosynthesize(photo, d_site_c, state_transition_excitation(rule, fl_conn, photo, d_site_c))).net_assimilation
                    sc_hi = LeafEnvironment(BigFloat(T_K), BigFloat(q0), BigFloat(28) + h, BigFloat(O2_Pa))
                    sc_lo = LeafEnvironment(BigFloat(T_K), BigFloat(q0), BigFloat(28) - h, BigFloat(O2_Pa))
                    fd_an_c = (photosynthesize(bf_case.cyto_leaf, sc_hi, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, sc_hi)).net_assimilation -
                               photosynthesize(bf_case.cyto_leaf, sc_lo, state_transition_excitation(bf_rule, bf_case.jb_connected, bf_case.cyto_leaf, sc_lo)).net_assimilation) / (2 * h)
                    @test isapprox(d_an_c.d, Float64(fd_an_c); rtol = 1e-11)
                end

                # Parameter AD: d(β_sat)/d(vcmax_ref) with homogeneous dual rate responses & LeafEnvironment{Float64}!
                to_dual_arr(a::Arrhenius{Float64}; dv = 0.0) = Arrhenius{CytoTestDual{Float64}}(
                    CytoTestDual(a.value_ref, dv), CytoTestDual(a.activation_energy),
                    CytoTestDual(a.temperature_ref), CytoTestDual(a.gas_constant))
                to_dual_pa(p::PeakedArrhenius{Float64}; dv = 0.0) = PeakedArrhenius(
                    to_dual_arr(p.arrhenius; dv), CytoTestDual(p.deactivation_energy), CytoTestDual(p.entropy))
                tr0 = photo.temperature
                dual_vc_photo = CytochromePhotosynthesis(;
                    temperature = CytochromeTemperatureResponse(;
                        vcmax = to_dual_pa(tr0.vcmax; dv = 1.0),
                        vqmax = Arrhenius{CytoTestDual{Float64}}(
                            CytoTestDual(tr0.vqmax.site_density * tr0.vqmax.turnover_rate.value_ref),
                            CytoTestDual(tr0.vqmax.turnover_rate.activation_energy),
                            CytoTestDual(tr0.vqmax.turnover_rate.temperature_ref),
                            CytoTestDual(tr0.vqmax.turnover_rate.gas_constant)),
                        respiration = to_dual_pa(tr0.respiration),
                        kc = to_dual_arr(tr0.kc),
                        ko = to_dual_arr(tr0.ko),
                        gamma_star = to_dual_arr(tr0.gamma_star),
                        linear_atp_per_electron = Arrhenius{CytoTestDual{Float64}}(
                            CytoTestDual(tr0.linear_atp_per_electron.value), CytoTestDual(0.0), CytoTestDual(298.15), CytoTestDual(8.314)),
                        cyclic_atp_per_electron = Arrhenius{CytoTestDual{Float64}}(
                            CytoTestDual(tr0.cyclic_atp_per_electron.value), CytoTestDual(0.0), CytoTestDual(298.15), CytoTestDual(8.314)),
                    ),
                    photo.electron_transport, photo.stoichiometry, photo.limitation, photo.compensation,
                )
                s900 = LeafEnvironment(T_K, 900.0, 28.0, O2_Pa)
                db_dvc = @inferred state_transition_partition(rule, fl_conn, dual_vc_photo, s900)
                s900_bf = LeafEnvironment(BigFloat(T_K), BigFloat(900), BigFloat(28), BigFloat(O2_Pa))
                vc0 = BigFloat(tr0.vcmax.arrhenius.value_ref)
                bf_hi = setup_cytochrome_stress(BigFloat; vcmax_factor = (vc0 + h) / vc0)
                bf_lo = setup_cytochrome_stress(BigFloat; vcmax_factor = (vc0 - h) / vc0)
                fd_dvc = (state_transition_partition(bf_rule, bf_hi.jb_connected, bf_hi.cyto_leaf, s900_bf) -
                          state_transition_partition(bf_rule, bf_lo.jb_connected, bf_lo.cyto_leaf, s900_bf)) / (2 * h)
                @test isapprox(db_dvc.d, Float64(fd_dvc); rtol = 1e-11)
            end

            # 6. PAM observation operator: BigFloat oracle, budget closure, PSI leakage theorem, & Float32 cancellation.
            a2_dark = FT(0.85) * beta0_ku2
            a1_dark = FT(0.85) * (one(FT) - beta0_ku2)
            dk_ideal = @inferred dark_pam_levels(obs_ideal, fl_conn, a2_dark, a1_dark)
            dk_pam = @inferred dark_pam_levels(obs_pam, fl_conn, a2_dark, a1_dark)
            @test isapprox(dk_ideal.fv_fm, psii_maximum_photochemical_yield(fl_pud, one(FT)); rtol)
            @test dk_pam.fv_fm < dk_ideal.fv_fm # PSI leakage depresses apparent dark Fv/Fm (~0.81 vs 0.882)
            # Johnson & Berry (2021) p. 117 check: at balanced (α2_dark, α1_dark) and ε_F1 / ε_F2 = 2, PSI is ~39.3% of Fo and ~6.9% of Fm.
            psi_fo_share = (a1_dark * (FT(0.05) / FT(15.1)) * FT(2)) / dk_pam.fo
            psi_fm_share = (a1_dark * (FT(0.05) / FT(15.1)) * FT(2)) / dk_pam.fm
            @test FT(0.385) < psi_fo_share < FT(0.40)
            @test FT(0.065) < psi_fm_share < FT(0.075)

            for q_abs in FT.((10, 250, 900, 1800))
                s = LeafEnvironment(T_K, q_abs, FT(28), O2_Pa)
                res_conn = fluorescence(fl_conn, photosystem_flux_inputs(photo, s))
                res_pud = fluorescence(fl_pud, photosystem_flux_inputs(photo, s))
                a2_s = FT(0.85) * photo.electron_transport.psii_fraction
                a1_s = FT(0.85) * (one(FT) - photo.electron_transport.psii_fraction)

                # Ideal PSII-only detector recovers internal PSII yields & puddle open fraction q₂.
                p_ideal = @inferred pam_fluorescence(obs_ideal, fl_conn, res_conn, a2_s, a1_s)
                p_pud = @inferred pam_fluorescence(obs_ideal, fl_pud, res_pud, a2_s, a1_s)
                @test pam_fluorescence_valid(p_ideal) && pam_fluorescence_valid(p_pud)
                @test isapprox(p_ideal.apparent_photochemical_yield, res_conn.psii_photochemical_yield; rtol)
                @test isapprox(p_ideal.apparent_npq, res_conn.npq; rtol, atol)
                @test isapprox(p_pud.qp, res_pud.psii_open_fraction; rtol)

                # PSI Leakage Invariance Theorem (K_P1 == K_X1, static cross-sections):
                # q_P is exact-invariant to ε_F1, whereas 1 - Fs/Fm', apparent closed lake 1 - q_L, and NPQ_app decrease monotonically!
                prev_phip, prev_closed_ql = p_ideal.apparent_photochemical_yield, one(FT) - p_ideal.ql
                for eps1 in FT.((0.5, 2.0, 25.0))
                    obs_e = PAMFluorescenceObserver(one(FT), eps1)
                    pe = @inferred pam_fluorescence(obs_e, fl_conn, res_conn, a2_s, a1_s)
                    bf_p = pam_bigfloat_oracle(obs_e, fl_conn, res_conn, a2_s, a1_s)
                    @test isapprox(pe.qp, p_ideal.qp; rtol)
                    @test pe.apparent_photochemical_yield < prev_phip
                    @test one(FT) - pe.ql < prev_closed_ql
                    prev_phip, prev_closed_ql = pe.apparent_photochemical_yield, one(FT) - pe.ql
                    @test isapprox(pe.apparent_photochemical_yield + pe.apparent_regulated_heat_yield +
                                   pe.apparent_residual_yield, one(FT); rtol)
                    for f in (:fs, :fm_dark, :fo_dark, :fm_prime, :fo_prime,
                              :apparent_photochemical_yield, :apparent_regulated_heat_yield,
                              :apparent_residual_yield, :apparent_npq, :qp, :ql)
                        @test isapprox(getproperty(pe, f), FT(getproperty(bf_p, f)); rtol, atol)
                    end

                    # Common absorptance scaling invariance of all 6 dimensionless PAM ratios (including extreme Float32 binades).
                    for a_scale in (FT(1e-30), FT(37), FT(1e30))
                        pe_scaled = pam_fluorescence(obs_e, fl_conn, res_conn, a_scale * a2_s, a_scale * a1_s)
                        @test pam_fluorescence_valid(pe_scaled)
                        for f in (:apparent_photochemical_yield, :apparent_regulated_heat_yield,
                                  :apparent_residual_yield, :apparent_npq, :qp, :ql)
                            @test isapprox(getproperty(pe, f), getproperty(pe_scaled, f); rtol)
                        end
                    end
                end

                # State-transition cross-section shift: even when K_N = 0, α2 < α2_dark produces positive apparent NPQ (qT)!
                s_jb = LeafEnvironment(FT(298.15), FT(250), FT(20), FT(20900))
                b_qt = state_transition_partition(rule, fl_conn, jb_photo, s_jb)
                exc_qt = state_transition_excitation(rule, fl_conn, jb_photo, s_jb)
                res_qt = fluorescence(fl_conn, photosystem_flux_inputs(jb_photo, s_jb, exc_qt))
                p_qt = pam_fluorescence(obs_pam, fl_conn, res_qt, FT(0.85) * b_qt, FT(0.85) * (one(FT) - b_qt), a2_dark, a1_dark)
                bf_qt = pam_bigfloat_oracle(obs_pam, fl_conn, res_qt, FT(0.85) * b_qt, FT(0.85) * (one(FT) - b_qt), a2_dark, a1_dark)
                @test res_qt.psii_quenching_rate <= 64 * eps(FT) * FT(15)
                @test p_qt.apparent_npq > FT(0.05) # apparent PAM NPQ captures qT antenna migration!
                @test isapprox(p_qt.apparent_npq, FT(bf_qt.apparent_npq); rtol)
            end

            # Float32 extreme PSI-background & rate-scale cancellation/overflow stress test.
            s_low_kn0 = LeafEnvironment(FT(298.15), FT(50), FT(20), FT(20900))
            exc_kn0 = state_transition_excitation(rule, fl_conn, jb_photo, s_low_kn0)
            res_kn0 = fluorescence(fl_conn, photosystem_flux_inputs(jb_photo, s_low_kn0, exc_kn0))
            b_kn0 = state_transition_partition(rule, fl_conn, jb_photo, s_low_kn0)
            obs_huge = PAMFluorescenceObserver(one(FT), FT(1000))
            p_huge = pam_fluorescence(obs_huge, fl_conn, res_kn0, b_kn0, one(FT) - b_kn0)
            bf_huge = pam_bigfloat_oracle(obs_huge, fl_conn, res_kn0, b_kn0, one(FT) - b_kn0)
            @test isapprox(p_huge.apparent_photochemical_yield, FT(bf_huge.apparent_photochemical_yield); rtol)
            @test isapprox(p_huge.qp, FT(bf_huge.qp); rtol)
            if iszero(res_kn0.psii_quenching_rate)
                @test iszero(p_huge.apparent_npq) && iszero(p_huge.apparent_regulated_heat_yield)
            end

            # Huge rate scales (including near floatmax(FT) where s2_closed_light + kp2 overflows without normalization).
            s_sat_pam = LeafEnvironment(FT(298.15), FT(1800), FT(20), FT(20900))
            p_ref = pam_fluorescence(obs_pam, fl_conn, fluorescence(fl_conn, photosystem_flux_inputs(jb_photo, s_sat_pam)), a2_dark, a1_dark)
            for big_scale in (FT(1e20), floatmax(FT) / FT(40))
                p2, p1 = fl_conn.psii, fl_conn.psi
                fl_big = CytochromeFluorescence(
                    PSIIRateParameters(p2.fluorescence_rate * big_scale, p2.constitutive_heat_rate * big_scale,
                        p2.maximum_photochemical_rate * big_scale, p2.exciton_transfer_rate * big_scale,
                        p2.sustained_quenching_rate * big_scale),
                    PSIRateParameters(p1.fluorescence_rate * big_scale, p1.constitutive_heat_rate * big_scale,
                        p1.maximum_photochemical_rate * big_scale, p1.closed_quenching_rate * big_scale),
                    fl_conn.yield_tolerance,
                )
                @test isapprox(state_transition_partition(rule, fl_big, jb_photo, s_mid),
                    state_transition_partition(rule, fl_conn, jb_photo, s_mid); rtol)
                res_big = fluorescence(fl_big, photosystem_flux_inputs(jb_photo, s_sat_pam))
                p_big = pam_fluorescence(obs_pam, fl_big, res_big, a2_dark, a1_dark)
                @test pam_fluorescence_valid(p_big)
                @test isapprox(p_big.apparent_photochemical_yield, p_ref.apparent_photochemical_yield; rtol)
                @test isapprox(p_big.qp, p_ref.qp; rtol)
                @test isapprox(p_big.ql, p_ref.ql; rtol)
            end

            # Non-lit & invalid PAM status propagation (zero throws).
            res_dark = fluorescence(fl_conn, photosystem_flux_inputs(jb_photo, dark_site, dark_exc))
            p_dark = @inferred pam_fluorescence(obs_pam, fl_conn, res_dark, a2_dark, a1_dark)
            @test !pam_fluorescence_valid(p_dark) && p_dark.status === FluorescenceDark && isnan(p_dark.fs)
            for (ba2, ba1) in ((zero(FT), a1_dark), (-one(FT), a1_dark), (a2_dark, -one(FT)), (FT(NaN), a1_dark))
                @test isnan(dark_pam_levels(obs_pam, fl_conn, ba2, ba1).fo)
                @test pam_fluorescence(obs_pam, fl_conn, res_kn0, ba2, ba1).status === FluorescenceInvalidInput
            end

            # 7. Zero heap allocations across all new scalar entry points.
            s_alloc = LeafEnvironment(T_K, FT(600), FT(28), O2_Pa)
            r_alloc = fluorescence(fl_conn, photosystem_flux_inputs(photo, s_alloc))
            @test additions_allocations(rule, obs_pam, fl_conn, photo, case.cyto_leaf, s_alloc, r_alloc, a2_dark, a1_dark) == 0
        end
    end
end
