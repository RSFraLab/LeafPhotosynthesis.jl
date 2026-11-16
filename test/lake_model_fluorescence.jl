using ForwardDiff

# Independent oracle: build explicit first-order rate constants in high precision
# (k_F = 1, k_D = k_DF, k_PM from Phi_PSIImax, k_NPQ = NPQ (k_F + k_D), k_P = q_L k_PM) and
# divide each by their sum, instead of the closed forms of the production code.
function lake_oracle(model, open_fraction, npq)
    kf = big(1)
    kd = BigFloat(model.heat_to_fluorescence_ratio)
    phi_max = BigFloat(model.maximum_photochemical_yield)
    kpm = phi_max * (kf + kd) / (1 - phi_max)
    kn = BigFloat(npq) * (kf + kd)
    kp = BigFloat(open_fraction) * kpm
    total = kf + kd + kn + kp
    return (; photochemical_yield = kp / total, fluorescence_yield = kf / total,
        constitutive_heat_yield = kd / total, regulated_heat_yield = kn / total)
end

lake_allocations(model, inputs) = (fluorescence(model, inputs); @allocated fluorescence(model, inputs))

@testset "Lake-model PSII fluorescence (Gu et al. 2019)" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            rtol = FT === Float32 ? FT(2e-5) : FT(1e-12)
            model = LakeModelFluorescence(FT(0.83), FT(19), eps(FT))
            @test model isa LakeModelFluorescence{FT}
            @test maximum_photochemical_yield(model) === FT(0.83)
            @test LakeModelFluorescence(; maximum_photochemical_yield = FT(0.83),
                heat_to_fluorescence_ratio = FT(19), yield_tolerance = eps(FT)) === model

            @testset "yields from qL and NPQ against explicit rate constants" begin
                for q in FT.((0, 0.05, 0.36, 0.79, 1)), npq in FT.((0, 0.22, 1.77, 12))
                    y = @inferred lake_model_yields(model, q, npq)
                    o = lake_oracle(model, q, npq)
                    @test y isa NamedTuple{(:photochemical_yield, :fluorescence_yield,
                        :constitutive_heat_yield, :regulated_heat_yield),NTuple{4,FT}}
                    for field in keys(o)
                        @test isapprox(getproperty(y, field), FT(getproperty(o, field)); rtol, atol = eps(FT))
                    end
                    @test isapprox(sum(values(y)), one(FT); rtol = 4 * eps(FT))
                end
                # PAM levels of a PSII-only lake: Fv/Fm = Phi_PSIImax and NPQ = Fm/Fm' - 1.
                fo = lake_model_yields(model, one(FT), zero(FT)).fluorescence_yield
                fm = lake_model_yields(model, zero(FT), zero(FT)).fluorescence_yield
                fmp = lake_model_yields(model, zero(FT), FT(1.5)).fluorescence_yield
                @test isapprox((fm - fo) / fm, FT(0.83); rtol = 10 * rtol)
                @test isapprox(fm / fmp - one(FT), FT(1.5); rtol = 10 * rtol)
                @test all(isnan, values(lake_model_yields(model, FT(1.1), zero(FT))))
                @test all(isnan, values(lake_model_yields(model, FT(0.5), FT(-0.1))))
                @test all(isnan, values(lake_model_yields(model, FT(NaN), zero(FT))))
            end

            @testset "NPQ and yields from electron flow and qL" begin
                excitation = FT(425)
                for q in FT.((0.05, 0.36, 0.79, 1)), npq in FT.((0, 0.22, 1.77, 12))
                    o = lake_oracle(model, q, npq)
                    flow = FT(o.photochemical_yield) * excitation
                    inputs = PSIIOpenFractionInputs(excitation, flow, q)
                    r = @inferred fluorescence(model, inputs)
                    @test r isa LakeModelResult{FT}
                    @test is_valid(r) && r.status == FluorescenceLit
                    @test iszero(lake_allocations(model, inputs))
                    # Recovering NPQ divides the yield deficit by Phi_PSII (1 - Phi_0); allow
                    # for its conditioning relative to the working precision.
                    loose = 8 * rtol
                    @test isapprox(r.npq, npq; rtol = loose, atol = loose)
                    @test r.psii_open_fraction === q
                    @test isapprox(r.fluorescence_yield, FT(o.fluorescence_yield); rtol = loose)
                    @test isapprox(r.constitutive_heat_yield, FT(o.constitutive_heat_yield); rtol = loose)
                    @test isapprox(r.regulated_heat_yield, FT(o.regulated_heat_yield); rtol = loose, atol = loose)
                    @test isapprox(r.photon_flux, excitation * r.fluorescence_yield; rtol = 4 * eps(FT))
                    @test isapprox(r.photochemical_yield + r.fluorescence_yield +
                        r.constitutive_heat_yield + r.regulated_heat_yield, one(FT); rtol = 8 * eps(FT))
                    # Gu et al. (2019), Eq. 20: energy balance of PSII in terms of the emitted flux.
                    @test isapprox(flow,
                        excitation - (one(FT) + r.npq) * (one(FT) + FT(19)) * r.photon_flux;
                        rtol = loose, atol = 8 * eps(FT) * excitation)
                    # Eq. 21: the electron flow follows from the emitted flux and qL alone, for
                    # any escape probability and PSII share of the observed flux.
                    for (escape, share) in ((one(FT), one(FT)), (FT(0.8), one(FT)), (FT(0.4), FT(0.75)))
                        observed = r.photon_flux * escape / share
                        back = @inferred electron_flow_from_fluorescence(
                            model, ObservedFluorescenceInputs(observed, escape, share, q))
                        @test back isa FT
                        @test isapprox(back, flow; rtol = loose)
                    end
                end
            end

            @testset "unquenched boundary" begin
                # A flow equal to the unquenched yield of its open fraction is valid, with
                # NPQ = 0 or a roundoff-sized positive value, for every open fraction.
                for q in FT(0.001):FT(0.001):one(FT)
                    phi = lake_model_yields(model, q, zero(FT)).photochemical_yield
                    r = fluorescence(model, PSIIOpenFractionInputs(one(FT), phi, q))
                    @test r.status == FluorescenceLit
                    @test zero(FT) <= r.npq < 64 * eps(FT) / q
                    below = fluorescence(model, PSIIOpenFractionInputs(one(FT), prevfloat(phi), q))
                    @test below.status == FluorescenceLit && below.npq >= zero(FT)
                end
                # Nothing overflows for the largest representable k_DF: the fluorescence yield
                # underflows to zero, the yields still sum to one, and zero flux gives zero flow.
                extreme = LakeModelFluorescence(FT(0.83), floatmax(FT) / 4, eps(FT))
                phi = lake_model_yields(extreme, FT(0.5), one(FT)).photochemical_yield
                r = fluorescence(extreme, PSIIOpenFractionInputs(one(FT), phi, FT(0.5)))
                @test r.status == FluorescenceLit
                @test isapprox(r.photochemical_yield + r.fluorescence_yield +
                    r.constitutive_heat_yield + r.regulated_heat_yield, one(FT); rtol = 8 * eps(FT))
                @test isapprox(sum(values(lake_model_yields(extreme, FT(0.5), one(FT)))), one(FT); rtol = 8 * eps(FT))
                @test iszero(electron_flow_from_fluorescence(extreme,
                    ObservedFluorescenceInputs(zero(FT), one(FT), one(FT), FT(0.5))))
                @test iszero(electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(FT(2), FT(0.8), one(FT), zero(FT))))
            end

            @testset "the open fraction of PAM fluorometry" begin
                # Kramer et al. (2004): qL = (Fm' - Fs) / (Fm' - Fo') * Fo' / Fs, evaluated on
                # the fluorescence levels of the model, returns the input open fraction.
                q, npq = FT(0.41), FT(1.1)
                fs = lake_model_yields(model, q, npq).fluorescence_yield
                fmp = lake_model_yields(model, zero(FT), npq).fluorescence_yield
                fop = lake_model_yields(model, one(FT), npq).fluorescence_yield
                @test isapprox((fmp - fs) / (fmp - fop) * fop / fs, q; rtol = 100 * rtol)
            end

            @testset "status codes" begin
                dark = fluorescence(model, PSIIOpenFractionInputs(zero(FT), zero(FT), one(FT)))
                @test dark.status == FluorescenceDark && is_valid(dark)
                @test iszero(dark.photon_flux) && isnan(dark.fluorescence_yield)
                for bad in (
                    PSIIOpenFractionInputs(FT(-1), zero(FT), one(FT)),
                    PSIIOpenFractionInputs(FT(100), FT(-1), one(FT)),
                    PSIIOpenFractionInputs(FT(100), FT(10), FT(1.2)),
                    PSIIOpenFractionInputs(FT(100), FT(10), FT(NaN)),
                    PSIIOpenFractionInputs(FT(NaN), FT(10), FT(0.5)),
                    PSIIOpenFractionInputs(zero(FT), FT(10), FT(0.5)),
                )
                    r = fluorescence(model, bad)
                    @test r.status == FluorescenceInvalidInput && !is_valid(r) && isnan(r.photon_flux)
                end
                # More electron flow than the open centers allow without quenching, closed
                # centers, and no electron flow in the light (NPQ undetermined or unbounded).
                limit = FT(0.5 * 0.83 / (0.17 + 0.5 * 0.83))
                for bad in (
                    PSIIOpenFractionInputs(FT(100), FT(100) * limit * FT(1.01), FT(0.5)),
                    PSIIOpenFractionInputs(FT(100), FT(10), zero(FT)),
                    PSIIOpenFractionInputs(FT(100), zero(FT), FT(0.5)),
                    PSIIOpenFractionInputs(FT(100), FT(100), one(FT)),
                )
                    r = fluorescence(model, bad)
                    @test r.status == FluorescenceIncompatibleRates && !is_valid(r)
                end
                # Roundoff above the unquenched yield is tolerated and flagged; NPQ is zero there.
                tolerant = LakeModelFluorescence(FT(0.83), FT(19), FT(1e-4))
                edge = fluorescence(tolerant, PSIIOpenFractionInputs(FT(100), FT(100) * (limit + FT(5e-5)), FT(0.5)))
                @test edge.status == FluorescenceLit && edge.boundary_roundoff && iszero(edge.npq)
                @test isapprox(edge.photochemical_yield + edge.fluorescence_yield +
                    edge.constitutive_heat_yield + edge.regulated_heat_yield, one(FT); rtol = 8 * eps(FT))
                @test isnan(electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(FT(1), zero(FT), one(FT), FT(0.5))))
                @test isnan(electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(FT(1), FT(1.2), one(FT), FT(0.5))))
                @test isnan(electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(FT(-1), FT(0.8), one(FT), FT(0.5))))
                @test iszero(electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(zero(FT), FT(0.8), one(FT), FT(0.5))))
            end

            @testset "construction domain and precision" begin
                @test_throws ArgumentError LakeModelFluorescence(zero(FT), FT(19), eps(FT))
                @test_throws ArgumentError LakeModelFluorescence(one(FT), FT(19), eps(FT))
                @test_throws ArgumentError LakeModelFluorescence(FT(0.83), FT(-1), eps(FT))
                @test_throws ArgumentError LakeModelFluorescence(FT(0.83), FT(Inf), eps(FT))
                @test_throws ArgumentError LakeModelFluorescence(FT(0.83), FT(19), FT(0.2))
                @test_throws ArgumentError LakeModelFluorescence(FT(0.83), FT(19), FT(-1e-6))
                @test_throws MethodError LakeModelFluorescence(1, 19, 0)
                Other = FT === Float32 ? Float64 : Float32
                @test_throws MethodError fluorescence(model, PSIIOpenFractionInputs(Other(100), Other(10), Other(0.5)))
                @test_throws MethodError lake_model_yields(model, Other(0.5), Other(0))
                @test_throws MethodError electron_flow_from_fluorescence(model,
                    ObservedFluorescenceInputs(Other(1), Other(0.8), Other(1), Other(0.5)))
                @test PSIIOpenFractionInputs(PSIIFluxInputs(FT(100), FT(10)), FT(0.5)) ===
                    PSIIOpenFractionInputs(FT(100), FT(10), FT(0.5))
            end

            @testset "batches" begin
                inputs = [PSIIOpenFractionInputs(FT(q2), FT(0.3) * FT(q2), FT(0.6)) for q2 in (0, 50, 400, 900)]
                results = fluorescence.(model, inputs)
                @test eltype(results) === LakeModelResult{FT}
                @test results[1].status == FluorescenceDark
                @test all(r -> r.status == FluorescenceLit, results[2:end])
                @test isnothing(fluorescence(NoFluorescence(), inputs[2]))
            end

            @testset "host coupling" begin
                # A C3 host whose low-light PSII yield does not exceed the lake-model maximum.
                host = ElectronTransportParameters(FT(0.5), FT(0.80), MinimumLimitation())
                @test validate_fluorescence_coupling(model, host) === nothing
                too_efficient = ElectronTransportParameters(FT(0.5), FT(0.85), MinimumLimitation())
                @test_throws ArgumentError validate_fluorescence_coupling(model, too_efficient)
            end

            @testset "lake limit of connected PSII units (Johnson & Berry 2021)" begin
                # An independent implementation of the same physics: with very fast exciton
                # transfer between PSII units the connected-units model is a lake, so it must
                # return the PSII yields and NPQ of the lake model for the same inputs.
                case = setup_flexas_2002(FT)
                psii = case.fl_cyto.psii
                transfer = FT === Float32 ? FT(1e5) : FT(1e9)
                lake_limit = CytochromeFluorescence(;
                    psii = PSIIRateParameters(psii.fluorescence_rate, psii.constitutive_heat_rate,
                        psii.maximum_photochemical_rate, transfer, psii.sustained_quenching_rate),
                    psi = case.fl_cyto.psi,
                    yield_tolerance = case.fl_cyto.yield_tolerance,
                )
                lake = LakeModelFluorescence(psii_maximum_photochemical_yield(case.fl_cyto),
                    psii.constitutive_heat_rate / psii.fluorescence_rate, case.fl_cyto.yield_tolerance)
                x = case.settings["experiment"]
                tol = FT === Float32 ? FT(2e-3) : FT(1e-7)
                for ppfd in FT.((50, 400, 1200)), co2 in FT.((8, 25))
                    q_abs = x["leaf_absorptance"] * ppfd
                    site = LeafEnvironment(x["temperature_K"], q_abs, co2, x["ambient_oxygen_Pa"])
                    inputs = fluorescence_inputs(case.cyto_photo, site, exact_dual_excitation(case.beta_dark, q_abs))
                    connected = fluorescence(lake_limit, inputs)
                    r = @inferred fluorescence(lake, inputs)
                    @test is_valid(connected) && is_valid(r)
                    @test isapprox(r.fluorescence_yield, connected.psii_fluorescence_yield; rtol = tol)
                    @test isapprox(r.regulated_heat_yield, connected.psii_regulated_heat_yield; rtol = tol, atol = tol)
                    @test isapprox(r.npq, connected.npq; rtol = tol, atol = tol)
                    @test isapprox(r.photon_flux, connected.psii_photon_flux; rtol = tol)
                    @test r.psii_open_fraction == connected.psii_open_fraction
                end
            end
        end
    end

    @testset "sensitivity to the open fraction" begin
        # At fixed electron flow, Phi_F = Phi_PSII (1 - Phi_PSIImax) / (qL Phi_PSIImax (1 + k_DF)),
        # so d ln Phi_F / d ln qL = -1, and d ln J / d ln qL = 1 in Eq. 21: a relative error in
        # qL passes undamped. At fixed NPQ instead (Eq. 14), d ln Phi_F / d ln qL = -Phi_PSII.
        model = LakeModelFluorescence(0.83, 19.0, eps())
        yield_at(q) = fluorescence(model, PSIIOpenFractionInputs(oftype(q, 400.0), oftype(q, 120.0), q)).fluorescence_yield
        q = 0.45
        @test isapprox(ForwardDiff.derivative(yield_at, q), -yield_at(q) / q; rtol = 1e-12)
        flow_at(q) = electron_flow_from_fluorescence(model,
            ObservedFluorescenceInputs(oftype(q, 2.0), oftype(q, 0.8), oftype(q, 1.0), q))
        @test isapprox(ForwardDiff.derivative(flow_at, q), flow_at(q) / q; rtol = 1e-12)
        at_fixed_npq(q) = lake_model_yields(model, q, oftype(q, 1.2))
        y = at_fixed_npq(q)
        @test isapprox(ForwardDiff.derivative(x -> at_fixed_npq(x).fluorescence_yield, q),
            -y.photochemical_yield * y.fluorescence_yield / q; rtol = 1e-12)
    end
end
