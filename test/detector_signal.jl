using ForwardDiff

# What a fluorometer reports against the PSII yield: the detector level of Johnson & Berry
# (2021, Eq. 38) for any PSII fluorescence model, and its use in the Flexas benchmark.
@testset "Detector signal: PSII yield plus a constant PSI level" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            case = setup_flexas_2002(FT)
            observer, ideal = case.pam_pam2000, case.pam_ideal
            a2, a1 = case.alpha2_dark, case.alpha1_dark

            @testset "detector_level" begin
                inputs = flexas_photosystem_inputs(case, FT(600))
                r = fluorescence(case.fl_cyto, inputs)
                pam = pam_fluorescence(observer, case.fl_cyto, r, a2, a1)
                level = @inferred detector_level(observer, a2, r.psii_fluorescence_yield, a1, r.psi_fluorescence_yield)
                @test level isa FT
                # The same expression as the PAM observation of the Cytochrome model.
                @test level == pam.fs
                @test detector_level(ideal, a2, r.psii_fluorescence_yield, a1, r.psi_fluorescence_yield) ==
                    pam_fluorescence(ideal, case.fl_cyto, r, a2, a1).fs
                # Dark levels from the rate constants: Fo with open centers, Fm with closed ones.
                psii, psi = case.fl_cyto.psii, case.fl_cyto.psi
                d2 = psii.fluorescence_rate + psii.constitutive_heat_rate + psii.sustained_quenching_rate
                psi_open = psi.fluorescence_rate /
                    (psi.fluorescence_rate + psi.constitutive_heat_rate + psi.maximum_photochemical_rate)
                psi_closed = psi.fluorescence_rate /
                    (psi.fluorescence_rate + psi.constitutive_heat_rate + psi.closed_quenching_rate)
                tol = 8 * eps(FT)
                @test isapprox(detector_level(observer, a2, psii.fluorescence_rate / (d2 + psii.maximum_photochemical_rate), a1, psi_open),
                    case.dark_pam.fo; rtol = tol)
                @test isapprox(detector_level(observer, a2, psii.fluorescence_rate / d2, a1, psi_closed),
                    case.dark_pam.fm; rtol = tol)
                # PSII alone: no PSI absorbed fraction or yield.
                @test detector_level(observer, a2, FT(0.01), zero(FT), zero(FT)) ==
                    (a2 * FT(0.01)) * observer.psii_sensitivity
                for bad in ((zero(FT), FT(0.01), a1, FT(0.003)), (a2, FT(-0.01), a1, FT(0.003)),
                    (a2, FT(0.01), FT(-0.1), FT(0.003)), (a2, FT(0.01), a1, FT(NaN)), (a2, FT(Inf), a1, FT(0.003)))
                    @test isnan(detector_level(observer, bad...))
                end
                Other = FT === Float32 ? Float64 : Float32
                @test_throws MethodError detector_level(observer, Other(0.4), Other(0.01), Other(0.4), Other(0.003))
                # Derivatives: the level is linear in each yield, and validity is decided on
                # values, so a zero PSI yield stays valid whatever the sign of its seed.
                @test ForwardDiff.derivative(y -> detector_level(observer, a2, y, a1, FT(0.003)), FT(0.01)) ≈
                    a2 * observer.psii_sensitivity
                for seed in (one(FT), -one(FT))
                    zero_yield = ForwardDiff.Dual{:seed}(zero(FT), seed)
                    level = detector_level(observer, ForwardDiff.Dual{:seed}(a2, zero(FT)),
                        ForwardDiff.Dual{:seed}(FT(0.01), zero(FT)), ForwardDiff.Dual{:seed}(a1, zero(FT)), zero_yield)
                    @test ForwardDiff.value(level) == (a2 * FT(0.01)) * observer.psii_sensitivity
                    @test ForwardDiff.partials(level, 1) ≈ seed * a1 * observer.psi_sensitivity
                end
                @test_throws MethodError detector_level(observer, ForwardDiff.Dual{:seed}(Other(0.4), one(Other)),
                    ForwardDiff.Dual{:seed}(Other(0.01), zero(Other)), ForwardDiff.Dual{:seed}(Other(0.4), zero(Other)),
                    ForwardDiff.Dual{:seed}(Other(0.003), zero(Other)))
                @test iszero((detector_level(observer, a2, FT(0.01), a1, FT(0.003));
                    @allocated detector_level(observer, a2, FT(0.01), a1, FT(0.003))))
            end

            tables = flexas_2002_curves(FT)
            psii_curves, detector_curves = tables["flexas_connectivity_fig2"], tables["flexas_connectivity_detector_fig2"]
            points, summary = tables["flexas_detector_fig2"], tables["flexas_detector_summary"]
            @test length(detector_curves) == length(psii_curves) == 84 && length(points) == 49
            for tbl in (detector_curves, points), row in tbl
                @test all(v -> v isa FT, values(row))
            end
            # With c the dark PSI level relative to the dark PSII level of a model, a constant
            # PSI level turns a PSII ratio R into (R + c) / (1 + c).
            through_detector(ratio, c) = (ratio + c) / (one(FT) + c)
            c_cytochrome = case.psi_dark_level / case.dark_psii.fo
            c_empirical = case.psi_dark_level / (a2 * case.empirical_fo_yield * observer.psii_sensitivity)
            tol = FT === Float32 ? FT(2e-5) : FT(1e-10)

            @testset "connectivity curves through the detector" begin
                columns = keys(psii_curves[1])[2:end]
                @test keys(detector_curves[1]) == keys(psii_curves[1])
                for (p, d) in zip(psii_curves, detector_curves), column in columns
                    @test isapprox(getproperty(d, column), through_detector(getproperty(p, column), c_cytochrome); rtol = tol)
                    # The detector signal is below Fo exactly where PSII is, and closer to it.
                    @test (getproperty(d, column) < one(FT)) == (getproperty(p, column) < one(FT))
                    @test abs(getproperty(d, column) - one(FT)) <= abs(getproperty(p, column) - one(FT))
                end
                peak(rows, column) = maximum(getproperty(r, column) for r in rows)
                @test FT(1.42) < peak(detector_curves, :g150_ku_2) < FT(1.46) < FT(1.7) < peak(psii_curves, :g150_ku_2)
                @test FT(1.9) < peak(detector_curves, :g150_ku_0p5) < FT(2.0)
                @test FT(1.1) < peak(detector_curves, :g150_ku_10) < FT(1.16)
            end

            @testset "PSII yield against detector signal at the digitized points" begin
                for r in points
                    @test isapprox(r.cytochrome_detector, through_detector(r.cytochrome_psii, c_cytochrome); rtol = 50 * tol)
                    @test isapprox(r.cotton_detector, through_detector(r.cotton_psii, c_empirical); rtol = 50 * tol)
                    @test isapprox(r.drought_detector, through_detector(r.drought_psii, c_empirical); rtol = 50 * tol)
                    # The PSI level is constant; its share of Fs is c / (R + c), 39 % at Fo,
                    # less where PSII fluorescence is above Fo and more where it is quenched.
                    @test isapprox(r.psi_share_of_fs, c_cytochrome / (r.cytochrome_psii + c_cytochrome); rtol = 50 * tol)
                    # Its share of Fm' exceeds the 7 % of the dark-adapted Fm as soon as NPQ acts.
                    @test r.psi_share_of_fm_prime > case.psi_dark_level / case.dark_pam.fm
                    @test r.psi_share_of_fm_prime < r.psi_share_of_fs
                end
                @test FT(0.385) < c_cytochrome / (one(FT) + c_cytochrome) < FT(0.39)
                @test FT(0.065) < case.psi_dark_level / case.dark_pam.fm < FT(0.075)
                @test FT(0.26) < minimum(r.psi_share_of_fs for r in points) < FT(0.29)
                @test FT(0.49) < maximum(r.psi_share_of_fs for r in points) < FT(0.53)
                @test FT(0.45) < maximum(r.psi_share_of_fm_prime for r in points) < FT(0.50)

                by_name = Dict(r.formulation => r for r in summary)
                scoreboard = Dict(r.formulation => r for r in tables["flexas_formulation_summary"])
                rtol = FT === Float32 ? FT(1e-5) : FT(1e-12)
                # The PSII-only columns reproduce the benchmark's scoreboard.
                @test isapprox(by_name["cytochrome"].psii_only_rmse, scoreboard["pure_psii_combined"].lit_fig2_rmse; rtol)
                @test isapprox(by_name["cytochrome"].detector_rmse, scoreboard["pam2000_combined"].lit_fig2_rmse; rtol)
                @test isapprox(by_name["cotton"].psii_only_rmse, scoreboard["cotton_combined"].lit_fig2_rmse; rtol)
                @test isapprox(by_name["drought"].psii_only_rmse, scoreboard["drought_combined"].lit_fig2_rmse; rtol)
                # Through the observer the Cytochrome model improves almost threefold. Adding a
                # prescribed PSI level to the empirical fits (a sensitivity case: their
                # constants come from fluorometer signals in which PSI was not corrected out)
                # changes little overall, better under drought and worse well watered.
                @test by_name["cytochrome"].detector_rmse < FT(0.4) * by_name["cytochrome"].psii_only_rmse
                @test FT(0.8) < by_name["cotton"].detector_rmse / by_name["cotton"].psii_only_rmse < one(FT)
                @test abs(by_name["drought"].detector_rmse - by_name["drought"].psii_only_rmse) < FT(0.005)
                for name in ("cotton", "drought")
                    @test by_name[name].detector_rmse > FT(1.8) * by_name["cytochrome"].detector_rmse
                    @test by_name[name].detector_rmse_g150 > by_name[name].psii_only_rmse_g150
                    @test by_name[name].detector_rmse_g10 < by_name[name].psii_only_rmse_g10
                end
            end
        end
    end
end
