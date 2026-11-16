# PSII antenna connectivity (exciton transfer rate KU2) and the PSI detector weight: how
# they set the size of the modeled fluorescence, and when it falls below the dark level Fo.
@testset "PSII connectivity, the PSI weight and the dark level Fo" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            tables = flexas_2002_curves(FT)
            curves, grid = tables["flexas_connectivity_fig2"], tables["flexas_connectivity_rmse"]
            by_name = Dict(r.formulation => r for r in tables["flexas_formulation_summary"])
            @test length(curves) == 84 && length(grid) == 13
            for tbl in (curves, grid), row in tbl
                @test all(v -> v isa FT || v isa Bool, values(row))
            end

            @testset "light responses of PSII Fs / Fo" begin
                # More exciton sharing gives less fluorescence at every irradiance: the same
                # electron flow and open fraction need more quenching.
                @test all(r -> r.g150_ku_0p5 >= r.g150_ku_2 >= r.g150_ku_10, curves)
                @test all(r -> r.g10_ku_0p5 >= r.g10_ku_2 >= r.g10_ku_10, curves)
                peak(field) = maximum(getproperty(r, field) for r in curves)
                @test FT(2.4) < peak(:g150_ku_0p5) < FT(2.7)
                @test FT(1.65) < peak(:g150_ku_2) < FT(1.8)
                @test FT(1.15) < peak(:g150_ku_10) < FT(1.3)
                # Well watered, the signal stays above Fo for every transfer rate; under severe
                # drought it falls below Fo, at lower light the more connected the units are.
                for field in (:g150_ku_0p5, :g150_ku_2, :g150_ku_10)
                    @test minimum(getproperty(r, field) for r in curves) > one(FT)
                end
                # First sampled irradiance below Fo on the 25-unit grid (the continuous crossings
                # are near 443, 351 and 243).
                crossing(field) = first(r.incident_par for r in curves[9:end] if getproperty(r, field) < one(FT))
                @test crossing(:g10_ku_0p5) == FT(450)
                @test crossing(:g10_ku_2) == FT(375)
                @test crossing(:g10_ku_10) == FT(250)
                @test all(getproperty(curves[end], f) < FT(0.6) for f in (:g10_ku_0p5, :g10_ku_2, :g10_ku_10))
            end

            @testset "error grid over KU2 and the PSI weight" begin
                columns = (:rmse_ku_0p5, :rmse_ku_1, :rmse_ku_1p5, :rmse_ku_2, :rmse_ku_3, :rmse_ku_5, :rmse_ku_10)
                @test keys(grid[1]) == (:psi_weight, :psi_share_of_fo, :apparent_fv_fm, :fv_fm_as_measured, columns...)
                at(weight) = only(r for r in grid if r.psi_weight == FT(weight))
                # The benchmark's own cells reproduce its scoreboard.
                tol = FT === Float32 ? FT(1e-5) : FT(1e-12)
                @test isapprox(at(2).rmse_ku_2, by_name["pam2000_combined"].lit_fig2_rmse; rtol = tol)
                @test isapprox(at(0).rmse_ku_2, by_name["pure_psii_combined"].lit_fig2_rmse; rtol = tol)
                @test iszero(at(0).psi_share_of_fo)
                @test FT(0.38) < at(2).psi_share_of_fo < FT(0.40)
                @test issorted([r.psi_share_of_fo for r in grid])
                # The two trade off: the weight with the smallest error falls as KU2 rises, and
                # along that valley the error stays within 0.02 from KU2 = 1 to 3 ns^-1.
                best = map(columns) do column
                    row = argmin(r -> getproperty(r, column), grid)
                    (; weight = row.psi_weight, rmse = getproperty(row, column), weight_as_measured = row.fv_fm_as_measured)
                end
                @test issorted([b.weight for b in best[1:6]]; rev = true)
                valley = [b.rmse for b in best[2:5]]
                @test maximum(valley) < FT(0.07) && maximum(valley) - minimum(valley) < FT(0.025)
                # Off the valley the error grows: without the PSI weight, or with half the
                # transfer rate, it is more than twice that of the carried-over pair.
                @test at(0).rmse_ku_2 > 2 * at(2).rmse_ku_2
                @test at(2).rmse_ku_1 > 2 * at(2).rmse_ku_2
                # Strongly connected units miss the data whatever the PSI weight.
                @test all(r -> r.rmse_ku_10 > FT(0.16), grid)

                # The PSI weight also lowers the dark-adapted apparent Fv/Fm: with c the PSI
                # level relative to PSII Fo, Fv/Fm = (Fm/Fo - 1) / (Fm/Fo + c) for PSII Fm/Fo.
                psii = setup_flexas_2002(FT).fl_cyto.psii
                d0 = psii.fluorescence_rate + psii.constitutive_heat_rate
                fm_over_fo = (d0 + psii.maximum_photochemical_rate) / d0
                for r in grid
                    c = r.psi_share_of_fo / (one(FT) - r.psi_share_of_fo)
                    @test isapprox(r.apparent_fv_fm, (fm_over_fo - one(FT)) / (fm_over_fo + c); rtol = 100 * tol)
                end
                @test issorted([r.apparent_fv_fm for r in grid]; rev = true)
                # Only the weights 2.0 and 2.5 of the grid give the 0.80-0.83 that Flexas et al.
                # measured on all leaves; there KU2 = 2 has the smallest error, and KU2 = 1.5
                # to 3 stay within 0.09 while 1 and 5 exceed 0.11.
                measured = [r for r in grid if r.fv_fm_as_measured]
                @test [r.psi_weight for r in measured] == FT[2.0, 2.5]
                for r in measured
                    errors = [getproperty(r, c) for c in columns]
                    @test argmin(errors) == 4
                    @test all(<(FT(0.09)), errors[3:5])
                    @test all(>(FT(0.11)), errors[[1, 2, 6, 7]])
                end
                # The low-KU2 end of the valley needs weights that Fv/Fm excludes.
                @test !best[2].weight_as_measured && !best[3].weight_as_measured
            end

            @testset "limits of the light-limited rise" begin
                # With an electron flow of q2 times the maximum yield per PSII photon (balanced
                # excitation, light-limited), separate units need no quenching and emit
                # Fs / Fo = q2 + (1 - q2) Fm / Fo, while a lake emits Fs / Fo = 1.
                case = setup_flexas_2002(FT)
                psii = case.fl_cyto.psii
                d0 = psii.fluorescence_rate + psii.constitutive_heat_rate
                fm_over_fo = (d0 + psii.maximum_photochemical_rate) / d0
                maximum_yield = psii.maximum_photochemical_rate / (d0 + psii.maximum_photochemical_rate)
                lake = LakeModelFluorescence(maximum_yield, psii.constitutive_heat_rate / psii.fluorescence_rate,
                    case.fl_cyto.yield_tolerance)
                loose = FT === Float32 ? FT(1e-4) : FT(1e-10)
                for q2 in FT.((0.9, 0.6, 0.3))
                    inputs = CoupledPhotosystemInputs(FT(100), FT(100), FT(100) * q2 * maximum_yield, FT(50), q2)
                    separate = flexas_fs_over_fo(case, with_transfer_rate(case, zero(FT)), case.pam_ideal, inputs)
                    connected = flexas_fs_over_fo(case, case.fl_cyto, case.pam_ideal, inputs)
                    @test isapprox(separate, q2 + (one(FT) - q2) * fm_over_fo; rtol = loose)
                    @test one(FT) < connected < separate
                    r = fluorescence(lake, inputs)
                    @test isapprox(r.fluorescence_yield,
                        lake_model_yields(lake, one(FT), zero(FT)).fluorescence_yield; rtol = loose)
                end
            end

            @testset "Johnson & Berry (2021) Case 9: the yield against Fo" begin
                # Their Fig. 8c (fixed chloroplast CO2, static antenna), normalized by the dark
                # PSII yield: the reference transfer rate ends below Fo at high light.
                base = setup_johnson_berry_2021(FT)
                x = base.settings["experiment"]
                psii = base.fl_ku2.psii
                d0 = psii.fluorescence_rate + psii.constitutive_heat_rate
                fo_yield = psii.fluorescence_rate / (d0 + psii.maximum_photochemical_rate)
                par = collect(x["incident_par_start"]:x["incident_par_step"]:x["incident_par_stop"])[2:end]
                function response(fl, co2)
                    map(par) do q_inc
                        q_abs = x["leaf_absorptance"] * q_inc
                        site = LeafEnvironment(x["temperature_K"], q_abs, co2, x["oxygen_Pa"])
                        r = fluorescence(fl, fluorescence_inputs(base.photo, site, exact_dual_excitation(base.beta_dark, q_abs)))
                        @test is_valid(r)
                        r
                    end
                end
                ratio(results) = [r.psii_fluorescence_yield / fo_yield for r in results]
                reference = response(base.fl_ku2, x["co2_Pa"])
                separate, strong = ratio(response(base.fl_puddle, x["co2_Pa"])), ratio(response(base.fl_ku10, x["co2_Pa"]))
                ku2 = ratio(reference)
                crossing(values) = par[findfirst(<(one(FT)), values)]
                @test FT(4.1) < maximum(separate) < FT(4.4) && minimum(separate) > one(FT)
                # `crossing` is the first sampled irradiance below Fo on the 20-unit grid.
                @test FT(1.5) < maximum(ku2) < FT(1.57) && crossing(ku2) == FT(1640)
                @test FT(1.13) < maximum(strong) < FT(1.19) && crossing(strong) == FT(940)
                @test FT(0.85) < ku2[end] < FT(0.91) && FT(0.69) < strong[end] < FT(0.75)
                # Whatever the transfer rate, the yield is below Fo exactly when the quenching
                # rate exceeds K_P2 (1 - Phi_P2) - (K_F2 + K_D2) Phi_P2.
                for r in reference
                    threshold = psii.maximum_photochemical_rate * (one(FT) - r.psii_photochemical_yield) -
                        d0 * r.psii_photochemical_yield
                    margin = FT === Float32 ? FT(1e-3) : FT(1e-9)
                    abs(r.psii_quenching_rate - threshold) > margin || continue
                    @test (r.psii_fluorescence_yield < fo_yield) == (r.psii_quenching_rate > threshold)
                end
                # A weaker carbon sink (lower chloroplast CO2) lowers the high-light yield
                # further; a stronger one keeps it above Fo.
                low, high = ratio(response(base.fl_ku2, FT(10))), ratio(response(base.fl_ku2, FT(40)))
                @test low[end] < ku2[end] < one(FT) < high[end]
                @test FT(0.55) < low[end] < FT(0.65) && minimum(high) > one(FT)
            end
        end
    end
end
