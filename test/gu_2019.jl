@testset "Gu et al. (2019) figures and the Flexas et al. (2002) comparison" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            ref = setup_gu_2019(FT)
            tables = gu_2019_curves(FT)
            vitis, sif = tables["gu_vitis_light_response"], tables["gu_sif_assimilation"]
            flexas, identity, summary = tables["gu_flexas_fig2"], tables["gu_flexas_identity"], tables["gu_flexas_summary"]
            inset = tables["gu_flexas_inset"]
            for tbl in (vitis, sif, flexas, inset, identity), row in tbl
                @test all(v -> v isa FT || v isa Bool, values(row))
            end
            @test all(row -> all(v -> v isa FT || v isa String, values(row)), summary)

            @testset "Fig. 4c: yields from the measured qL and NPQ of Fig. 1c" begin
                @test length(vitis) == 8
                phi_max = ref.paper_lake.maximum_photochemical_yield
                k_df = ref.paper_lake.heat_to_fluorescence_ratio
                for row in vitis
                    # Eq. 16 reproduces the plotted Phi_PSII within the digitization error.
                    @test abs(row.photochemical_yield - row.paper_photochemical_yield) < FT(0.002)
                    # Eq. 21 returns that yield from the fluorescence yield and qL alone.
                    @test isapprox(row.yield_recovered_from_fluorescence, row.photochemical_yield; rtol = FT(1e-5))
                    # The plotted Phi_SIF lies 26-33 % above Eq. 14 and matches Eq. 14 without
                    # (1 - Phi_PSIImax) in the bracket, which does not conserve energy with the
                    # Phi_PSII of the same panel.
                    @test FT(1.25) < row.paper_fluorescence_yield / row.fluorescence_yield < FT(1.34)
                    combination = (one(FT) - phi_max) * row.measured_npq + phi_max * row.measured_open_fraction
                    without_term = (one(FT) - phi_max) / ((one(FT) + k_df) * combination)
                    @test isapprox(row.paper_fluorescence_yield, without_term; rtol = FT(0.015))
                    plotted_total = row.paper_photochemical_yield + (one(FT) + k_df) *
                        (one(FT) + row.measured_npq) * row.paper_fluorescence_yield
                    @test plotted_total > FT(1.05)
                end
                # The muted light response of the fluorescence yield, against Phi_PSII.
                yields = [r.fluorescence_yield for r in vitis]
                @test maximum(yields) / minimum(yields) < FT(1.25)
                @test vitis[1].photochemical_yield / vitis[end].photochemical_yield > FT(1.9)
            end

            @testset "Fig. 2a: assimilation against fluorescence" begin
                f2 = ref.settings["figure2"]
                dig = ref.settings["digitized_figure2a"]
                last_row = sif[end]
                @test last_row.sif == dig["sif"]
                modeled = (last_row.gpp_rsb_025, last_row.gpp_rsb_050, last_row.gpp_rsb_075,
                    last_row.gpp_rsb_100, last_row.gpp_rsb_125)
                for (m, d) in zip(modeled, FT.(dig["gpp"]))
                    @test abs(m - d) < FT(0.1)
                end
                @test iszero(sif[1].linear_electron_flow)
                @test isapprox(sif[1].gpp_rsb_100, zero(FT); atol = FT(1e-4))
                # Independent oracle for the package's gas exchange: Eq. 4 at
                # C_c = C_a - (GPP - R_d) r is a quadratic in GPP, solved here in high precision.
                function oracle(flow, resistance)
                    ca, gamma, rd = big(f2["ambient_co2_Pa"]), big(f2["gamma_star_Pa"]), big(f2["respiration"])
                    ec, a = big(f2["electrons_per_carboxylation"]), big(f2["oxygenation_cost"])
                    j, r = big(flow), big(resistance)
                    supply = ca + rd * r
                    b = ec * (supply + a * gamma) + j * r
                    c = j * (supply - gamma)
                    return 2c / (b + sqrt(b^2 - 4 * ec * r * c))
                end
                tol = FT === Float32 ? FT(2e-4) : FT(1e-9)
                keys_and_resistances = zip((:gpp_rsb_025, :gpp_rsb_050, :gpp_rsb_075, :gpp_rsb_100, :gpp_rsb_125),
                    FT.(f2["boundary_stomatal_resistances"]))
                for row in sif[2:end], (key, r_sb) in keys_and_resistances
                    expected = oracle(row.linear_electron_flow, r_sb + f2["wall_resistance"])
                    @test isapprox(getproperty(row, key), FT(expected); rtol = tol)
                end
            end

            @testset "Flexas et al. (2002) Figure 2" begin
                @test length(flexas) == 49 && length(inset) == 6
                @test all(r -> !isnan(r.lake_measured_open_fraction), flexas)
                by_name = Dict(r.formulation => r for r in summary)
                jb, model_q, measured_q = by_name["johnson_berry_pam2000"],
                    by_name["lake_model_open_fraction"], by_name["lake_measured_open_fraction"]
                @test jb.lit_fig2_rmse < FT(0.065)
                @test jb.lit_fig2_rmse < model_q.lit_fig2_rmse < measured_q.lit_fig2_rmse
                @test jb.inset_rmse < model_q.inset_rmse < measured_q.inset_rmse
                @test FT(0.23) < model_q.lit_fig2_rmse < FT(0.26)
                @test FT(0.40) < measured_q.lit_fig2_rmse < FT(0.46)
                # PSII alone against a signal that contains PSI scores worse for both models.
                @test by_name["johnson_berry_psii"].lit_fig2_rmse > jb.lit_fig2_rmse
                @test by_name["lake_model_open_fraction_psii"].lit_fig2_rmse > model_q.lit_fig2_rmse

                # The observer adds a PSI level that is constant here (K_X1 = K_P1): with c the
                # PSI share of dark Fo relative to PSII, detector Fs/Fo = (PSII Fs/Fo + c) / (1 + c).
                case = ref.flexas
                c = (case.dark_pam.fo - case.dark_psii.fo) / case.dark_psii.fo
                @test FT(0.6) < c < FT(0.67)
                tol = FT === Float32 ? FT(2e-5) : FT(1e-10)
                @test all(r -> isapprox(r.lake_model_open_fraction,
                    (r.lake_model_open_fraction_psii + c) / (one(FT) + c); rtol = tol), flexas)
                @test all(r -> isapprox(r.johnson_berry_pam2000,
                    (r.johnson_berry_psii + c) / (one(FT) + c); rtol = 50 * tol), flexas)

                # Connectivity experiment: the same electron flow and the same open-center
                # population q2 in a lake instead of connected units. Where light limits, the
                # Johnson & Berry flow is q2 Phi_PSIImax eta_0 / eta(C_c) per PSII photon (eta:
                # PSI/PSII electron ratio, eta_0 its value at the CO2 of the dark-balanced
                # partition), for which a lake emits Fs/Fo = eta_0 / eta(C_c), close to one.
                light_limited = [r for r in flexas if r.conductance_mmol == FT(150) && r.incident_par < FT(130)]
                @test length(light_limited) == 4
                @test all(r -> abs(r.lake_model_open_fraction_psii - one(FT)) < FT(0.01), light_limited)
                @test all(r -> r.observed_fs_over_fo > FT(1.05), light_limited)
                @test maximum(r.lake_model_open_fraction_psii for r in flexas) < FT(1.01)
                # The lake needs more quenching than connected units for the same flow and q2
                # (both are intrinsic PSII Stern-Volmer indices).
                @test all(r -> r.lake_model_npq > r.johnson_berry_psii_npq, flexas)
                # The two antenna models share the rate physics: given the lake qL that a
                # PSII-only fluorometer would report for the connected units, instead of the
                # population q2, the lake relations return the connected-units result.
                loose = FT === Float32 ? FT(2e-3) : FT(1e-8)
                for r in flexas
                    @test r.model_pam_lake_open_fraction < r.model_open_fraction
                    lake = fluorescence(ref.host_lake, PSIIOpenFractionInputs(
                        one(FT), r.psii_photochemical_yield, r.model_pam_lake_open_fraction))
                    @test isapprox(lake.npq, r.johnson_berry_psii_npq; rtol = loose, atol = loose)
                    @test isapprox(fs_over_fo(ref.host_lake, lake), r.johnson_berry_psii; rtol = loose)
                end

                # Figure 2 Inset, one irradiance: every formulation shares the electron flow and
                # falls with drought, as observed.
                @test issorted([r.conductance_mmol for r in inset])
                dry, wet = inset[1], inset[end]
                for field in (:observed_fs_over_fo, :johnson_berry_pam2000, :lake_model_open_fraction,
                    :lake_measured_open_fraction)
                    @test issorted([getproperty(r, field) for r in inset])
                end
                drop(field) = getproperty(wet, field) / getproperty(dry, field)
                @test FT(1.55) < drop(:observed_fs_over_fo) < FT(1.65)
                @test FT(1.65) < drop(:johnson_berry_pam2000) < FT(1.75)
                @test FT(1.35) < drop(:lake_model_open_fraction) < FT(1.45)
                @test drop(:lake_measured_open_fraction) > FT(2.5)
                # But the level is missed: no lake variant rises above Fo at this irradiance.
                @test all(r -> r.lake_model_open_fraction < one(FT) && r.lake_measured_open_fraction < one(FT), inset)
                @test wet.observed_fs_over_fo > FT(1.4)
                # Conditional on the apparent electron flow of the benchmark, the lake relations
                # need an open fraction below the transferred one at every stage, falling with
                # drought; with the transferred one Eq. 21 returns from the observed Fs/Fo
                # 1.8 to 3.4 times that apparent flow (a comparison between two models: the
                # electron transport measured by Flexas et al. is not digitized here).
                @test issorted([r.implied_open_fraction for r in inset])
                @test all(r -> r.implied_open_fraction < FT(0.6) * r.measured_open_fraction, inset)
                @test all(r -> FT(1.8) < r.electron_flow_from_observed_fs / r.apparent_electron_flow < FT(3.4), inset)
                @test all(r -> zero(FT) < r.implied_open_fraction < one(FT), flexas)
                @test all(r -> r.apparent_photochemical_yield < r.psii_photochemical_yield, flexas)
            end

            @testset "Fig. 1c inputs against the well-watered Flexas curve" begin
                inside = [r for r in identity if !r.outside_g150_range]
                @test length(inside) == 7
                rmse(field) = flexas_rmse(FT, [r.lake_measured_inputs for r in inside],
                    [getproperty(r, field) for r in inside])
                @test rmse(:observed_g150) < FT(0.05)
                @test rmse(:observed_g150) < rmse(:observed_g80)
                peak = maximum(r.lake_measured_inputs for r in identity)
                @test FT(1.40) < peak < FT(1.46)
            end
        end
    end
end
