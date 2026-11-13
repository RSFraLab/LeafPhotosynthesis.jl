include("../examples/flexas_2002.jl")

@testset "Flexas et al. (2002) water-stress fluorescence benchmark" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            ref = setup_flexas_2002(FT)
            # Walz PAM-2000 (>700 nm, eps_F1 / eps_F2 = 2.0) dark-adapted apparent Fv/Fm
            # matches Flexas et al. (2002, p. 234) measured range 0.80-0.83 without tuning.
            @test ref.dark_pam.fv_fm isa FT
            @test FT(0.80) < ref.dark_pam.fv_fm < FT(0.83)
            @test ref.dark_psii.fv_fm > FT(0.88)

            tables = flexas_2002_curves(FT)
            fig2 = tables["flexas_mechanistic_fig2"]
            emp = tables["flexas_empirical_fig2"]
            inset = tables["flexas_inset_conductance"]
            ctrl = tables["flexas_controls"]

            for tbl in (fig2, emp, inset, ctrl), row in tbl
                @test all(v -> v isa FT || v isa Bool || v === nothing, values(row))
            end

            # 1. Dark basal normalization Fs(0) / Fo == 1 for PAM-2000 and dark gaps for actinic-only yields.
            @test fig2[1].dark
            @test fig2[1].g150_fs_over_fo == one(FT)
            @test fig2[1].g10_fs_over_fo == one(FT)
            @test emp[1].cotton_g150 === nothing
            @test emp[1].drought_g10 === nothing
            @test ctrl[1].dtt_unquenched === nothing

            # 2. Mechanistic Cytochrome + PAM-2000 reproduces all three Fig. 2 drought signatures:
            #    (a) monotonic decline of peak Fs/Fo across g = 150 -> 80 -> 62 -> 55 -> 39 -> 10,
            #    (b) monotonic leftward shift of peak PPFD,
            #    (c) high-light sub-Fo quenching (Fs/Fo < 1.0) under severe water stress.
            series_keys = (
                :g150_fs_over_fo,
                :g80_fs_over_fo,
                :g62_fs_over_fo,
                :g55_fs_over_fo,
                :g39_fs_over_fo,
                :g10_fs_over_fo,
            )
            peaks = map(series_keys) do k
                vals = [getproperty(r, k) for r in fig2]
                val, idx = findmax(vals)
                (ppfd = fig2[idx].incident_par, peak = val)
            end
            for i in 1:(length(peaks) - 1)
                @test peaks[i].peak > peaks[i + 1].peak
                @test peaks[i].ppfd >= peaks[i + 1].ppfd
            end
            @test peaks[1].ppfd > peaks[end].ppfd
            @test fig2[end].g150_fs_over_fo > one(FT)
            @test fig2[end].g39_fs_over_fo < one(FT)
            @test fig2[end].g10_fs_over_fo < one(FT)

            # 3. Single-variable KN(x) peak-invariance & Cotton high-light sign inversion:
            lit_emp = emp[2:end]
            cot_hi_peak = maximum(r.cotton_g150 for r in lit_emp)
            cot_lo_peak = maximum(r.cotton_g10 for r in lit_emp)
            drt_hi_peak = maximum(r.drought_g150 for r in lit_emp)
            drt_lo_peak = maximum(r.drought_g10 for r in lit_emp)
            @test abs(cot_hi_peak - cot_lo_peak) < FT(0.02)
            @test abs(drt_hi_peak - drt_lo_peak) < FT(0.02)
            # At high irradiance, empirical Cotton KN(x) predicts higher Fs/Fo under severe drought (wrong sign).
            @test emp[end].cotton_g10 > emp[end].cotton_g150

            # 4. Figure 2 Inset (PPFD = 750): positive curvilinear rise vs. g_sw for Cytochrome + PAM-2000,
            #    whereas Cotton KN(x) decreases with increasing g_sw.
            @test inset[1].pam2000_combined < one(FT)
            @test inset[end].pam2000_combined > FT(1.38)
            @test issorted([r.pam2000_combined for r in inset])
            @test issorted([r.pam2000_stomatal_only for r in inset])
            @test inset[1].empirical_cotton > inset[end].empirical_cotton
            @test count(r -> r.observed_flexas !== nothing, inset) == 6

            # 5. Fig. 3C high-CO2 / O2-depleted reversion & Fig. 5 DTT unquenched surge.
            @test ctrl[end].severe_ambient_pam < one(FT)
            @test ctrl[end].high_co2_reversion_pam > FT(1.30)
            @test ctrl[end].dtt_unquenched > FT(2.30)
        end
    end
end
