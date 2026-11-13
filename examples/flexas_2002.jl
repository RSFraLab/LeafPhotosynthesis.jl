# Flexas et al. (2002, Physiologia Plantarum 114:231-240) water-stress benchmark:
# steady-state normalized chlorophyll fluorescence (Fs / Fo) light-response curves across
# six drought stomatal conductances (Fig. 2), Fs / Fo vs. conductance at PPFD = 750 (Fig. 2 Inset),
# high-CO2 / low-O2 reversion (Fig. 3C), and DTT xanthophyll inhibition (Fig. 5).
# Run: julia --project=. examples/flexas_2002.jl [output_directory]
using LeafPhotosynthesis
include("setup/flexas_2002.jl")

function flexas_2002_curves(::Type{FT}) where {FT<:AbstractFloat}
    ref = setup_flexas_2002(FT)
    x = ref.settings["experiment"]
    dig = ref.settings["digitized_flexas_2002"]

    ppfd_grid = collect(x["ppfd_start"]:x["ppfd_step"]:x["ppfd_stop"])
    g_levels = FT.(x["fig2_conductances_mmol"])
    cases_combined = [setup_flexas_2002(FT; g_sw_mmol = g, enable_nonstomatal = true) for g in g_levels]
    cases_stomatal_only = [setup_flexas_2002(FT; g_sw_mmol = g, enable_nonstomatal = false) for g in g_levels]

    # Full conductance x irradiance evaluation grid across all fluorescence formulations
    # and internal thylakoid/stromal diagnostics.
    grid_by_g = map(eachindex(g_levels)) do idx
        g_mmol = g_levels[idx]
        c_comb = cases_combined[idx]
        c_stom = cases_stomatal_only[idx]
        map(ppfd_grid) do ppfd
            p_comb = evaluate_flexas_point(c_comb, ppfd)
            p_stom = evaluate_flexas_point(c_stom, ppfd)
            (;
                conductance_mmol = g_mmol,
                incident_par = ppfd,
                pam2000_combined = p_comb.pam2000_fs_over_fo,
                pam2000_stomatal_only = p_stom.pam2000_fs_over_fo,
                pure_psii_combined = p_comb.psii_fs_over_fo,
                qt_pam2000_combined = p_comb.qt_pam2000_fs_over_fo,
                cotton_combined = p_comb.cotton_fs_over_fo,
                drought_combined = p_comb.drought_fs_over_fo,
                cotton_stomatal_only = p_stom.cotton_fs_over_fo,
                drought_stomatal_only = p_stom.drought_fs_over_fo,
                dtt_unquenched = p_comb.unquenched_fs_over_fo,
                psii_open_fraction = p_comb.psii_open_fraction,
                psii_photochemical_yield = p_comb.psii_photochemical_yield,
                pam_apparent_photochemical_yield = p_comb.pam_apparent_photochemical_yield,
                psii_quenching_rate = p_comb.psii_quenching_rate,
                pam2000_npq = p_comb.pam2000_npq,
                net_assimilation = p_comb.net_assimilation,
                biochemical_co2 = p_comb.biochemical_co2,
                dark = p_comb.dark,
            )
        end
    end
    full_grid_rows = reduce(vcat, grid_by_g)

    # 1. Figure 2 main panel: Fs / Fo vs. PPFD across all 6 water-stress stomatal conductances.
    mechanistic_rows = map(eachindex(ppfd_grid)) do i
        (;
            incident_par = ppfd_grid[i],
            g150_fs_over_fo = grid_by_g[1][i].pam2000_combined,
            g80_fs_over_fo = grid_by_g[2][i].pam2000_combined,
            g62_fs_over_fo = grid_by_g[3][i].pam2000_combined,
            g55_fs_over_fo = grid_by_g[4][i].pam2000_combined,
            g39_fs_over_fo = grid_by_g[5][i].pam2000_combined,
            g10_fs_over_fo = grid_by_g[6][i].pam2000_combined,
            dark = grid_by_g[1][i].dark,
        )
    end

    # 2. Empirical van der Tol (2014) Cotton & Drought KN(x) comparison on well-watered (g = 150)
    #    vs. severe drought (g = 10), using pure stomatal closure to isolate KN(x) behavior.
    empirical_rows = map(eachindex(ppfd_grid)) do i
        hi_g = grid_by_g[1][i]
        lo_g = grid_by_g[6][i]
        (;
            incident_par = ppfd_grid[i],
            cotton_g150 = hi_g.cotton_stomatal_only,
            cotton_g10 = lo_g.cotton_stomatal_only,
            drought_g150 = hi_g.drought_stomatal_only,
            drought_g10 = lo_g.drought_stomatal_only,
            dark = hi_g.dark,
        )
    end

    # 3. Figure 2 Inset: Fs / Fo vs. stomatal conductance g at PPFD = 750 micromol m^-2 s^-1.
    obs_g = FT.(dig["inset_conductances_mmol"])
    obs_fsfo = FT.(dig["inset_fs_over_fo"])
    obs_lookup = Dict(g => v for (g, v) in zip(obs_g, obs_fsfo))
    g_sweep = sort!(unique!(vcat(
        collect(x["inset_g_start_mmol"]:x["inset_g_step_mmol"]:x["inset_g_stop_mmol"]),
        obs_g,
    )))
    inset_ppfd = x["inset_ppfd"]
    inset_rows = map(g_sweep) do g_mmol
        c_comb = setup_flexas_2002(FT; g_sw_mmol = g_mmol, enable_nonstomatal = true)
        c_stom = setup_flexas_2002(FT; g_sw_mmol = g_mmol, enable_nonstomatal = false)
        p_comb = evaluate_flexas_point(c_comb, inset_ppfd)
        p_stom = evaluate_flexas_point(c_stom, inset_ppfd)
        (;
            conductance_mmol = g_mmol,
            pam2000_combined = p_comb.pam2000_fs_over_fo,
            pam2000_stomatal_only = p_stom.pam2000_fs_over_fo,
            pure_psii_combined = p_comb.psii_fs_over_fo,
            qt_pam2000_combined = p_comb.qt_pam2000_fs_over_fo,
            cotton_combined = p_comb.cotton_fs_over_fo,
            drought_combined = p_comb.drought_fs_over_fo,
            empirical_cotton = p_stom.cotton_fs_over_fo,
            empirical_drought = p_stom.drought_fs_over_fo,
            observed_flexas = get(obs_lookup, g_mmol, nothing),
            net_assimilation = p_comb.net_assimilation,
            pam2000_npq = p_comb.pam2000_npq,
        )
    end

    # 4. Mechanistic controls at severe water stress (g = 39 mmol m^-2 s^-1):
    #    ambient air vs. Fig. 3C high-CO2 / O2-depleted reversion, qT state transitions,
    #    and Fig. 5 DTT xanthophyll-inhibited unquenched limit (KN = 0).
    severe_case = setup_flexas_2002(FT; g_sw_mmol = x["controls_conductance_mmol"], enable_nonstomatal = true)
    controls_grid = collect(x["ppfd_start"]:x["ppfd_step"]:x["controls_ppfd_stop"])
    controls_rows = map(controls_grid) do ppfd
        amb = evaluate_flexas_point(severe_case, ppfd)
        rev = evaluate_flexas_point(
            severe_case,
            ppfd;
            co2_Pa = x["high_co2_Pa"],
            oxygen_Pa = x["low_oxygen_Pa"],
        )
        (;
            incident_par = ppfd,
            severe_ambient_pam = amb.pam2000_fs_over_fo,
            high_co2_reversion_pam = rev.pam2000_fs_over_fo,
            severe_qt_pam = amb.qt_pam2000_fs_over_fo,
            dtt_unquenched = amb.unquenched_fs_over_fo,
            dark = amb.dark,
        )
    end

    # 5. Subpixel-digitized Figure 2 full light-response trajectories (all 6 curves)
    #    evaluated directly at each exact observed PPFD in Julia (zero browser interpolation).
    digitized_by_g = map(eachindex(g_levels)) do idx
        g_mmol = g_levels[idx]
        prefix = "g" * string(round(Int, Float64(g_mmol)))
        obs_ppfd = FT.(dig[prefix * "_ppfd"])
        obs_y = FT.(dig[prefix * "_fs_over_fo"])
        c_comb = cases_combined[idx]
        c_stom = cases_stomatal_only[idx]
        map(zip(obs_ppfd, obs_y)) do (ppfd, y_obs)
            p_comb = evaluate_flexas_point(c_comb, ppfd)
            p_stom = evaluate_flexas_point(c_stom, ppfd)
            (;
                conductance_mmol = g_mmol,
                incident_par = ppfd,
                observed_fs_over_fo = y_obs,
                pam2000_combined = p_comb.pam2000_fs_over_fo,
                pam2000_stomatal_only = p_stom.pam2000_fs_over_fo,
                pure_psii_combined = p_comb.psii_fs_over_fo,
                qt_pam2000_combined = p_comb.qt_pam2000_fs_over_fo,
                cotton_combined = p_comb.cotton_fs_over_fo,
                drought_combined = p_comb.drought_fs_over_fo,
                dark = p_comb.dark,
            )
        end
    end
    digitized_rows = reduce(vcat, digitized_by_g)

    # 6. Per-conductance peak & RMSE summary (ascending g_sw: 10 -> 150 mmol m^-2 s^-1).
    asc_indices = sortperm(g_levels)
    conductance_summary_rows = map(asc_indices) do idx
        g_mmol = g_levels[idx]
        c_comb = cases_combined[idx]
        g_rows = grid_by_g[idx]
        d_lit = [r for r in digitized_by_g[idx] if !r.dark]
        obs_lit = [r.observed_fs_over_fo for r in d_lit]
        obs_peak_val, obs_peak_i = findmax(obs_lit)
        pk_comb = flexas_lit_peak(FT, g_rows, :pam2000_combined)
        pk_stom = flexas_lit_peak(FT, g_rows, :pam2000_stomatal_only)
        pk_psii = flexas_lit_peak(FT, g_rows, :pure_psii_combined)
        pk_qt = flexas_lit_peak(FT, g_rows, :qt_pam2000_combined)
        pk_cot = flexas_lit_peak(FT, g_rows, :cotton_combined)
        pk_drt = flexas_lit_peak(FT, g_rows, :drought_combined)
        (;
            conductance_mmol = g_mmol,
            vcmax_fraction = c_comb.f_vc,
            vqmax_fraction = c_comb.f_vq,
            observed_peak_fs_over_fo = obs_peak_val,
            observed_peak_ppfd = d_lit[obs_peak_i].incident_par,
            observed_inset_fs_over_fo = obs_lookup[g_mmol],
            pam2000_combined_peak = pk_comb.peak,
            pam2000_combined_peak_ppfd = pk_comb.ppfd,
            pam2000_stomatal_peak = pk_stom.peak,
            pam2000_stomatal_peak_ppfd = pk_stom.ppfd,
            pure_psii_peak = pk_psii.peak,
            pure_psii_peak_ppfd = pk_psii.ppfd,
            qt_pam2000_peak = pk_qt.peak,
            qt_pam2000_peak_ppfd = pk_qt.ppfd,
            cotton_peak = pk_cot.peak,
            cotton_peak_ppfd = pk_cot.ppfd,
            drought_peak = pk_drt.peak,
            drought_peak_ppfd = pk_drt.ppfd,
            pam2000_combined_rmse = flexas_rmse(FT, [r.pam2000_combined for r in d_lit], obs_lit),
            pam2000_stomatal_rmse = flexas_rmse(FT, [r.pam2000_stomatal_only for r in d_lit], obs_lit),
            pure_psii_rmse = flexas_rmse(FT, [r.pure_psii_combined for r in d_lit], obs_lit),
            qt_pam2000_rmse = flexas_rmse(FT, [r.qt_pam2000_combined for r in d_lit], obs_lit),
            cotton_rmse = flexas_rmse(FT, [r.cotton_combined for r in d_lit], obs_lit),
            drought_rmse = flexas_rmse(FT, [r.drought_combined for r in d_lit], obs_lit),
        )
    end

    # 7. Overall formulation scoreboard computed in Julia across all lit Fig. 2 points (n = 49)
    #    and all Fig. 2 Inset points (n = 6).
    lit_digitized = [r for r in digitized_rows if !r.dark]
    inset_obs_rows = [r for r in inset_rows if r.observed_flexas !== nothing]
    paper_stop_idx = findlast(<=(x["paper_domain_ppfd_stop"]), ppfd_grid)
    formulation_specs = [
        (:pam2000_combined, ref.dark_pam.fv_fm),
        (:pam2000_stomatal_only, ref.dark_pam.fv_fm),
        (:pure_psii_combined, ref.dark_psii.fv_fm),
        (:qt_pam2000_combined, ref.dark_pam.fv_fm),
        (:cotton_combined, ref.empirical_fv_fm),
        (:drought_combined, ref.empirical_fv_fm),
    ]
    formulation_summary_rows = map(formulation_specs) do (field, fv_fm)
        pk_hi = flexas_lit_peak(FT, grid_by_g[1], field)
        pk_lo = flexas_lit_peak(FT, grid_by_g[6], field)
        (;
            formulation = String(field),
            dark_fv_fm = fv_fm,
            lit_fig2_rmse = flexas_rmse(
                FT,
                [getproperty(r, field) for r in lit_digitized],
                [r.observed_fs_over_fo for r in lit_digitized],
            ),
            inset_rmse = flexas_rmse(
                FT,
                [getproperty(r, field) for r in inset_obs_rows],
                [r.observed_flexas for r in inset_obs_rows],
            ),
            peak_g150 = pk_hi.peak,
            peak_g10 = pk_lo.peak,
            peak_suppression = FT(pk_hi.peak - pk_lo.peak),
            peak_ppfd_g150 = pk_hi.ppfd,
            peak_ppfd_g10 = pk_lo.ppfd,
            inset_g150 = FT(getproperty(inset_obs_rows[end], field)),
            inset_g10 = FT(getproperty(inset_obs_rows[1], field)),
            high_light_g10 = FT(getproperty(grid_by_g[6][paper_stop_idx], field)),
        )
    end

    return Dict(
        "flexas_mechanistic_fig2" => mechanistic_rows,
        "flexas_empirical_fig2" => empirical_rows,
        "flexas_inset_conductance" => inset_rows,
        "flexas_controls" => controls_rows,
        "flexas_full_grid" => full_grid_rows,
        "flexas_digitized_fig2" => digitized_rows,
        "flexas_conductance_summary" => conductance_summary_rows,
        "flexas_formulation_summary" => formulation_summary_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "flexas_2002") : only(ARGS)
    for FT in (Float32, Float64)
        write_tables(joinpath(directory, string(FT)), flexas_2002_curves(FT))
    end
    println("Flexas et al. (2002) water-stress fluorescence curves written to ", directory)
end
