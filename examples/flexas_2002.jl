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

    # 1. Figure 2 main panel: Fs / Fo vs. PPFD across all 6 water-stress stomatal conductances.
    mechanistic_rows = map(ppfd_grid) do ppfd
        pts = [evaluate_flexas_point(c, ppfd) for c in cases_combined]
        (;
            incident_par = ppfd,
            g150_fs_over_fo = pts[1].pam2000_fs_over_fo,
            g80_fs_over_fo = pts[2].pam2000_fs_over_fo,
            g62_fs_over_fo = pts[3].pam2000_fs_over_fo,
            g55_fs_over_fo = pts[4].pam2000_fs_over_fo,
            g39_fs_over_fo = pts[5].pam2000_fs_over_fo,
            g10_fs_over_fo = pts[6].pam2000_fs_over_fo,
            dark = pts[1].dark,
        )
    end

    # 2. Empirical van der Tol (2014) Cotton & Drought KN(x) comparison on well-watered (g = 150)
    #    vs. severe drought (g = 10), using pure stomatal closure to isolate KN(x) behavior.
    empirical_rows = map(ppfd_grid) do ppfd
        hi_g = evaluate_flexas_point(cases_stomatal_only[1], ppfd)
        lo_g = evaluate_flexas_point(cases_stomatal_only[6], ppfd)
        (;
            incident_par = ppfd,
            cotton_g150 = hi_g.cotton_fs_over_fo,
            cotton_g10 = lo_g.cotton_fs_over_fo,
            drought_g150 = hi_g.drought_fs_over_fo,
            drought_g10 = lo_g.drought_fs_over_fo,
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

    return Dict(
        "flexas_mechanistic_fig2" => mechanistic_rows,
        "flexas_empirical_fig2" => empirical_rows,
        "flexas_inset_conductance" => inset_rows,
        "flexas_controls" => controls_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "flexas_2002") : only(ARGS)
    for FT in (Float32, Float64)
        write_tables(joinpath(directory, string(FT)), flexas_2002_curves(FT))
    end
    println("Flexas et al. (2002) water-stress fluorescence curves written to ", directory)
end
