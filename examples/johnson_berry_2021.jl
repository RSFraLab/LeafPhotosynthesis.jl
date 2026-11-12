# Johnson & Berry (2021) forward benchmarks: PSII exciton connectivity (Case 9),
# analytical state transitions qT (Case 10), coordinated Cyt b6f + Rubisco capacity
# scaling (Case 12), and apparent PAM fluorometry with PSI detector leakage (Eqs. 38-42).
# Run: julia --project=. examples/johnson_berry_2021.jl [output_directory]
using LeafPhotosynthesis
include("setup/johnson_berry_2021.jl")

function johnson_berry_2021_curves(::Type{FT}) where {FT}
    base = setup_johnson_berry_2021(FT)
    x = base.settings["experiment"]
    absorptance = x["leaf_absorptance"]
    T, co2, O2 = x["temperature_K"], x["co2_Pa"], x["oxygen_Pa"]
    par_grid = x["incident_par_start"]:x["incident_par_step"]:x["incident_par_stop"]
    scales = FT.(x["capacity_scales"])
    scaled_cases = [setup_johnson_berry_2021(FT; capacity_scale = s) for s in scales]

    alpha2_dark = absorptance * base.beta_dark
    alpha1_dark = absorptance * (one(FT) - base.beta_dark)

    connectivity_rows = map(par_grid) do q_inc
        q_abs = absorptance * q_inc
        site = LeafEnvironment(T, q_abs, co2, O2)
        static_exc = exact_dual_excitation(base.beta_dark, q_abs)
        inp = photosystem_flux_inputs(base.photo, site, static_exc)
        r_pud = fluorescence(base.fl_puddle, inp)
        r_ku2 = fluorescence(base.fl_ku2, inp)
        r_ku10 = fluorescence(base.fl_ku10, inp)
        all(fluorescence_valid, (r_pud, r_ku2, r_ku10)) ||
            error("Invalid Case 9 fluorescence at incident_par = $q_inc")
        dark = r_ku2.status == FluorescenceDark
        (;
            incident_par = q_inc,
            absorbed_par = q_abs,
            puddle_psii_yield = dark ? nothing : r_pud.psii_fluorescence_yield,
            ku2_psii_yield = dark ? nothing : r_ku2.psii_fluorescence_yield,
            ku10_psii_yield = dark ? nothing : r_ku10.psii_fluorescence_yield,
            dark,
        )
    end

    state_transition_rows = map(par_grid) do q_inc
        q_abs = absorptance * q_inc
        site = LeafEnvironment(T, q_abs, co2, O2)
        static_exc = exact_dual_excitation(base.beta_dark, q_abs)
        qt_exc_ku2 = state_transition_excitation(base.rule, base.fl_ku2, base.photo, site)
        qt_exc_ku10 = state_transition_excitation(base.rule, base.fl_ku10, base.photo, site)
        beta_ku2 = state_transition_partition(base.rule, base.fl_ku2, base.photo, site)
        beta_ku10 = state_transition_partition(base.rule, base.fl_ku10, base.photo, site)

        r_stat_ku2 = fluorescence(base.fl_ku2, photosystem_flux_inputs(base.photo, site, static_exc))
        r_qt_ku2 = fluorescence(base.fl_ku2, photosystem_flux_inputs(base.photo, site, qt_exc_ku2))
        r_stat_ku10 = fluorescence(base.fl_ku10, photosystem_flux_inputs(base.photo, site, static_exc))
        r_qt_ku10 = fluorescence(base.fl_ku10, photosystem_flux_inputs(base.photo, site, qt_exc_ku10))
        all(fluorescence_valid, (r_stat_ku2, r_qt_ku2, r_stat_ku10, r_qt_ku10)) ||
            error("Invalid Case 10 fluorescence at incident_par = $q_inc")
        dark = r_qt_ku2.status == FluorescenceDark
        (;
            incident_par = q_inc,
            absorbed_par = q_abs,
            beta_ku2,
            beta_ku10,
            static_ku2_yield = dark ? nothing : r_stat_ku2.psii_fluorescence_yield,
            qt_ku2_yield = dark ? nothing : r_qt_ku2.psii_fluorescence_yield,
            static_ku10_yield = dark ? nothing : r_stat_ku10.psii_fluorescence_yield,
            qt_ku10_yield = dark ? nothing : r_qt_ku10.psii_fluorescence_yield,
            dark,
        )
    end

    capacity_rows = map(par_grid) do q_inc
        q_abs = absorptance * q_inc
        site = LeafEnvironment(T, q_abs, co2, O2)
        res = map(scaled_cases) do sc
            exc = state_transition_excitation(sc.rule, sc.fl_ku2, sc.photo, site)
            r = fluorescence(sc.fl_ku2, photosystem_flux_inputs(sc.photo, site, exc))
            fluorescence_valid(r) || error("Invalid Case 12 fluorescence at incident_par = $q_inc")
            r
        end
        dark = res[1].status == FluorescenceDark
        (;
            incident_par = q_inc,
            absorbed_par = q_abs,
            cap_05x_yield = dark ? nothing : res[1].psii_fluorescence_yield,
            cap_10x_yield = dark ? nothing : res[2].psii_fluorescence_yield,
            cap_15x_yield = dark ? nothing : res[3].psii_fluorescence_yield,
            cap_25x_yield = dark ? nothing : res[4].psii_fluorescence_yield,
            dark,
        )
    end

    pam_rows = map(par_grid) do q_inc
        q_abs = absorptance * q_inc
        site = LeafEnvironment(T, q_abs, co2, O2)
        static_exc = exact_dual_excitation(base.beta_dark, q_abs)
        r = fluorescence(base.fl_ku2, photosystem_flux_inputs(base.photo, site, static_exc))
        fluorescence_valid(r) || error("Invalid PAM benchmark fluorescence at incident_par = $q_inc")
        dark = r.status == FluorescenceDark
        p_ideal = pam_fluorescence(base.pam_ideal, base.fl_ku2, r, alpha2_dark, alpha1_dark)
        p_leak = pam_fluorescence(base.pam_leaked, base.fl_ku2, r, alpha2_dark, alpha1_dark)
        dark || (pam_fluorescence_valid(p_ideal) && pam_fluorescence_valid(p_leak)) ||
            error("Invalid PAM observation at incident_par = $q_inc")
        (;
            incident_par = q_inc,
            absorbed_par = q_abs,
            true_psii_photochemical_yield = dark ? nothing : r.psii_photochemical_yield,
            pam_apparent_photochemical_yield = dark ? nothing : p_leak.apparent_photochemical_yield,
            ideal_closed_lake_1_minus_ql = dark ? nothing : one(FT) - p_ideal.ql,
            pam_closed_lake_1_minus_ql = dark ? nothing : one(FT) - p_leak.ql,
            true_closed_psii_fraction = dark ? nothing : one(FT) - r.psii_open_fraction,
            pam_closed_puddle_1_minus_qp = dark ? nothing : one(FT) - p_leak.qp,
            ideal_closed_puddle_1_minus_qp = dark ? nothing : one(FT) - p_ideal.qp,
            pam_apparent_npq = dark ? nothing : p_leak.apparent_npq,
            true_internal_npq = dark ? nothing : r.npq,
            dark,
        )
    end

    return Dict(
        "jb_connectivity" => connectivity_rows,
        "jb_state_transitions" => state_transition_rows,
        "jb_capacity" => capacity_rows,
        "jb_pam" => pam_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "johnson_berry_2021") : only(ARGS)
    for FT in (Float32, Float64)
        write_tables(joinpath(directory, string(FT)), johnson_berry_2021_curves(FT))
    end
    println("Johnson & Berry (2021) benchmark curves written to ", directory)
end
