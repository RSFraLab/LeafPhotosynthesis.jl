# Stress-response comparison: empirical PSII rate competition (van der Tol et al. 2014)
# vs. mechanistic Cytochrome b6f + coupled PSII/PSI fluorescence (Johnson & Berry 2021).
# Run: julia --project=. examples/cytochrome_stress.jl [output_directory]
include("setup/cytochrome_stress.jl")

function cytochrome_stress_curves(::Type{FT}) where {FT}
    base = setup_cytochrome_stress(FT)
    x = base.settings["experiment"]
    T, O2 = x["temperature_K"], x["oxygen_Pa"]

    # 1. Light response at fixed biochemical-site CO2: biphasic peak & PSI baseline.
    light_rows = [
        evaluate_fluorescence_comparison(base, LeafEnvironment(T, par, x["light_sweep_co2_Pa"], O2))
        for par in x["par_start"]:x["par_step"]:x["par_stop"]
    ]

    # 2. Stomatal / drought CO2 drawdown at high (PAR = 1200) and moderate (PAR = 250) light.
    drought_rows = map(x["co2_drawdown_start"]:x["co2_drawdown_step"]:x["co2_drawdown_stop"]) do co2
        hi = evaluate_fluorescence_comparison(base, LeafEnvironment(T, x["high_light_par"], co2, O2))
        mod = evaluate_fluorescence_comparison(base, LeafEnvironment(T, x["moderate_light_par"], co2, O2))
        (;
            biochemical_co2 = co2,
            high_jb_psii_yield = hi.jb_psii_yield,
            high_jb_total_yield = hi.jb_total_yield,
            high_jb_npq = hi.jb_npq,
            high_vdt_cotton_yield = hi.vdt_cotton_yield,
            high_vdt_drought_yield = hi.vdt_drought_yield,
            mod_jb_psii_yield = mod.jb_psii_yield,
            mod_jb_total_yield = mod.jb_total_yield,
            mod_jb_eta = mod.jb_psi_to_psii_ratio,
            mod_vdt_cotton_yield = mod.vdt_cotton_yield,
            mod_vdt_drought_yield = mod.vdt_drought_yield,
        )
    end

    # 3. Coupled Medlyn stomatal equilibrium: Rubisco sink loss (Vcmax) vs. Cyt-b6f / Jmax loss.
    gas_env = GasExchangeEnvironment(
        T,
        x["capacity_stress_par"],
        x["ambient_co2_Pa"],
        O2,
        x["pressure_Pa"],
        x["relative_humidity"],
        x["vpd_Pa"],
    )
    capacity_rows = map(x["capacity_fraction_start"]:x["capacity_fraction_step"]:x["capacity_fraction_stop"]) do factor
        rubisco_case = setup_cytochrome_stress(FT; vcmax_factor = factor)
        thylakoid_case = setup_cytochrome_stress(FT; transport_factor = factor)

        sol_jb_vc = equilibrium_gas_exchange(rubisco_case.cyto_leaf, gas_env, rubisco_case.diffusion, rubisco_case.solver)
        sol_jb_vq = equilibrium_gas_exchange(thylakoid_case.cyto_leaf, gas_env, thylakoid_case.diffusion, thylakoid_case.solver)
        sol_vdt_vc = equilibrium_gas_exchange(rubisco_case.fvcb_leaf, gas_env, rubisco_case.diffusion, rubisco_case.solver)
        sol_vdt_jm = equilibrium_gas_exchange(thylakoid_case.fvcb_leaf, gas_env, thylakoid_case.diffusion, thylakoid_case.solver)
        all(s -> s.converged, (sol_jb_vc, sol_jb_vq, sol_vdt_vc, sol_vdt_jm)) ||
            error("Equilibrium stress solve failed at capacity_fraction = $factor")

        fl_jb_vc = fluorescence(rubisco_case.jb_connected,
            photosystem_flux_inputs(rubisco_case.cyto_leaf.photosynthesis,
                LeafEnvironment(T, gas_env.absorbed_par, sol_jb_vc.biochemical_co2, O2)))
        fl_jb_vq = fluorescence(thylakoid_case.jb_connected,
            photosystem_flux_inputs(thylakoid_case.cyto_leaf.photosynthesis,
                LeafEnvironment(T, gas_env.absorbed_par, sol_jb_vq.biochemical_co2, O2)))
        fl_vdt_vc = fluorescence(rubisco_case.cotton_fl,
            psii_flux_inputs(rubisco_case.fvcb_leaf.photosynthesis,
                LeafEnvironment(T, gas_env.absorbed_par, sol_vdt_vc.biochemical_co2, O2)))
        fl_vdt_jm = fluorescence(thylakoid_case.cotton_fl,
            psii_flux_inputs(thylakoid_case.fvcb_leaf.photosynthesis,
                LeafEnvironment(T, gas_env.absorbed_par, sol_vdt_jm.biochemical_co2, O2)))
        all(fluorescence_valid, (fl_jb_vc, fl_jb_vq, fl_vdt_vc, fl_vdt_jm)) ||
            error("Fluorescence invalid at capacity_fraction = $factor")

        (;
            capacity_fraction = factor,
            jb_rubisco_psii_yield = fl_jb_vc.psii_fluorescence_yield,
            jb_b6f_psii_yield = fl_jb_vq.psii_fluorescence_yield,
            jb_rubisco_total_yield = fl_jb_vc.total_fluorescence_yield,
            jb_b6f_total_yield = fl_jb_vq.total_fluorescence_yield,
            jb_rubisco_npq = fl_jb_vc.npq,
            jb_b6f_npq = fl_jb_vq.npq,
            jb_rubisco_q2 = fl_jb_vc.psii_open_fraction,
            jb_b6f_q2 = fl_jb_vq.psii_open_fraction,
            vdt_rubisco_yield = fl_vdt_vc.fluorescence_yield,
            vdt_jmax_yield = fl_vdt_jm.fluorescence_yield,
        )
    end

    return Dict(
        "cytochrome_light" => light_rows,
        "cytochrome_drought" => drought_rows,
        "cytochrome_capacity" => capacity_rows,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "cytochrome_stress") : only(ARGS)
    for FT in (Float32, Float64)
        write_tables(joinpath(directory, string(FT)), cytochrome_stress_curves(FT))
    end
    println("Cytochrome stress comparison curves written to ", directory)
end
