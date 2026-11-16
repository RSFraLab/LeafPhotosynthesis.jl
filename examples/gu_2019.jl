# Gu et al. (2019, New Phytologist 223:1179-1191) lake-model relations between PSII
# fluorescence, the open reaction-center fraction qL, NPQ and linear electron flow:
# their Vitis vinifera light response (Figs. 1c and 4c), the fluorescence-to-assimilation
# curves of Fig. 2a, and a comparison with the Johnson & Berry (2021) Cytochrome b6f
# fluorescence model on the Flexas et al. (2002) water-stress Figure 2.
# Run: julia --project=. examples/gu_2019.jl [output_directory]
using LeafPhotosynthesis
include("setup/gu_2019.jl")

# One irradiance of the Flexas et al. (2002) benchmark at the drought stage of `case`.
# The Cytochrome b6f photosynthesis of the benchmark supplies the electron flow, and the
# PAM-2000 observer of the benchmark (PSII plus PSI emission) the detector signal. Three
# fluorescence formulations are compared with the measured Fs/Fo:
# - Johnson & Berry: connected PSII units (KU2) with their open fraction q2 and feedback quenching;
# - the lake model for PSII with that same open fraction q2, through the same observer;
# - the lake relations between the detector's own PAM quantities (apparent Fv/Fm and
#   photochemical yield of the benchmark), with the qL(PAR) of Gu et al. Fig. 1c, which is
#   such a PAM quantity, at every drought stage.
function gu_flexas_point(ref, case, ppfd::FT, observed::FT) where {FT<:AbstractFloat}
    jb = evaluate_flexas_point(case, ppfd)
    inputs = flexas_photosystem_inputs(case, ppfd)
    r_jb = fluorescence(case.fl_cyto, inputs)
    pam = pam_fluorescence(case.pam_pam2000, case.fl_cyto, r_jb, case.alpha2_dark, case.alpha1_dark)
    r_model = fluorescence(ref.host_lake, inputs)
    all(is_valid, (r_jb, pam, r_model)) ||
        error("Electron flow outside a fluorescence model's domain at g = $(case.g_sw), PPFD = $ppfd")

    # PSII lake emission as the PAM-2000 observer reports it, with the PSI emission of the
    # Cytochrome leaf.
    host_fo_yield = lake_model_yields(ref.host_lake, one(FT), zero(FT)).fluorescence_yield
    lake_model = FT(detector_fs_over_fo(
        case, r_model.fluorescence_yield, host_fo_yield, r_jb.psi_fluorescence_yield))

    # Lake relations on the detector's PAM quantities. `excitation` only sets the scale of
    # the fluxes: yields and Fs/Fo do not depend on it.
    apparent = ref.apparent_lake
    excitation = inputs.psii_excitation
    apparent_flow = pam.apparent_photochemical_yield * excitation
    q_measured = measured_open_fraction(ref.settings, ppfd)
    r_measured = fluorescence(apparent, PSIIOpenFractionInputs(excitation, apparent_flow, q_measured))
    # The measured Fs/Fo as an emitted flux, and what Eq. 21 makes of it: the open fraction
    # that reproduces the observation at the apparent electron flow of the benchmark, and
    # the electron flow at the open fraction of Gu et al. Fig. 1c.
    fo_yield = lake_model_yields(apparent, one(FT), zero(FT)).fluorescence_yield
    flow_at(q) = electron_flow_from_fluorescence(apparent,
        ObservedFluorescenceInputs(observed * fo_yield * excitation, one(FT), one(FT), q))
    # The lake qL a PSII-only fluorometer would report for the connected-units model.
    ideal_pam = pam_fluorescence(case.pam_ideal, case.fl_cyto, r_jb, case.alpha2_dark, case.alpha1_dark)

    return (;
        conductance_mmol = case.g_sw,
        incident_par = ppfd,
        observed_fs_over_fo = observed,
        johnson_berry_pam2000 = jb.pam2000_fs_over_fo,
        lake_model_open_fraction = lake_model,
        # NaN where the transferred qL allows less photochemistry than the benchmark has.
        lake_measured_open_fraction = is_valid(r_measured) ? FT(fs_over_fo(apparent, r_measured)) : FT(NaN),
        johnson_berry_psii = jb.psii_fs_over_fo,
        lake_model_open_fraction_psii = FT(fs_over_fo(ref.host_lake, r_model)),
        psii_photochemical_yield = r_model.photochemical_yield,
        apparent_photochemical_yield = pam.apparent_photochemical_yield,
        model_open_fraction = inputs.psii_open_fraction,
        model_pam_lake_open_fraction = ideal_pam.ql,
        measured_open_fraction = q_measured,
        implied_open_fraction = FT(apparent_flow / flow_at(one(FT))),
        johnson_berry_psii_npq = r_jb.npq,
        lake_model_npq = r_model.npq,
        johnson_berry_apparent_npq = pam.apparent_npq,
        linear_electron_flow = inputs.linear_electron_flow,
        apparent_electron_flow = FT(apparent_flow),
        electron_flow_from_observed_fs = FT(flow_at(q_measured)),
    )
end

function gu_2019_curves(::Type{FT}) where {FT<:AbstractFloat}
    ref = setup_gu_2019(FT)
    s = ref.settings
    fig1c, fig4c = s["digitized_figure1c"], s["digitized_figure4c"]

    # 1. Figs. 1c and 4c: yields from the measured qL and NPQ (Eqs. 14 and 16), and the
    #    photochemical yield recovered from the fluorescence yield with qL alone (Eq. 21).
    vitis_rows = map(eachindex(fig1c["par"])) do i
        par, q, npq = FT(fig1c["par"][i]), FT(fig1c["open_fraction"][i]), FT(fig1c["npq"][i])
        y = lake_model_yields(ref.paper_lake, q, npq)
        recovered = electron_flow_from_fluorescence(
            ref.paper_lake,
            ObservedFluorescenceInputs(y.fluorescence_yield, one(FT), one(FT), q),
        )
        (;
            incident_par = par,
            measured_open_fraction = q,
            measured_npq = npq,
            photochemical_yield = y.photochemical_yield,
            fluorescence_yield = y.fluorescence_yield,
            fs_over_fo = fs_over_fo(ref.paper_lake, y),
            yield_recovered_from_fluorescence = recovered,
            paper_photochemical_yield = FT(fig4c["photochemical_yield"][i]),
            paper_fluorescence_yield = FT(fig4c["fluorescence_yield"][i]),
        )
    end

    # 2. Fig. 2a: electron flow from the observed fluorescence (Eq. 21), then gross
    #    assimilation through series diffusion (Eq. 24 at r_ch = 0) by the package's
    #    prescribed-conductance gas exchange.
    f2 = s["figure2"]
    sif_grid = collect(f2["sif_start"]:f2["sif_step"]:f2["sif_stop"])
    resistances = FT.(f2["boundary_stomatal_resistances"])
    sif_rows = map(sif_grid) do sif
        flow = electron_flow_from_fluorescence(
            ref.paper_lake,
            ObservedFluorescenceInputs(sif, f2["escape_probability"], f2["psii_fraction_of_sif"], f2["open_fraction"]),
        )
        gpp = map(resistances) do r_sb
            leaf = prescribed_electron_flow_exchange(ref.paper_lake, flow, r_sb, s)
            exchange = gas_exchange(leaf.photosynthesis, leaf.environment, leaf.state, leaf.diffusion)
            is_valid(exchange) || error("Invalid gas exchange at SIF = $sif, r_sb = $r_sb")
            exchange.net_assimilation + exchange.respiration
        end
        (;
            sif,
            linear_electron_flow = flow,
            gpp_rsb_025 = gpp[1],
            gpp_rsb_050 = gpp[2],
            gpp_rsb_075 = gpp[3],
            gpp_rsb_100 = gpp[4],
            gpp_rsb_125 = gpp[5],
        )
    end

    # 3. Flexas et al. (2002) Figure 2 at every digitized lit point, and its Inset: the six
    #    drought stages at one irradiance.
    fx = ref.flexas.settings
    x, dig = fx["experiment"], fx["digitized_flexas_2002"]
    conductances = FT.(x["fig2_conductances_mmol"])
    cases = [setup_flexas_2002(FT; g_sw_mmol = g, enable_nonstomatal = true) for g in conductances]
    flexas_rows = reduce(vcat, map(cases) do case
        prefix = "g" * string(round(Int, Float64(case.g_sw)))
        lit = [(p, y) for (p, y) in zip(FT.(dig[prefix * "_ppfd"]), FT.(dig[prefix * "_fs_over_fo"])) if p > zero(FT)]
        [gu_flexas_point(ref, case, ppfd, observed) for (ppfd, observed) in lit]
    end)
    inset = Dict(g => y for (g, y) in zip(FT.(dig["inset_conductances_mmol"]), FT.(dig["inset_fs_over_fo"])))
    inset_rows = [gu_flexas_point(ref, case, x["inset_ppfd"], inset[case.g_sw])
        for case in sort(cases; by = c -> c.g_sw)]

    # 4. The Fig. 1c light response (measured qL and NPQ, no photosynthesis model, no PSI)
    #    against the two least-stressed Flexas curves, linearly interpolated to the Fig. 1c
    #    irradiances inside each curve's range.
    identity_rows = map(vitis_rows) do row
        (;
            row.incident_par,
            lake_measured_inputs = row.fs_over_fo,
            observed_g150 = interpolate_observation(FT, dig, "g150", row.incident_par),
            observed_g80 = interpolate_observation(FT, dig, "g80", row.incident_par),
            outside_g150_range = row.incident_par > FT(dig["g150_ppfd"][end]),
        )
    end

    # 5. Root-mean-square error of Fs/Fo over the 49 lit Figure 2 points, per conductance,
    #    and over the six Inset points. The `_psii` formulations are PSII emission alone
    #    against the measured signal, which contains PSI.
    fields = (:johnson_berry_pam2000, :lake_model_open_fraction, :lake_measured_open_fraction,
        :johnson_berry_psii, :lake_model_open_fraction_psii)
    function score(rows, field)
        valid = [r for r in rows if !isnan(getproperty(r, field))]
        return flexas_rmse(FT, [getproperty(r, field) for r in valid], [r.observed_fs_over_fo for r in valid])
    end
    summary_rows = map(fields) do field
        by_g = [score([r for r in flexas_rows if r.conductance_mmol == g], field) for g in conductances]
        (;
            formulation = String(field),
            lit_fig2_rmse = score(flexas_rows, field),
            inset_rmse = score(inset_rows, field),
            rmse_g150 = by_g[1], rmse_g80 = by_g[2], rmse_g62 = by_g[3],
            rmse_g55 = by_g[4], rmse_g39 = by_g[5], rmse_g10 = by_g[6],
        )
    end

    return Dict(
        "gu_vitis_light_response" => vitis_rows,
        "gu_sif_assimilation" => sif_rows,
        "gu_flexas_fig2" => flexas_rows,
        "gu_flexas_inset" => inset_rows,
        "gu_flexas_identity" => identity_rows,
        "gu_flexas_summary" => collect(summary_rows),
    )
end

# Observed Fs/Fo of one digitized Flexas curve at `par`, held at the last point beyond it.
function interpolate_observation(::Type{FT}, dig, prefix, par) where {FT<:AbstractFloat}
    xs, ys = FT.(dig[prefix * "_ppfd"]), FT.(dig[prefix * "_fs_over_fo"])
    par >= xs[end] && return ys[end]
    i = findlast(<=(par), xs)
    return FT(ys[i] + (ys[i + 1] - ys[i]) * (par - xs[i]) / (xs[i + 1] - xs[i]))
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "gu_2019") : only(ARGS)
    for FT in (Float32, Float64)
        write_tables(joinpath(directory, string(FT)), gu_2019_curves(FT))
    end
    println("Gu et al. (2019) lake-model fluorescence tables written to ", directory)
end
