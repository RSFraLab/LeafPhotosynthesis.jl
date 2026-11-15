# Owner's answer 9 (2026-11-14): the stand-alone `oxygen_exchange(::CoreOxygenExchange,
# ::OxygenFluxInputs)` takes carboxylation and oxygenation, as the coupled form does, so both
# give the same budget, bit for bit, for the same turnover at any electron stoichiometry.

@testset "Stand-alone and coupled oxygen budgets agree bit for bit" begin
    @test !hasfield(CoreOxygenExchange{Float64}, :electrons_per_oxygen)
    @test fieldnames(OxygenFluxInputs) == (:carboxylation, :rubisco_oxygenation, :respiration)
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
        leaf, _ = core_leaf(FT, p, settings)
        ox = leaf.oxygen
        c3 = leaf.photosynthesis
        cyto = setup_johnson_berry_2021(FT).photo
        models = (
            "C3, classical" => c3,
            "C3, e_c = 4.5" => C3Photosynthesis(c3.temperature, c3.electron_transport,
                C3Stoichiometry(FT(4.5), FT(10.5) / FT(4.5)), c3.limitation, c3.assimilation_sign),
            "Cytochrome, classical" => cyto,
            "Cytochrome, e_c = 4.5" => CytochromePhotosynthesis(cyto.temperature, cyto.electron_transport,
                CytochromeStoichiometry(FT(4.5), FT(10.5), FT(3.5), FT(8)), cyto.limitation,
                cyto.assimilation_sign),
        )
        for (label, model) in models, co2 in FT.((2, 5, 20, 40, 100)), par in FT.((0, 200, 1600))
            site = LeafEnvironment(FT(298.15), par, co2, FT(21000))
            coupled = oxygen_exchange(ox, model, site)
            turnover = biochemical_turnover(model, site)
            oxygenation = turnover.photorespiratory_co2 / ox.co2_per_oxygenation
            alone = @inferred oxygen_exchange(ox, OxygenFluxInputs(turnover.carboxylation, oxygenation,
                turnover.respiration))
            budget = (; coupled.gross_production, coupled.photorespiratory_uptake,
                coupled.respiratory_uptake, coupled.net_release)
            @test alone === budget
            @test oxygen_exchange(ox, OxygenFluxInputs(coupled.carboxylation, coupled.rubisco_oxygenation,
                turnover.respiration)) === budget
            @test all(isfinite, values(alone))
            # The stand-alone form keeps net O₂ = net CO₂ (RQ = 1) at every stoichiometry.
            carbon = photosynthesize(model, site)
            scale = turnover.carboxylation + turnover.photorespiratory_co2 + turnover.respiration
            @test isapprox(alone.net_release, carbon.net_assimilation; atol = 16eps(FT) * scale)
        end
    end
end
