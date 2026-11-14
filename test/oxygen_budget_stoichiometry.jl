# The coupled oxygen budget is computed from carboxylation and oxygenation, so with a
# respiratory quotient of one net O2 release equals net CO2 uptake for every valid electron
# stoichiometry, not only for e_c = 4, a = 2.

@testset "Oxygen budget for any electron stoichiometry" begin
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
        leaf, _ = core_leaf(FT, p, settings)
        base = leaf.photosynthesis
        for stoichiometry in (C3Stoichiometry(FT(4), FT(2)), C3Stoichiometry(FT(4.5), FT(10.5) / FT(4.5)),
                C3Stoichiometry(FT(5), FT(1)))
            model = C3Photosynthesis(base.temperature, base.electron_transport, stoichiometry,
                base.limitation, base.assimilation_sign)
            for co2 in FT.((2, 5, 20, 40, 100)), par in FT.((0, 200, 1600))
                site = LeafEnvironment(FT(298.15), par, co2, FT(21000))
                ox = oxygen_exchange(leaf.oxygen, model, site)
                carbon = photosynthesize(model, site)
                turnover = biochemical_turnover(model, site)
                scale = turnover.carboxylation + turnover.photorespiratory_co2 + turnover.respiration
                @test isapprox(ox.net_release, carbon.net_assimilation; atol = 16eps(FT) * scale)
                @test isapprox(ox.gross_production, ox.carboxylation + ox.rubisco_oxygenation;
                    rtol = 4eps(FT), atol = 4eps(FT) * scale)
            end
        end
    end
end
