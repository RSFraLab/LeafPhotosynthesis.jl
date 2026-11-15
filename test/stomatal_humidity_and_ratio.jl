# Owner's answers 7 and 8 (2026-11-14): each stomatal law checks only the humidity driver it
# reads, and the H₂O/CO₂ diffusivity ratio is stored once, in `SeriesDiffusion`, and used by
# both the Medlyn law and the diffusion path.

function humidity_case(::Type{FT}) where {FT}
    p = read_bonan_parameters(FT, parameter_path)
    settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
    leaf, diffusion = core_leaf(FT, p, settings)
    env = core_environment(FT, p, settings)
    stress = setup_cytochrome_stress(FT)
    return (; leaf, diffusion, env, solver = core_solver(FT, settings), stress)
end

with_humidity(env, rh, vpd) = GasExchangeEnvironment(env.temperature, env.absorbed_par,
    env.ambient_co2, env.oxygen, env.pressure, oftype(env.pressure, rh), oftype(env.pressure, vpd))

@testset "Each stomatal law checks only the humidity driver it reads" begin
    for FT in (Float32, Float64)
        c = humidity_case(FT)
        # Ball–Berry reads relative humidity and ignores the VPD.
        bb_leaf = c.leaf                                   # PrognosticStomata(BallBerry)
        bb_steady = SteadyStateStomata(bb_leaf.stomata.target)
        reference = with_humidity(c.env, c.env.relative_humidity, c.env.vpd)
        state = StomatalState(FT(0.2))
        t_ref = leaf_tendencies(bb_leaf, reference, state, c.diffusion)
        e_ref = equilibrium_gas_exchange(bb_leaf.photosynthesis, bb_steady, reference, c.diffusion, c.solver)
        @test is_valid(e_ref) && isfinite(t_ref.conductance_target)
        for vpd in (NaN, Inf, -Inf)
            env = with_humidity(c.env, c.env.relative_humidity, vpd)
            t = leaf_tendencies(bb_leaf, env, state, c.diffusion)
            @test t.conductance_target === t_ref.conductance_target
            @test t.conductance_tendency === t_ref.conductance_tendency
            e = equilibrium_gas_exchange(bb_leaf.photosynthesis, bb_steady, env, c.diffusion, c.solver)
            @test e.status == EquilibriumConverged && e === e_ref
        end
        for rh in (NaN, -0.1, 1.5)
            env = with_humidity(c.env, rh, c.env.vpd)
            @test isnan(leaf_tendencies(bb_leaf, env, state, c.diffusion).conductance_target)
            @test equilibrium_gas_exchange(bb_leaf.photosynthesis, bb_steady, env, c.diffusion,
                c.solver).status == EquilibriumInvalidInput
        end
        # Medlyn reads the VPD and ignores relative humidity.
        cyto = c.stress.cyto_leaf                          # SteadyStateStomata(Medlyn)
        md_prognostic = LeafModel(; photosynthesis = cyto.photosynthesis,
            stomata = PrognosticStomata(cyto.stomata.target, FT(300), FT(600)), oxygen = cyto.oxygen)
        d = c.stress.diffusion
        m_ref = equilibrium_gas_exchange(cyto, reference, d, c.solver)
        tm_ref = leaf_tendencies(md_prognostic, reference, state, d)
        @test is_valid(m_ref) && isfinite(tm_ref.conductance_target)
        for rh in (NaN, -0.1, 1.5, Inf)
            env = with_humidity(c.env, rh, c.env.vpd)
            m = equilibrium_gas_exchange(cyto, env, d, c.solver)
            @test m.status == EquilibriumConverged && m === m_ref
            tm = leaf_tendencies(md_prognostic, env, state, d)
            @test tm.conductance_target === tm_ref.conductance_target
            @test tm.conductance_tendency === tm_ref.conductance_tendency
        end
        for vpd in (NaN, Inf)
            env = with_humidity(c.env, c.env.relative_humidity, vpd)
            @test equilibrium_gas_exchange(cyto, env, d, c.solver).status == EquilibriumInvalidInput
            @test isnan(leaf_tendencies(md_prognostic, env, state, d).conductance_target)
        end
        # A model that supplies only the driver of its law: the other one NaN.
        bb_only = StomatalDrivers{FT}(10.0, 30.0, 101325.0, 0.6, NaN)
        md_only = StomatalDrivers{FT}(10.0, 30.0, 101325.0, NaN, 1200.0)
        @test isfinite(target_conductance(bb_leaf.stomata.target, bb_only, c.diffusion))
        @test isfinite(target_conductance(cyto.stomata.target, md_only, d))
        @test isnan(target_conductance(bb_leaf.stomata.target, md_only, c.diffusion))
        @test isnan(target_conductance(cyto.stomata.target, bb_only, d))
    end
end

@testset "One H₂O/CO₂ diffusivity ratio, read from the diffusion" begin
    @test !hasfield(Medlyn{Float64}, :diffusivity_ratio)
    @test hasfield(SeriesDiffusion{Float64}, :diffusivity_ratio)
    for FT in (Float32, Float64)
        c = humidity_case(FT)
        cyto = c.stress.cyto_leaf
        target = cyto.stomata.target
        d16 = c.stress.diffusion
        d14 = SeriesDiffusion(d16.boundary_co2, d16.mesophyll_co2, FT(1.4), d16.mole_to_micromole)
        drivers = StomatalDrivers{FT}(10.0, 30.0, 101325.0, 0.6, 1200.0)
        # The Medlyn target scales with the diffusion's ratio: (g - g0) ∝ r.
        g16 = target_conductance(target, drivers, d16)
        g14 = target_conductance(target, drivers, d14)
        @test isapprox((g14 - target.minimum_conductance) / (g16 - target.minimum_conductance),
            FT(1.4) / FT(1.6); rtol = 8eps(FT))
        # Ball–Berry does not use the ratio.
        bb = c.leaf.stomata.target
        @test target_conductance(bb, drivers, d14) === target_conductance(bb, drivers, d16)
        # An invalid ratio is invalid for the law that reads it.
        dnan = SeriesDiffusion(d16.boundary_co2, d16.mesophyll_co2, FT(NaN), d16.mole_to_micromole)
        @test isnan(target_conductance(target, drivers, dnan))
        @test isfinite(target_conductance(bb, drivers, dnan))
        # The equilibrium uses the one ratio in the target and in the diffusion path: the
        # solved conductance is the target at the solved drivers with that same diffusion.
        for d in (d16, d14)
            sol = equilibrium_gas_exchange(cyto, c.env, d, c.solver)
            @test is_valid(sol)
            x = sol.exchange
            at_solution = StomatalDrivers(x.net_assimilation, x.surface_co2, c.env.pressure,
                c.env.relative_humidity, c.env.vpd)
            @test isapprox(sol.stomatal_conductance, target_conductance(target, at_solution, d);
                rtol = 64eps(FT))
            @test isapprox(x.conductance, inv(inv(d.boundary_co2) + d.diffusivity_ratio / sol.stomatal_conductance +
                inv(d.mesophyll_co2)); rtol = 16eps(FT))
        end
        @test equilibrium_gas_exchange(cyto, c.env, d14, c.solver).stomatal_conductance !=
            equilibrium_gas_exchange(cyto, c.env, d16, c.solver).stomatal_conductance
        # The stand-alone, prognostic and LeafModel forms take the same diffusion.
        prognostic = PrognosticStomata(target, FT(300), FT(600))
        @test target_conductance(prognostic, drivers, d14) === g14
        @test target_conductance(cyto, drivers, d14) === g14
        @test stomatal_tendency(prognostic, StomatalState(g14), drivers, d14) == zero(FT)
    end
end
