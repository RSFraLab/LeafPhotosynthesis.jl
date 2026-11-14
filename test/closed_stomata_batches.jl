# Batch checks of the closed-stomata, zero-minimum-conductance and invalid-driver paths.
# Runs on the CPU in the "equilibrium" group (identity) and on CUDA from gpu/runtests.jl, so
# the new branches compile and run inside device kernels and agree with the CPU.

function check_closed_stomata_batches(to_device, synchronize)
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
        leaf, diffusion = core_leaf(FT, p, settings)
        env = core_environment(FT, p, settings)
        solver = core_solver(FT, settings)
        photo = leaf.photosynthesis
        bb = leaf.stomata.target
        forcing(par, co2) = GasExchangeEnvironment(env.temperature, FT(par), FT(co2), env.oxygen,
            env.pressure, env.relative_humidity, env.vpd)
        envs = [forcing(par, co2) for par in (0, 1, 200, 1600), co2 in (2, 40)][:]
        states = [StomatalState(FT(gs)) for gs in (0, 0.2, -0.1, NaN)]
        cases = [(e, s) for e in envs for s in states]
        cpu_exchange = [gas_exchange(leaf, e, s, diffusion) for (e, s) in cases]
        cpu_tendency = [leaf_tendencies(leaf, e, s, diffusion) for (e, s) in cases]
        cpu_water = [transpiration_sensitivity(WaterVaporDiffusion(FT(1.4)), s, e.vpd, e.pressure)
            for (e, s) in cases]
        d_envs, d_states = to_device(first.(cases)), to_device(last.(cases))
        dev_exchange = Array(gas_exchange.(Ref(leaf), d_envs, d_states, Ref(diffusion)))
        dev_tendency = Array(leaf_tendencies.(Ref(leaf), d_envs, d_states, Ref(diffusion)))
        d_vpd = to_device([e.vpd for (e, _) in cases])
        d_pressure = to_device([e.pressure for (e, _) in cases])
        dev_water = Array(transpiration_sensitivity.(Ref(WaterVaporDiffusion(FT(1.4))), d_states,
            d_vpd, d_pressure))
        targets = (BallBerry{FT}(0.0, bb.slope, bb.micromole_to_mole), BallBerry{FT}(0.0, 1.0, bb.micromole_to_mole))
        cpu_eq = [equilibrium_gas_exchange(photo, SteadyStateStomata(t), e, diffusion, solver)
            for t in targets for e in envs]
        dev_eq = vcat((Array(equilibrium_gas_exchange.(Ref(photo), Ref(SteadyStateStomata(t)),
            to_device(envs), Ref(diffusion), Ref(solver))) for t in targets)...)
        drivers = [StomatalDrivers{FT}(a, cs, 101325.0, rh, vpd)
            for a in (-1.0, 0.0, 10.0, NaN), cs in (0.0, 30.0, -1.0), rh in (0.6, 1.5), vpd in (1200.0, NaN)][:]
        cpu_target = [target_conductance(bb, d) for d in drivers]
        dev_target = Array(target_conductance.(Ref(bb), to_device(drivers)))
        synchronize()
        same(a::Real, b::Real) = isnan(a) ? isnan(b) : isapprox(a, b; rtol = 64eps(FT), atol = 64eps(FT))
        same(a::NamedTuple, b::NamedTuple) = all(map(same, a, b))
        @test all(map(same, cpu_exchange, dev_exchange))
        @test all(t -> same(t[1].gas_exchange, t[2].gas_exchange) &&
            same(t[1].conductance_target, t[2].conductance_target) &&
            same(t[1].conductance_tendency, t[2].conductance_tendency), zip(cpu_tendency, dev_tendency))
        @test all(map(same, cpu_water, dev_water))
        @test all(map(same, cpu_target, dev_target))
        @test all(t -> t[1].status === t[2].status && same(t[1].biochemical_co2, t[2].biochemical_co2) &&
            same(t[1].stomatal_conductance, t[2].stomatal_conductance), zip(cpu_eq, dev_eq))
        # The paths are exercised: closed stomata in light, no steady state in darkness,
        # invalid drivers, the closed-stomata status and the invalid-bracket status.
        @test any(x -> isfinite(x.biochemical_co2) && iszero(x.conductance), dev_exchange)
        @test any(x -> isnan(x.biochemical_co2), dev_exchange)
        @test any(isnan, dev_target) && any(isfinite, dev_target)
        @test any(x -> x.status == EquilibriumStomataClosed, dev_eq)
        @test any(x -> x.status == EquilibriumInvalidBracket, dev_eq)
        @test any(x -> x.status == EquilibriumConverged, dev_eq)
    end
end
