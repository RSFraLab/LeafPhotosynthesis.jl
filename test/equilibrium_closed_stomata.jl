# Zero minimum conductance in the equilibrium solver, the reused-exchange guard of
# assimilation_sensitivity, and the domain of transpiration_sensitivity.

function equilibrium_closed_case(::Type{FT}) where {FT}
    p = read_bonan_parameters(FT, parameter_path)
    settings = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
    leaf, diffusion = core_leaf(FT, p, settings)
    env = core_environment(FT, p, settings)
    return (; leaf, photo = leaf.photosynthesis, diffusion, env, solver = core_solver(FT, settings))
end

with_forcing(env; par = env.absorbed_par, ambient_co2 = env.ambient_co2) = GasExchangeEnvironment(
    env.temperature, par, ambient_co2, env.oxygen, env.pressure, env.relative_humidity, env.vpd)

@testset "Equilibrium with zero minimum conductance" begin
    for FT in (Float32, Float64)
        c = equilibrium_closed_case(FT)
        photo, env, d, solver = c.photo, c.env, c.diffusion, c.solver
        bb = c.leaf.stomata.target
        targets = (
            BallBerry{FT}(0.0, bb.slope, bb.micromole_to_mole),
            Medlyn{FT}(0.0, 126.49, d.diffusivity_ratio, 1e-6, 1.0),
        )
        closed_exchange = gas_exchange(photo, env, StomatalState(zero(FT)), d)
        gamma = closed_exchange.biochemical_co2
        for target in targets
            steady = SteadyStateStomata(target)
            # Open stomata at the solution: the ordinary root with g_sw = gain A / Cs.
            open = @inferred equilibrium_gas_exchange(photo, steady, env, d, solver)
            @test open.status == EquilibriumConverged && is_valid(open)
            @test open.exchange.biochemical_co2 > gamma && open.stomatal_conductance > zero(FT)
            drivers = StomatalDrivers(open.exchange.net_assimilation, open.exchange.surface_co2,
                env.pressure, env.relative_humidity, env.vpd)
            @test isapprox(open.stomatal_conductance, target_conductance(target, drivers); rtol = 64eps(FT))
            @test abs(open.exchange.diffusion_residual) <= 2 * (solver.flux_atol + solver.flux_rtol * abs(open.exchange.net_assimilation))
            # Ambient CO2 below the compensation point: the closed-stomata limit at Γ.
            low = @inferred equilibrium_gas_exchange(photo, steady, with_forcing(env; ambient_co2 = FT(2)), d, solver)
            @test low.status == EquilibriumStomataClosed && is_valid(low)
            @test low.stomatal_conductance == zero(FT) && low.exchange.conductance == zero(FT)
            @test isapprox(low.exchange.biochemical_co2, gamma; atol = 4 * solver.co2_atol, rtol = 4 * solver.co2_rtol + 64eps(FT))
            @test abs(low.exchange.net_assimilation) <= 4 * (solver.flux_atol + solver.flux_rtol * low.exchange.respiration)
            @test low.exchange.diffusion_residual == -low.exchange.net_assimilation
            # No steady state of a respiring closed leaf: darkness, or light below the
            # light compensation point.
            for par in (zero(FT), FT(1))
                none = equilibrium_gas_exchange(photo, steady, with_forcing(env; par), d, solver)
                @test none.status == EquilibriumInvalidBracket && !is_valid(none)
                @test isnan(none.exchange.net_assimilation)
            end
        end
        # A slope too small to open stomata (Ca − Γ ≤ (P/u) r Ca / gain): closed limit.
        shy = SteadyStateStomata(BallBerry{FT}(0.0, 1.0, bb.micromole_to_mole))
        shut = equilibrium_gas_exchange(photo, shy, env, d, solver)
        @test shut.status == EquilibriumStomataClosed && is_valid(shut)
        @test isapprox(shut.exchange.biochemical_co2, gamma; atol = 4 * solver.co2_atol, rtol = 4 * solver.co2_rtol + 64eps(FT))
        # Limit g0 → 0⁺ of the open root.
        tiny = equilibrium_gas_exchange(photo, SteadyStateStomata(BallBerry{FT}(1e-7, bb.slope, bb.micromole_to_mole)), env, d, solver)
        open = equilibrium_gas_exchange(photo, SteadyStateStomata(first(targets)), env, d, solver)
        @test isapprox(tiny.exchange.biochemical_co2, open.exchange.biochemical_co2; rtol = FT(1e-4))
    end
end

@testset "A reused exchange must belong to the passed excitation and state" begin
    for FT in (Float32, Float64)
        c = equilibrium_closed_case(FT)
        photo, env, d = c.photo, c.env, c.diffusion
        state = StomatalState(FT(0.2))
        other = gas_exchange(photo, env, state, d, PrescribedPSIIExcitation(FT(100)))
        @test all(isnan, values(assimilation_sensitivity(photo, env, state, d, BroadbandExcitation(), other)))
        own = gas_exchange(photo, env, state, d)
        @test assimilation_sensitivity(photo, env, state, d, BroadbandExcitation(), own) ==
            assimilation_sensitivity(photo, env, state, d)
        sol = equilibrium_gas_exchange(photo, SteadyStateStomata(c.leaf.stomata.target), env, d, c.solver)
        @test all(isfinite, values(assimilation_sensitivity(photo, env,
            StomatalState(sol.stomatal_conductance), d, BroadbandExcitation(), sol)))
        @test all(isnan, values(assimilation_sensitivity(photo, env,
            StomatalState(2 * sol.stomatal_conductance), d, BroadbandExcitation(), sol)))
    end
end

@testset "transpiration_sensitivity domain and closed stomata" begin
    for FT in (Float32, Float64)
        w = WaterVaporDiffusion(FT(1.4))
        vpd, pressure = FT(1400), FT(101325)
        closed = @inferred transpiration_sensitivity(w, StomatalState(zero(FT)), vpd, pressure)
        @test closed.transpiration == zero(FT) && closed.dtranspiration_dgs == vpd / pressure
        open = transpiration_sensitivity(w, StomatalState(FT(1e-6)), vpd, pressure)
        @test isapprox(open.dtranspiration_dgs, vpd / pressure; rtol = FT(1e-5))
        @test transpiration_sensitivity(WaterVaporDiffusion(FT(Inf)), StomatalState(FT(0.2)), vpd, pressure).transpiration ==
            FT(0.2) * vpd / pressure
        for (wd, gs, v, p) in (
            (w, FT(-0.1), vpd, pressure), (w, FT(NaN), vpd, pressure), (w, FT(Inf), vpd, pressure),
            (w, FT(0.2), FT(NaN), pressure), (w, FT(0.2), FT(Inf), pressure),
            (w, FT(0.2), vpd, FT(-1)), (w, FT(0.2), vpd, zero(FT)), (w, FT(0.2), vpd, FT(Inf)),
            (WaterVaporDiffusion(FT(-1)), FT(0.2), vpd, pressure),
            (WaterVaporDiffusion(zero(FT)), FT(0.2), vpd, pressure),
            (WaterVaporDiffusion(FT(NaN)), FT(0.2), vpd, pressure),
        )
            result = @inferred transpiration_sensitivity(wd, StomatalState(gs), v, p)
            @test all(x -> x isa FT && isnan(x), values(result))
        end
    end
end

@testset "Closed-stomata CPU batches" begin
    check_closed_stomata_batches(identity, () -> nothing)
end
