# Independent high-precision demand and bracketed diffusion/stomatal closure.
# Freeze the input FT thermal rates so root errors are separated from exp rounding.
function reference_demand(p, env)
    rates = biochemical_rates(p.temperature, env.temperature)
    r = map(BigFloat, (rates.vcmax, rates.kc, rates.ko, rates.gamma_star, rates.respiration))
    vc, kc, ko, gamma, rd = r
    electrons = BigFloat(electron_transport(p.electron_transport, env.absorbed_par, rates.jmax))
    ec, oc = BigFloat(p.stoichiometry.electrons_per_carboxylation), BigFloat(p.stoichiometry.oxygenation_cost)
    function demand(c)
        a, b = vc / (c + kc * (1 + BigFloat(env.oxygen) / ko)), electrons / (ec * (c + oc * gamma))
        capacity = if p.limitation isa MinimumLimitation
            min(a, b)
        else
            theta = BigFloat(p.limitation.curvature)
            (a + b - sqrt((a + b)^2 - 4theta*a*b)) / (2theta)
        end
        excess = p.compensation isa SignedAssimilation ? c - gamma : max(c - gamma, 0)
        return capacity * excess - rd
    end
    return demand
end

function reference_equilibrium(p, target, env, diffusion; fixed_gs = nothing)
    demand = reference_demand(p, env)
    ca, pressure, conversion = BigFloat(env.co2), BigFloat(env.pressure), BigFloat(diffusion.mole_to_micromole)
    rb, rm, ratio = inv(BigFloat(diffusion.boundary_co2)), inv(BigFloat(diffusion.mesophyll_co2)), BigFloat(diffusion.diffusivity_ratio)
    g0 = isnothing(fixed_gs) ? BigFloat(target.minimum_conductance) : BigFloat(fixed_gs)
    function residual(c)
        a = demand(c)
        cs = ca - a * pressure / conversion * rb
        if isnothing(fixed_gs) && a > 0 && cs <= 0
            return -one(c)
        end
        gs = g0
        if isnothing(fixed_gs) && a > 0
            slope = BigFloat(target.slope)
            factor = target isa BallBerry ? slope * BigFloat(env.relative_humidity) :
                BigFloat(target.diffusivity_ratio) * (1 + slope / sqrt(max(BigFloat(env.vpd), BigFloat(target.vpd_floor))))
            gs += factor * a * BigFloat(target.micromole_to_mole) * pressure / cs
        end
        return ca - c - a * pressure / conversion * (rb + ratio / gs + rm)
    end
    lo, hi = BigFloat(0), ca + max(-demand(BigFloat(0)), 0) * pressure / conversion * (rb + ratio / g0 + rm)
    for _ in 1:200
        mid = (lo + hi) / 2
        if residual(mid) > 0
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2, demand((lo + hi) / 2)
end

function solver_allocations(p, s, env, d, solver)
    equilibrium_gas_exchange(p, s, env, d, solver)
    return @allocated equilibrium_gas_exchange(p, s, env, d, solver)
end

function test_targets(FT, leaf, settings, d)
    s = settings["stomata"]
    return (leaf.stomata.target, Medlyn{FT}(; minimum_conductance = s["minimum_conductance"],
        slope = s["medlyn_slope"], diffusivity_ratio = d.diffusivity_ratio,
        micromole_to_mole = s["micromole_to_mole"], vpd_floor = s["vpd_floor"]))
end

@testset "Equilibrium and analytic sensitivities" begin
    for FT in (Float32, Float64)
        p = read_bonan_parameters(FT, parameter_path)
        settings = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "examples", "parameters", "core_exchange.toml"))
        leaf, d = core_leaf(FT, p, settings)
        env = core_environment(FT, p, settings)
        solver = core_solver(FT, settings)
        bisector = core_solver(FT, settings; method = Bisection())
        base = leaf.photosynthesis
        targets = test_targets(FT, leaf, settings, d)
        @test solver_allocations(base, SteadyStateStomata(first(targets)), env, d, solver) == 0
        @test isbitstype(typeof(solver))
        @testset "$FT physical grid" begin
            for target in targets, limitation in (MinimumLimitation(), QuadraticLimitation(FT(0.95))),
                compensation in (SignedAssimilation(), NonnegativeAssimilation())
                model = C3Photosynthesis(; temperature = base.temperature, electron_transport = base.electron_transport,
                    stoichiometry = base.stoichiometry, limitation, compensation)
                for light in FT.((0, 1, 1500)), ca in FT.((0, 4, 40, 100)),
                    boundary in FT.((0.005, 1, Inf)), vpd in FT.((0, 1200, 10000))
                    forcing = GasExchangeEnvironment(env.temperature, light, ca, env.oxygen,
                        env.pressure, env.relative_humidity, vpd)
                    diffusion = SeriesDiffusion(boundary, d.mesophyll_co2, d.diffusivity_ratio, d.mole_to_micromole)
                    steady = SteadyStateStomata(target)
                    result = @inferred equilibrium_gas_exchange(model, steady, forcing, diffusion, solver)
                    @test result.converged
                    @test result.status == EquilibriumConverged
                    @test result.biochemical_co2 isa FT
                    @test result.stomatal_conductance isa FT
                    @test result.biochemical_co2 >= zero(FT)
                    @test all(x -> x isa FT, values(result.carbon))
                    cc, assimilation = reference_equilibrium(model, target, forcing, diffusion)
                    tol_c = solver.co2_atol + solver.co2_rtol * max(abs(ca), abs(result.biochemical_co2))
                    tol_a = solver.flux_atol + solver.flux_rtol * abs(result.carbon.net_assimilation)
                    @test isapprox(result.biochemical_co2, cc; atol = 4tol_c, rtol = 64eps(FT))
                    @test isapprox(result.carbon.net_assimilation, assimilation; atol = 4tol_a, rtol = 64eps(FT))
                    @test abs(result.diffusion_residual) <= tol_a * FT(1.01)
                    drivers = StomatalDrivers(result.carbon.net_assimilation, result.surface_co2,
                        forcing.pressure, forcing.relative_humidity, forcing.vpd)
                    @test isapprox(result.stomatal_conductance, target_conductance(target, drivers); rtol = 16eps(FT))
                    b = equilibrium_gas_exchange(model, steady, forcing, diffusion, bisector)
                    @test b.converged
                    @test isapprox(result.biochemical_co2, b.biochemical_co2; atol = 8tol_c, rtol = 64eps(FT))
                    if iszero(light)
                        @test result.iterations == 0
                        @test result.biochemical_co2 >= ca
                    end
                end
            end
        end
        @testset "$FT solver failures and warm starts" begin
            steady = SteadyStateStomata(first(targets))
            solution = equilibrium_gas_exchange(base, steady, env, d, solver)
            warm = equilibrium_gas_exchange(base, steady, env, d, solver; initial_co2 = solution.biochemical_co2)
            @test warm.converged && warm.iterations == 1
            for guess in FT.((NaN, -1, Inf, 1e9))
                @test equilibrium_gas_exchange(base, steady, env, d, solver; initial_co2 = guess) == solution
            end
            short = EquilibriumSolver(solver.method, solver.co2_atol, solver.co2_rtol,
                solver.flux_atol, solver.flux_rtol, 1)
            failed = equilibrium_gas_exchange(base, steady, env, d, short)
            @test !failed.converged && failed.status == EquilibriumIterationLimit
            @test isfinite(failed.carbon.net_assimilation)
            tight = EquilibriumSolver(; method = Bisection(), co2_atol = floatmin(FT), co2_rtol = zero(FT),
                flux_atol = floatmin(FT), flux_rtol = zero(FT), max_iterations = 200)
            # At an exactly representable residual zero, convergence is legitimate;
            # otherwise explicit precision failure is required, never tolerance relaxation.
            limited = equilibrium_gas_exchange(base, steady, env, d, tight)
            @test limited.status in (EquilibriumConverged, EquilibriumPrecisionLimit)
            @test limited.converged == (limited.co2_residual == zero(FT))
            for ca in FT.((-1, NaN, Inf))
                bad = GasExchangeEnvironment(env.temperature, env.absorbed_par, ca, env.oxygen,
                    env.pressure, env.relative_humidity, env.vpd)
                result = equilibrium_gas_exchange(base, steady, bad, d, solver)
                @test result.status == EquilibriumInvalidInput && !result.converged
                @test isnan(result.carbon.net_assimilation)
            end
            zero_target = BallBerry(zero(FT), first(targets).slope, first(targets).micromole_to_mole)
            @test equilibrium_gas_exchange(base, SteadyStateStomata(zero_target), env, d, solver).status == EquilibriumInvalidInput
            for target in targets
                drivers = StomatalDrivers(zero(FT), zero(FT), env.pressure, env.relative_humidity, env.vpd)
                @test target_conductance(target, drivers) == target.minimum_conductance
            end
            @test_throws ArgumentError SafeguardedNewton(; bisection_interval = 0)
        end
        @testset "$FT analytic derivatives" begin
            # Passing a converged exchange also supports smooth/clipped demand,
            # even though the no-exchange analytical shortcut is minimum/signed.
            for compensation in (SignedAssimilation(), NonnegativeAssimilation())
                model = C3Photosynthesis(; temperature = base.temperature, electron_transport = base.electron_transport,
                    stoichiometry = base.stoichiometry, limitation = QuadraticLimitation(FT(0.95)), compensation)
                exchange = equilibrium_gas_exchange(model, SteadyStateStomata(first(targets)), env, d, solver)
                @test exchange.converged
                gs = exchange.stomatal_conductance
                sensitivity = @inferred assimilation_sensitivity(model, env, StomatalState(gs), d, exchange)
                h = BigFloat(gs) * BigFloat("1e-12")
                cp, ap = reference_equilibrium(model, nothing, env, d; fixed_gs = BigFloat(gs) + h)
                cm, am = reference_equilibrium(model, nothing, env, d; fixed_gs = BigFloat(gs) - h)
                @test isapprox(sensitivity.dnet_assimilation_dgs, (ap - am) / (2h);
                    rtol = 16solver.co2_rtol, atol = 16solver.co2_atol)
                @test isapprox(sensitivity.dco2_dgs, (cp - cm) / (2h);
                    rtol = 16solver.co2_rtol, atol = 16solver.co2_atol)
            end
            for limitation in (MinimumLimitation(), QuadraticLimitation(FT(0.95)), QuadraticLimitation(one(FT))),
                compensation in (SignedAssimilation(), NonnegativeAssimilation()), light in FT.((0, 1, 1500))
                model = C3Photosynthesis(; temperature = base.temperature, electron_transport = base.electron_transport,
                    stoichiometry = base.stoichiometry, limitation, compensation)
                forcing = GasExchangeEnvironment(env.temperature, light, env.co2, env.oxygen,
                    env.pressure, env.relative_humidity, env.vpd)
                demand = reference_demand(model, forcing)
                r = biochemical_rates(model.temperature, forcing.temperature)
                j = electron_transport(model.electron_transport, light, r.jmax)
                for c in FT.((1, 10, 40, 100))
                    h = BigFloat("1e-15")
                    expected = (demand(BigFloat(c) + h) - demand(BigFloat(c) - h)) / (2h)
                    derivative = @inferred carbon_assimilation_derivative(r, c, env.oxygen, j,
                        model.stoichiometry, limitation, compensation)
                    @test derivative isa FT
                    @test isapprox(derivative, expected; rtol = 64eps(FT), atol = 64eps(FT))
                end
            end
            for light in FT.((0, 1500)), ca in FT.((0, 40)), gs in FT.((0.01, 0.1, 1))
                forcing = GasExchangeEnvironment(env.temperature, light, ca, env.oxygen,
                    env.pressure, env.relative_humidity, env.vpd)
                state = StomatalState(gs)
                sensitivity = @inferred assimilation_sensitivity(base, forcing, state, d)
                h = BigFloat(gs) * BigFloat("1e-12")
                cp, ap = reference_equilibrium(base, nothing, forcing, d; fixed_gs = BigFloat(gs) + h)
                cm, am = reference_equilibrium(base, nothing, forcing, d; fixed_gs = BigFloat(gs) - h)
                @test all(x -> x isa FT, values(sensitivity))
                @test isapprox(sensitivity.dco2_dgs, (cp - cm) / (2h); rtol = 128eps(FT), atol = 128eps(FT))
                @test isapprox(sensitivity.dnet_assimilation_dgs, (ap - am) / (2h); rtol = 128eps(FT), atol = 128eps(FT))
                @test sensitivity.dassimilation_dgs == sensitivity.dnet_assimilation_dgs
                iszero(light) && @test sensitivity.dnet_assimilation_dgs == zero(FT)
                for boundary in FT.((0.01, 1.4, Inf)), vpd in FT.((-100, 0, 1400))
                    water = @inferred transpiration_sensitivity(WaterVaporDiffusion(boundary), state, vpd, env.pressure)
                    flux(g) = inv(inv(BigFloat(boundary)) + inv(g)) * BigFloat(vpd) / BigFloat(env.pressure)
                    expected = (flux(BigFloat(gs) + h) - flux(BigFloat(gs) - h)) / (2h)
                    @test all(x -> x isa FT, values(water))
                    @test isapprox(water.dtranspiration_dgs, expected; rtol = 16eps(FT), atol = eps(FT))
                    @test isapprox(water.transpiration, flux(BigFloat(gs)); rtol = 8eps(FT))
                end
            end
        end
        @testset "$FT prognostic approaches coupled equilibrium" begin
            for target in targets, multiplier in FT.((0.5, 2))
                dynamic = PrognosticStomata(; target, opening_timescale = leaf.stomata.opening_timescale,
                    closing_timescale = leaf.stomata.closing_timescale)
                model = LeafModel(; photosynthesis = base, stomata = dynamic)
                equilibrium = equilibrium_gas_exchange(base, SteadyStateStomata(target), env, d, solver)
                state = StomatalState(multiplier * equilibrium.stomatal_conductance)
                dt = settings["experiment"]["step"]
                duration = 100 * max(dynamic.opening_timescale, dynamic.closing_timescale)
                for _ in dt:dt:duration
                    result = leaf_tendencies(model, env, state, d)
                    state = StomatalState(state.conductance + dt * result.conductance_tendency)
                end
                @test isapprox(state.conductance, equilibrium.stomatal_conductance; rtol = FT === Float32 ? FT(2e-4) : FT(1e-8))
            end
        end
    end
end
