# One result shape and one validity predicate, the typed warm start of the equilibrium
# solver, and the abstract component types that composite constructors require.

at_ambient(env, ca) = GasExchangeEnvironment(env.temperature, env.absorbed_par, ca, env.oxygen,
    env.pressure, env.relative_humidity, env.vpd)

function shapes_case(::Type{FT}) where {FT}
    case = setup_cytochrome_stress(FT)
    x = case.settings["experiment"]
    env = GasExchangeEnvironment{FT}(x["temperature_K"], x["high_light_par"],
        x["ambient_co2_Pa"], x["oxygen_Pa"], x["pressure_Pa"], x["relative_humidity"], x["vpd_Pa"])
    return (; case, env, site = LeafEnvironment{FT}(env.temperature, env.absorbed_par, FT(30), env.oxygen))
end

@testset "One result shape: gas exchange and equilibrium share the photosynthesize fields" begin
    for FT in (Float32, Float64)
        (; case, env, site) = shapes_case(FT)
        state = StomatalState(FT(0.2))
        for leaf in (case.fvcb_leaf, case.cyto_leaf)
            p = leaf.photosynthesis
            leaf_result = photosynthesize(p, site)
            exchange = gas_exchange(p, env, state, case.diffusion)
            # The gas-exchange result starts with the fields of `photosynthesize`, flat.
            @test keys(exchange)[1:length(leaf_result)] == keys(leaf_result)
            @test keys(exchange)[length(leaf_result)+1:end] ==
                (:surface_co2, :intercellular_co2, :biochemical_co2, :conductance, :diffusion_residual)
            sol = equilibrium_gas_exchange(leaf, env, case.diffusion, case.solver)
            @test sol.status == EquilibriumConverged
            @test keys(sol.exchange) == keys(exchange)
            @test sol.exchange isa NamedTuple && all(x -> x isa FT, values(sol.exchange))
            # The solution's carbon and electron fields, including the Cytochrome diagnostics,
            # are those of the leaf at the solved biochemical-site CO₂, bit for bit.
            at_solution = photosynthesize(p, LeafEnvironment(env.temperature, env.absorbed_par,
                sol.exchange.biochemical_co2, env.oxygen))
            @test NamedTuple{keys(leaf_result)}(sol.exchange) == at_solution
            # An invalid solution has the same type, with NaN in every exchange field.
            bad = equilibrium_gas_exchange(leaf, at_ambient(env, FT(NaN)), case.diffusion, case.solver)
            @test typeof(bad) == typeof(sol)
            @test bad.status == EquilibriumInvalidInput && all(isnan, values(bad.exchange))
            # A solution is accepted as a reused exchange, through its `exchange`.
            solved = StomatalState(sol.stomatal_conductance)
            @test assimilation_sensitivity(p, env, solved, case.diffusion, BroadbandExcitation(), sol) ==
                assimilation_sensitivity(p, env, solved, case.diffusion, BroadbandExcitation(), sol.exchange)
            @test keys(assimilation_sensitivity(p, env, state, case.diffusion)) ==
                (:dco2_dgs, :dnet_assimilation_dgs)
        end
    end
end

@testset "One validity predicate: is_valid for every result type" begin
    for FT in (Float32, Float64)
        (; case, env, site) = shapes_case(FT)
        state = StomatalState(FT(0.2))
        p = case.cyto_leaf.photosynthesis
        bad_site = LeafEnvironment{FT}(-one(FT), site.absorbed_par, FT(30), site.oxygen)
        @test is_valid(photosynthesize(p, site)) && !is_valid(photosynthesize(p, bad_site))
        @test is_valid(biochemical_turnover(p, site)) && !is_valid(biochemical_turnover(p, bad_site))
        @test is_valid(oxygen_exchange(case.cyto_leaf, site)) && !is_valid(oxygen_exchange(case.cyto_leaf, bad_site))
        @test is_valid(gas_exchange(p, env, state, case.diffusion))
        @test !is_valid(gas_exchange(p, env, StomatalState(FT(-1)), case.diffusion))
        @test is_valid(assimilation_sensitivity(p, env, state, case.diffusion))
        @test !is_valid(assimilation_sensitivity(p, env, StomatalState(FT(-1)), case.diffusion))
        water = WaterVaporDiffusion(FT(1.4))
        @test is_valid(transpiration_sensitivity(water, state, FT(1400), env.pressure))
        @test !is_valid(transpiration_sensitivity(water, state, FT(1400), -env.pressure))
        # Nested results and `nothing` parts: invalid stomatal drivers make the tendencies invalid.
        prognostic = LeafModel(; photosynthesis = p,
            stomata = PrognosticStomata(case.cyto_leaf.stomata.target, FT(600), FT(900)))
        @test is_valid(leaf_tendencies(prognostic, env, state, case.diffusion))
        @test leaf_tendencies(prognostic, env, state, case.diffusion).oxygen === nothing
        # (answer 7: the Medlyn target reads the VPD, not the relative humidity)
        humid = GasExchangeEnvironment{FT}(env.temperature, env.absorbed_par, FT(41), env.oxygen,
            env.pressure, env.relative_humidity, FT(NaN))
        @test !is_valid(leaf_tendencies(prognostic, humid, state, case.diffusion))
        ignored = GasExchangeEnvironment{FT}(env.temperature, env.absorbed_par, FT(41), env.oxygen,
            env.pressure, FT(1.5), env.vpd)
        @test is_valid(leaf_tendencies(prognostic, ignored, state, case.diffusion))
        @test is_valid((; a = one(FT), b = (; c = one(FT), d = nothing)))
        @test !is_valid((; a = one(FT), b = (; c = FT(NaN), d = nothing)))
        # Equilibrium statuses.
        sol = equilibrium_gas_exchange(case.cyto_leaf, env, case.diffusion, case.solver)
        @test is_valid(sol)
        short = EquilibriumSolver(case.solver.method, case.solver.co2_atol, case.solver.co2_rtol,
            case.solver.flux_atol, case.solver.flux_rtol, 1)
        limited = equilibrium_gas_exchange(case.cyto_leaf, env, case.diffusion, short)
        @test limited.status == EquilibriumIterationLimit && !is_valid(limited)
        closed = LeafPhotosynthesis.EquilibriumSolution(sol.exchange, zero(FT), zero(FT), zero(FT), 0, EquilibriumStomataClosed)
        @test is_valid(closed)
        # Fluorescence: lit and dark are valid, a PAM result only when lit.
        inputs = fluorescence_inputs(p, site)
        lit = fluorescence(case.jb_connected, inputs)
        dark = fluorescence(case.jb_connected, CoupledPhotosystemInputs{FT}(0.0, 0.0, 0.0, 0.0, 1.0))
        invalid = fluorescence(case.jb_connected, CoupledPhotosystemInputs{FT}(-1.0, 0.0, 0.0, 0.0, 1.0))
        @test is_valid(lit) && is_valid(dark) && dark.status == FluorescenceDark && !is_valid(invalid)
        c3_inputs = fluorescence_inputs(case.fvcb_leaf.photosynthesis, site)
        @test is_valid(fluorescence(case.cotton_fl, c3_inputs))
        @test is_valid(fluorescence(case.cotton_fl, PSIIFluxInputs{FT}(0.0, 0.0)))
        @test !is_valid(fluorescence(case.cotton_fl, PSIIFluxInputs{FT}(-1.0, 0.0)))
        observer = PAMFluorescenceObserver{FT}(1.0, 0.0)
        a2, a1 = FT(0.5), FT(0.45)
        @test is_valid(pam_fluorescence(observer, case.jb_connected, lit, a2, a1))
        @test !is_valid(pam_fluorescence(observer, case.jb_connected, dark, a2, a1))
    end
end

@testset "Warm start: typed, positional and broadcast per leaf" begin
    for FT in (Float32, Float64)
        (; case, env) = shapes_case(FT)
        leaf, d, solver = case.fvcb_leaf, case.diffusion, case.solver
        envs = [at_ambient(env, FT(ca)) for ca in (30, 41, 60, 90)]
        cold = equilibrium_gas_exchange.(leaf, envs, d, solver)
        @test all(is_valid, cold) && all(s -> s.iterations > 1, cold)
        previous = [s.exchange.biochemical_co2 for s in cold]
        warm = equilibrium_gas_exchange.(leaf, envs, d, solver, BroadbandExcitation(), previous)
        @test all(is_valid, warm)
        # Each leaf starts from its own value: one iteration confirms it.
        @test all(s -> s.iterations == 1, warm)
        @test all(isapprox(w.exchange.biochemical_co2, c.exchange.biochemical_co2;
            atol = 4solver.co2_atol, rtol = 4solver.co2_rtol) for (w, c) in zip(warm, cold))
        # One shared warm start broadcasts as a scalar.
        shared = equilibrium_gas_exchange.(leaf, envs, d, solver, BroadbandExcitation(), previous[2])
        @test shared[2].iterations == 1 && all(is_valid, shared)
        # `nothing` is the cold start; a dual over the model's precision is accepted.
        @test equilibrium_gas_exchange(leaf, env, d, solver, BroadbandExcitation(), nothing) ==
            equilibrium_gas_exchange(leaf, env, d, solver)
        # Another number type is a MethodError, never silently ignored.
        other = FT === Float32 ? Float64 : Float32
        @test_throws MethodError equilibrium_gas_exchange(leaf, env, d, solver, BroadbandExcitation(), 30)
        @test_throws MethodError equilibrium_gas_exchange(leaf, env, d, solver, BroadbandExcitation(), other(30))
        @test_throws MethodError equilibrium_gas_exchange(leaf, env, d, solver, BroadbandExcitation(), "30")
        @test_throws MethodError equilibrium_gas_exchange(leaf.photosynthesis, leaf.stomata, env, d, solver,
            BroadbandExcitation(), 30)
        @test_throws MethodError equilibrium_gas_exchange.(leaf, envs, d, solver, BroadbandExcitation(),
            other.(previous))
    end
end

@testset "Composite constructors require their component types" begin
    case = setup_cytochrome_stress(Float64)
    p = case.fvcb_leaf.photosynthesis
    @test_throws MethodError C3Photosynthesis(1, 2, 3, 4, 5)
    @test_throws MethodError C3Photosynthesis(p.temperature, p.electron_transport, p.stoichiometry,
        SignedAssimilation(), MinimumLimitation())
    # The sign treatment is the field `assimilation_sign`; the former keyword is not accepted.
    @test_throws UndefKeywordError C3Photosynthesis(; p.temperature, p.electron_transport,
        p.stoichiometry, p.limitation, compensation = SignedAssimilation())
    @test_throws MethodError CytochromePhotosynthesis(1, 2, 3, 4, 5)
    @test_throws MethodError C3TemperatureResponse(; vcmax = 1.0, jmax = 1.0, respiration = 1.0,
        kc = 1.0, ko = 1.0, gamma_star = 1.0)
    law = ConstantResponse(1.0)
    @test C3TemperatureResponse(law, law, law, law, law, PrescribedCompensationPoint(law)) isa
        C3TemperatureResponse
    @test_throws MethodError C3TemperatureResponse(law, law, law, law, law, MinimumLimitation())
    @test_throws MethodError CytochromeTemperatureResponse(law, law, law, law, law, law, law, 1.0)
    @test_throws MethodError ElectronTransportParameters(0.5, 0.8, SignedAssimilation())
    @test_throws MethodError CytochromeElectronParameters(0.5, 0.8, 1.0)
    @test_throws MethodError CytochromeCapacity(1.0, 2.0)
    @test_throws MethodError PrescribedCompensationPoint(1.0)
    @test_throws MethodError OxygenScaledCompensationPoint(1.0, 21000.0)
    # Every package component is a subtype of its abstract component type.
    @test MinimumLimitation() isa LeafPhotosynthesis.AbstractColimitation
    @test QuadraticLimitation(0.9) isa LeafPhotosynthesis.AbstractColimitation
    @test NonnegativeAssimilation() isa LeafPhotosynthesis.AbstractAssimilationSign
    @test CytochromeCapacity(1.0, law) isa LeafPhotosynthesis.AbstractTemperatureResponse
    @test PeakedArrhenius(Arrhenius(1.0, 6e4, 298.15, 8.314), 2e5, 650.0) isa
        LeafPhotosynthesis.AbstractTemperatureResponse
end

@testset "LeafModel forms of the stomatal and fluorescence functions" begin
    for FT in (Float32, Float64)
        (; case, env, site) = shapes_case(FT)
        steady = case.cyto_leaf
        drivers = StomatalDrivers{FT}(10.0, 30.0, env.pressure, env.relative_humidity, env.vpd)
        @test target_conductance(steady, drivers, case.diffusion) ==
            target_conductance(steady.stomata.target, drivers, case.diffusion)
        @test fluorescence_inputs(steady, site) == fluorescence_inputs(steady.photosynthesis, site)
        # A steady-state leaf has no tendency.
        @test_throws MethodError leaf_tendencies(steady, env, StomatalState(FT(0.2)), case.diffusion)
    end
end

@testset "CO₂ field names say which CO₂" begin
    for FT in (Float32, Float64)
        site = LeafEnvironment{FT}(; temperature = 298.15, absorbed_par = 800.0,
            biochemical_co2 = 28.0, oxygen = 21000.0)
        env = GasExchangeEnvironment{FT}(; temperature = 298.15, absorbed_par = 800.0,
            ambient_co2 = 40.0, oxygen = 21000.0, pressure = 101325.0, relative_humidity = 0.6,
            vpd = 1200.0)
        @test site == LeafEnvironment{FT}(298.15, 800.0, 28.0, 21000.0)
        @test env == GasExchangeEnvironment{FT}(298.15, 800.0, 40.0, 21000.0, 101325.0, 0.6, 1200.0)
        @test site.biochemical_co2 == FT(28) && env.ambient_co2 == FT(40)
        @test !hasfield(LeafEnvironment, :co2) && !hasfield(GasExchangeEnvironment, :co2)
        @test_throws UndefKeywordError LeafEnvironment{FT}(; temperature = 298.15,
            absorbed_par = 800.0, co2 = 28.0, oxygen = 21000.0)
    end
end

@testset "Temperature-law constructors and integer solver arguments" begin
    @test_throws ArgumentError Arrhenius(NaN, 6e4, 298.15, 8.314)
    @test_throws ArgumentError Arrhenius(1.0, Inf, 298.15, 8.314)
    @test_throws ArgumentError Arrhenius(1.0, 6e4, 0.0, 8.314)
    @test_throws ArgumentError Arrhenius(1.0, 6e4, 298.15, -8.314)
    law = Arrhenius(1.0, 6e4, 298.15, 8.314)
    @test_throws ArgumentError PeakedArrhenius(law, NaN, 650.0)
    @test_throws ArgumentError PeakedArrhenius(law, 2e5, Inf)
    @test PeakedArrhenius(law, 2e5, 650.0) isa PeakedArrhenius{Float64}
    @test SafeguardedNewton(Int32(4)) == SafeguardedNewton(4)
    @test EquilibriumSolver(Bisection(), 1e-6, 0.0, 1e-6, 0.0, Int32(30)).max_iterations === 30
    @test_throws MethodError EquilibriumSolver(Bisection(), 1e-6, 0.0, 1e-6, 0.0, 30.0)
end
