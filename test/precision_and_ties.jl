# Number types, the precision rule, primal decisions at ties, the no-throw rule for values
# that a kernel may see, and the constructor contract.
using ForwardDiff
using ForwardDiff: Dual
using LeafPhotosynthesis: EquilibriumSolution, EquilibriumInvalidInput, EquilibriumStomataClosed,
    FluorescenceLit, CytochromeFluorescenceResult, PAMFluorescenceResult, RateCompetitionResult,
    dark_balanced_psii_partition, state_transition_partition, quenching_rate

function precision_case(::Type{FT}) where {FT}
    bio = read_bonan_parameters(FT, joinpath(dirname(@__DIR__), "examples", "parameters", "bonan_11_1.toml"))
    core = read_bonan_parameters(FT, joinpath(dirname(@__DIR__), "examples", "parameters", "core_exchange.toml"))
    leaf, diffusion = core_leaf(FT, bio, core)
    gas = core_environment(FT, bio, core)
    site = LeafEnvironment(gas.temperature, gas.absorbed_par, FT(28), gas.oxygen)
    steady = LeafModel(; photosynthesis = leaf.photosynthesis,
        stomata = SteadyStateStomata(leaf.stomata.target), oxygen = leaf.oxygen)
    return (; leaf, steady, diffusion, gas, site, solver = core_solver(FT, core),
        state = StomatalState(FT(0.2)), jb = setup_johnson_berry_2021(FT))
end

@testset "Result number type from computed values (dual curvature)" begin
    c = precision_case(Float64)
    p = c.leaf.photosynthesis
    with_curvature(θ) = C3Photosynthesis(p.temperature, p.electron_transport, p.stoichiometry,
        QuadraticLimitation(θ), p.assimilation_sign)
    steady(θ) = LeafModel(; photosynthesis = with_curvature(θ),
        stomata = SteadyStateStomata(c.leaf.stomata.target), oxygen = c.leaf.oxygen)
    for f in (
        θ -> photosynthesize(with_curvature(θ), c.site).net_assimilation,
        θ -> biochemical_turnover(with_curvature(θ), c.site).carboxylation,
        θ -> oxygen_exchange(c.leaf.oxygen, with_curvature(θ), c.site).net_release,
        θ -> equilibrium_gas_exchange(steady(θ), c.gas, c.diffusion, c.solver).exchange.net_assimilation,
    )
        derivative = ForwardDiff.derivative(f, 0.98)
        h = 1e-6
        @test isfinite(derivative)
        @test isapprox(derivative, (f(0.98 + h) - f(0.98 - h)) / 2h; rtol = 1e-5)
    end
    # A Float32 model with a Float64 component is rejected on the host instead of computing
    # in Float64 behind a Float32 result.
    c32 = precision_case(Float32)
    p32 = c32.leaf.photosynthesis
    @test_throws MethodError C3Photosynthesis(p32.temperature, p32.electron_transport,
        p32.stoichiometry, QuadraticLimitation(0.98), p32.assimilation_sign)
    @test_throws MethodError ElectronTransportParameters(0.5f0, 0.85f0, QuadraticLimitation(0.7))
    @test_throws MethodError PrognosticStomata(BallBerry(0.01, 9.0, 1e-6), 600f0, 900f0)
    @test_throws MethodError LeafModel(p32, c.leaf.stomata, c.leaf.oxygen)
    @test_throws MethodError C3TemperatureResponse(p32.temperature.vcmax, p.temperature.jmax,
        p32.temperature.respiration, p32.temperature.kc, p32.temperature.ko, p32.temperature.gamma_star)
    @test all(x -> x isa Float32, values(photosynthesize(c32.leaf, c32.site)))
end

@testset "One precision rule: MethodError at the entry point" begin
    c64, c32 = precision_case(Float64), precision_case(Float32)
    @test_throws MethodError photosynthesize(c64.leaf, c32.site)
    @test_throws MethodError photosynthesize(c32.leaf, c64.site)
    @test_throws MethodError photosynthesize(c64.leaf, c64.site, PrescribedPSIIExcitation(400f0))
    @test_throws MethodError gas_exchange(c64.leaf, c64.gas, StomatalState(0.2f0), c64.diffusion)
    @test_throws MethodError leaf_tendencies(c64.leaf, c64.gas, c64.state, c32.diffusion)
    @test_throws MethodError equilibrium_gas_exchange(c64.steady, c64.gas, c64.diffusion, c32.solver)
    @test_throws MethodError assimilation_sensitivity(c64.leaf.photosynthesis, c64.gas, c32.state, c64.diffusion)
    @test_throws MethodError transpiration_sensitivity(WaterVaporDiffusion(1.4), c32.state, 1400.0, 101325.0)
    @test_throws MethodError target_conductance(c64.leaf.stomata.target,
        StomatalDrivers(10f0, 30f0, 101325f0, 0.6f0, 1200f0), c64.diffusion)
    @test_throws MethodError target_conductance(c64.leaf.stomata.target,
        StomatalDrivers(10.0, 30.0, 101325.0, 0.6, 1200.0), c32.diffusion)
    @test_throws MethodError temperature_response(c32.leaf.photosynthesis.temperature.kc, 298.15)
    jb64, jb32 = c64.jb, c32.jb
    site32 = LeafEnvironment(298.15f0, 800f0, 25f0, 21000f0)
    @test_throws MethodError state_transition_partition(jb64.rule, jb64.fl_ku2, jb64.photo, site32)
    @test_throws MethodError dark_balanced_psii_partition(jb64.fl_ku2, jb64.photo, 298.15f0, 25f0, 21000f0)
    @test_throws MethodError fluorescence(jb64.fl_ku2, fluorescence_inputs(jb32.photo, site32))
    # Duals over the same float type are accepted, and the check folds away.
    d = Dual{:t}(800.0, 1.0)
    @test isfinite(ForwardDiff.value(photosynthesize(c64.leaf,
        LeafEnvironment(298.15 + zero(d), d, 28.0 + zero(d), 21000.0 + zero(d))).net_assimilation))
    @test LeafPhotosynthesis.same_precision(c64.leaf, LeafEnvironment{typeof(d)}(298.15, d, 28.0, 21000.0))
    @test !LeafPhotosynthesis.same_precision(c64.leaf, c32.site)
end

@testset "Validity and kinks decided on primal values (dual seeds ±1 at ties)" begin
    o = setup_oxygen_response(Float64)
    leaf = o.scaled
    steady = LeafModel(; photosynthesis = leaf.photosynthesis,
        stomata = SteadyStateStomata(leaf.stomata.target), oxygen = leaf.oxygen)
    e = o.environment
    zero_o2(O) = GasExchangeEnvironment(e.temperature + zero(O), e.absorbed_par + zero(O),
        e.ambient_co2 + zero(O), O, e.pressure + zero(O), e.relative_humidity + zero(O), e.vpd + zero(O))
    statuses = map((0.0, Dual{:t}(0.0, 1.0), Dual{:t}(0.0, -1.0))) do O
        equilibrium_gas_exchange(steady, zero_o2(O), o.diffusion, o.solver).status
    end
    @test statuses[1] == statuses[2] == statuses[3] == EquilibriumInvalidInput
    for x in (1.0, Dual{:t}(1.0, -1.0), Dual{:t}(1.0, 1.0))
        @test_throws ArgumentError CytochromeElectronParameters(x, 0.9, SerialLimitation())
    end
    for seed in (1.0, -1.0)
        tie = colimit(Dual{:t}(1.0, seed), Dual{:t}(1.0, 2seed), MinimumLimitation())
        @test ForwardDiff.partials(tie)[1] == seed        # first capacity at a tie
        bb = BallBerry(0.01, 9.0, 1e-6)
        g = target_conductance(bb, StomatalDrivers(Dual{:t}(0.0, seed), 38.0, 1e5, 0.7, 1000.0),
            SeriesDiffusion(1.0, 0.2, 1.6, 1e6))
        @test ForwardDiff.partials(g)[1] == 0              # g0 branch at A = 0
        r = biochemical_rates(leaf.photosynthesis.temperature, 298.15, 21000.0)
        cc = Dual{:t}(r.gamma_star, seed)
        net = carbon_assimilation(r, cc, 21000.0, 100.0, leaf.photosynthesis.stoichiometry,
            MinimumLimitation(), NonnegativeAssimilation()).net_assimilation
        @test ForwardDiff.partials(net)[1] == 0            # clipped side at Cc = Γ*
    end
end

@testset "Nothing callable in a kernel throws for a value" begin
    for FT in (Float32, Float64)
        c = precision_case(FT)
        medlyn = Medlyn{FT}(0.01, 126.5, 1e-6, 1.0)
        @test isfinite(target_conductance(medlyn, StomatalDrivers{FT}(10.0, 30.0, 101325.0, 0.6, -5.0), c.diffusion))
        bad_env = GasExchangeEnvironment(c.gas.temperature, c.gas.absorbed_par, c.gas.ambient_co2,
            c.gas.oxygen, c.gas.pressure, c.gas.relative_humidity, FT(NaN))
        @test equilibrium_gas_exchange(c.leaf.photosynthesis, SteadyStateStomata(medlyn), bad_env,
            c.diffusion, c.solver).status == EquilibriumInvalidInput
        @test isnan(quenching_rate(SaturationQuenching(FT(5), FT(2), FT(10)), FT(-0.1)))
        # Per-leaf quantities store any value; their entry points return NaN.
        @test SeriesDiffusion{FT}(-1.0, -1.0, 0.0, 1e6) isa SeriesDiffusion{FT}
        @test WaterVaporDiffusion(FT(-1)) isa WaterVaporDiffusion{FT}
        @test StomatalState(FT(-1)) isa StomatalState{FT}
        # Fixed parameters are validated on the host.
        @test_throws ArgumentError Medlyn{FT}(0.01, 126.5, 1e-6, -1.0)
        @test_throws ArgumentError Medlyn{FT}(0.01, 126.5, 1e-6, 0.0)
        @test_throws ArgumentError Medlyn{FT}(-0.01, 126.5, 1e-6, 1.0)
        @test_throws ArgumentError BallBerry{FT}(-0.01, -9.0, 1e-6)
        @test_throws ArgumentError BallBerry{FT}(0.01, 9.0, 0.0)
        @test_throws ArgumentError C3Stoichiometry{FT}(0.0, -2.0)
        @test_throws ArgumentError CytochromeStoichiometry{FT}(4.0, 8.0, 0.0, 7.0)
        @test_throws ArgumentError CoreOxygenExchange{FT}(1.5, 0.5, 0.0)
        @test_throws ArgumentError CoreOxygenExchange{FT}(1.5, 0.0, 1.0)
        @test BallBerry{FT}(0.0, 9.0, 1e-6).minimum_conductance == zero(FT)
    end
end

@testset "Constructor contract: keyword forms, no discrete values" begin
    @test_throws MethodError Arrhenius{Float64}(1, 2, 3, 4)
    @test_throws MethodError Arrhenius{Float64}(true, 2.0, 3.0, 4.0)
    @test_throws MethodError LeafEnvironment{Float64}(300, 800.0, 30.0, 21000.0)
    @test_throws MethodError StomatalState{Float32}(1)
    @test_throws MethodError StomatalState{Float64}(; conductance = 1//5)
    @test StomatalState(; conductance = 0.2) === StomatalState(0.2)
    @test StomatalState{Float32}(; conductance = 0.2) === StomatalState(0.2f0)
    @test PrescribedPSIIExcitation(; photons = 100.0) === PrescribedPSIIExcitation(100.0)
    @test WaterVaporDiffusion(; boundary_water = 1.4) === WaterVaporDiffusion(1.4)
    gamma = ConstantResponse(4.0)
    @test PrescribedCompensationPoint(; response = gamma) === PrescribedCompensationPoint(gamma)
    @test OxygenScaledCompensationPoint(; response = gamma, reference_oxygen = 21000.0) ===
        OxygenScaledCompensationPoint(gamma, 21000.0)
    target = BallBerry(0.01, 9.0, 1e-6)
    @test SteadyStateStomata(; target) === SteadyStateStomata(target)
    @test RateCompetitionResult(; status = FluorescenceLit, boundary_roundoff = false,
        photon_flux = 1.0, photochemical_yield = 0.5, fluorescence_yield = 0.02,
        constitutive_heat_yield = 0.3, regulated_heat_yield = 0.18, photochemical_rate = 4.0,
        quenching_rate = 1.0) isa RateCompetitionResult{Float64}
    @test PAMFluorescenceResult(; status = FluorescenceLit, fs = 1.0, fm_dark = 2.0,
        fo_dark = 0.5, fm_prime = 1.5, fo_prime = 0.4, apparent_photochemical_yield = 0.3,
        apparent_regulated_heat_yield = 0.2, apparent_residual_yield = 0.5,
        apparent_npq = 0.3, qp = 0.6, ql = 0.4) isa PAMFluorescenceResult{Float64}
    names = fieldnames(CytochromeFluorescenceResult)
    values_ = (FluorescenceLit, false, ntuple(_ -> 0.5, length(names) - 2)...)
    @test CytochromeFluorescenceResult(; NamedTuple{names}(values_)...) isa CytochromeFluorescenceResult{Float64}
    exchange = (; net_assimilation = 1.0, electron_transport = 1.0, surface_co2 = 1.0,
        intercellular_co2 = 1.0, biochemical_co2 = 1.0, conductance = 1.0,
        diffusion_residual = 0.0)
    @test EquilibriumSolution(; exchange, stomatal_conductance = 1.0, co2_residual = 0.0,
        bracket_width = 0.0, iterations = 1, status = EquilibriumStomataClosed) isa
        EquilibriumSolution{Float64}
end

# Specialize on the function and every argument type (Julia 1.10 does not for a plain vararg).
closed_path_allocations(f::F, args::Vararg{Any,N}) where {F,N} = (f(args...); @allocated f(args...))

@testset "Closed-stomata paths: inferred and allocation-free, floats and duals" begin
    for FT in (Float32, Float64)
        c = precision_case(FT)
        p, d, solver = c.leaf.photosynthesis, c.diffusion, c.solver
        bb = c.leaf.stomata.target
        d0 = SteadyStateStomata(BallBerry{FT}(0.0, bb.slope, bb.micromole_to_mole))
        shut = SteadyStateStomata(BallBerry{FT}(0.0, 1.0, bb.micromole_to_mole))
        w = WaterVaporDiffusion(FT(1.4))
        for g in (c.gas, LeafPhotosynthesis.GasExchangeEnvironment(c.gas.temperature,
                ForwardDiff.Dual{:t}(c.gas.absorbed_par, one(FT)), c.gas.ambient_co2, c.gas.oxygen,
                c.gas.pressure, c.gas.relative_humidity, c.gas.vpd))
            closed = StomatalState(zero(typeof(g.absorbed_par)))
            calls = (
                (gas_exchange, (p, g, closed, d)),
                (assimilation_sensitivity, (p, g, closed, d)),
                (leaf_tendencies, (c.leaf, g, closed, d)),
                (equilibrium_gas_exchange, (p, d0, g, d, solver)),
                (equilibrium_gas_exchange, (p, shut, g, d, solver)),
                (transpiration_sensitivity, (w, closed, g.vpd, g.pressure)),
            )
            for (f, args) in calls
                @inferred f(args...)
                @test closed_path_allocations(f, args...) == 0
            end
        end
    end
end
