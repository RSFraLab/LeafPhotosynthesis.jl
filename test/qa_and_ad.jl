using LeafPhotosynthesis
using Test

const _HAS_EXTENDED_QA_DEPS = all(
    pkg -> Base.find_package(pkg) !== nothing,
    ("Aqua", "Documenter", "ExplicitImports", "ForwardDiff", "OffsetArrays"),
)
if dirname(something(Base.active_project(), "")) == @__DIR__
    @test _HAS_EXTENDED_QA_DEPS
end
if _HAS_EXTENDED_QA_DEPS
    using Aqua
    using Documenter
    using ExplicitImports
    using ForwardDiff
    using OffsetArrays
end

function _bigfloat_deriv(f, x::Float64; h::BigFloat = BigFloat(10)^(-24))
    return setprecision(BigFloat, 256) do
        xb = BigFloat(x)
        Float64((f(xb + h) - f(xb - h)) / (2 * h))
    end
end

function _bigfloat_right_deriv(f, x::Float64; h::BigFloat = BigFloat(10)^(-24))
    return setprecision(BigFloat, 256) do
        xb = BigFloat(x)
        Float64((f(xb + h) - f(xb)) / h)
    end
end

function _bigfloat_left_deriv(f, x::Float64; h::BigFloat = BigFloat(10)^(-24))
    return setprecision(BigFloat, 256) do
        xb = BigFloat(x)
        Float64((f(xb) - f(xb - h)) / h)
    end
end

function _bigfloat_grad(f, x::AbstractVector{Float64}; h::BigFloat = BigFloat(10)^(-24))
    return setprecision(BigFloat, 256) do
        xb = BigFloat.(x)
        out = Vector{Float64}(undef, length(x))
        for i in eachindex(x)
            xp = copy(xb)
            xm = copy(xb)
            xp[i] += h
            xm[i] -= h
            out[i] = Float64((f(xp) - f(xm)) / (2 * h))
        end
        out
    end
end

@testset "Package QA, Generic Indexing, Broadcast, & Automatic Differentiation" begin
    if _HAS_EXTENDED_QA_DEPS
        @testset "Aqua.jl hygiene" begin
            Aqua.test_all(LeafPhotosynthesis)
        end

        if VERSION >= v"1.11"
            @testset "ExplicitImports.jl hygiene" begin
                @test check_no_implicit_imports(LeafPhotosynthesis) === nothing
                @test check_no_stale_explicit_imports(LeafPhotosynthesis) === nothing
                @test check_all_explicit_imports_via_owners(LeafPhotosynthesis) === nothing
                @test check_all_qualified_accesses_via_owners(LeafPhotosynthesis) === nothing
                @test check_no_self_qualified_accesses(LeafPhotosynthesis) === nothing
            end

            @testset "Public & exported docstring completeness" begin
                @test isempty(Docs.undocumented_names(LeafPhotosynthesis))
            end
        end

        @testset "Documenter.jl doctests" begin
            DocMeta.setdocmeta!(
                LeafPhotosynthesis,
                :DocTestSetup,
                :(using LeafPhotosynthesis);
                recursive = true,
            )
            doctest(LeafPhotosynthesis; manual = false)
        end
    end

    @testset "Inner-constructor invariants, scalar rejection, layout, & broadcastable" begin
        # Positional inner-constructor domain validation
        @test_throws ArgumentError ElectronTransportParameters(1.2, 0.85, MinimumLimitation())
        @test_throws ArgumentError ElectronTransportParameters(-0.1, 0.85, MinimumLimitation())
        @test_throws ArgumentError ElectronTransportParameters(0.5, NaN, MinimumLimitation())

        @test_throws ArgumentError CytochromeElectronParameters(0.0, 0.95, SerialLimitation())
        @test_throws ArgumentError CytochromeElectronParameters(1.0, 0.95, SerialLimitation())
        @test_throws ArgumentError CytochromeElectronParameters(0.5, NaN, SerialLimitation())

        med_valid = Medlyn(0.01, 4.0, 1.6, 1e-6, 50.0)
        @test_throws ArgumentError PrognosticStomata(med_valid, 0.0, 600.0)
        @test_throws ArgumentError PrognosticStomata(med_valid, 300.0, -5.0)
        @test_throws ArgumentError PrognosticStomata(med_valid, Inf, 600.0)

        @test_throws ArgumentError SafeguardedNewton(0)
        @test_throws ArgumentError SafeguardedNewton(-3)
        @test_throws ArgumentError EquilibriumSolver(SafeguardedNewton(4), 1e-6, 0.0, 1e-6, 0.0, 0)
        @test_throws ArgumentError EquilibriumSolver(SafeguardedNewton(4), -1e-6, 0.0, 1e-6, 0.0, 30)
        @test_throws ArgumentError EquilibriumSolver(SafeguardedNewton(4), 0.0, 0.0, 1e-6, 0.0, 30)
        @test_throws ArgumentError EquilibriumSolver(SafeguardedNewton(4), 1e-6, 0.0, NaN, 0.0, 30)

        # Discrete scalar rejection via MethodError across positional, keyword, and {FT} constructors
        for bad in (1, true, 1 // 2, π)
            BT = typeof(bad)
            @test_throws MethodError Arrhenius(bad, bad, bad, bad)
            @test_throws MethodError Arrhenius(; value_ref = bad, activation_energy = bad, temperature_ref = bad, gas_constant = bad)
            @test_throws MethodError Arrhenius{BT}(1.0, 1.0, 1.0, 1.0)
            @test_throws MethodError PeakedArrhenius{BT}(Arrhenius(60.0, 65000.0, 298.15, 8.314), 200000.0, 650.0)
            @test_throws MethodError PeakedArrhenius(Arrhenius(60.0, 65000.0, 298.15, 8.314), bad, bad)
            @test_throws MethodError ConstantResponse(bad)
            @test_throws MethodError C3Rates(bad, bad, bad, bad, bad, bad)
            @test_throws MethodError C3Rates{BT}(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
            @test_throws MethodError QuadraticLimitation(bad)
            @test_throws MethodError PrescribedPSIIExcitation(bad)
            @test_throws MethodError PrescribedPSIIExcitation{BT}(100.0)
            @test_throws MethodError PrescribedDualExcitation(bad, bad)
            @test_throws MethodError PrescribedDualExcitation{BT}(100.0, 100.0)
            @test_throws MethodError ElectronTransportParameters(bad, bad, MinimumLimitation())
            @test_throws MethodError C3Stoichiometry(bad, bad)
            @test_throws MethodError C3Stoichiometry(; electrons_per_carboxylation = bad, oxygenation_cost = bad)
            @test_throws MethodError C3Stoichiometry{BT}(4.0, 4.0)
            @test_throws MethodError OxygenScaledCompensationPoint(ConstantResponse(4.0), bad)
            @test_throws MethodError LeafEnvironment(bad, bad, bad, bad)
            @test_throws MethodError LeafEnvironment(; temperature = bad, absorbed_par = bad, co2 = bad, oxygen = bad)
            @test_throws MethodError LeafEnvironment{BT}(298.15, 800.0, 28.0, 21000.0)
            @test_throws MethodError CoreOxygenExchange(bad, bad, bad, bad)
            @test_throws MethodError CoreOxygenExchange{BT}(4.0, 1.5, 0.5, 1.0)
            @test_throws MethodError OxygenFluxInputs(bad, bad, bad)
            @test_throws MethodError OxygenFluxInputs{BT}(10.0, 2.0, 1.0)
            @test_throws MethodError ConstantQuenching(bad)
            @test_throws MethodError SaturationQuenching(bad, bad, bad)
            @test_throws MethodError PSIIFluxInputs(bad, bad)
            @test_throws MethodError PSIIFluxInputs{BT}(100.0, 40.0)
            @test_throws MethodError RateCompetitionResult(FluorescenceLit, false, ntuple(_ -> bad, 7)...)
            @test_throws MethodError RateCompetitionResult{BT}(FluorescenceLit, false, ntuple(_ -> 0.5, 7)...)
            @test_throws MethodError StomatalDrivers(bad, bad, bad, bad, bad)
            @test_throws MethodError StomatalDrivers{BT}(10.0, 35.0, 101325.0, 0.6, 1200.0)
            @test_throws MethodError BallBerry(bad, bad, bad)
            @test_throws MethodError BallBerry{BT}(0.01, 9.0, 1e-6)
            @test_throws MethodError Medlyn(bad, bad, bad, bad, bad)
            @test_throws MethodError Medlyn{BT}(0.01, 4.0, 1.6, 1e-6, 50.0)
            @test_throws MethodError PrognosticStomata(med_valid, bad, bad)
            @test_throws MethodError StomatalState(bad)
            @test_throws MethodError StomatalState{BT}(0.2)
            @test_throws MethodError SeriesDiffusion(bad, bad, bad, bad)
            @test_throws MethodError SeriesDiffusion(; boundary_co2 = bad, mesophyll_co2 = bad, diffusivity_ratio = bad, mole_to_micromole = bad)
            @test_throws MethodError SeriesDiffusion{BT}(1.4, 0.2, 1.6, 1e6)
            @test_throws MethodError GasExchangeEnvironment(bad, bad, bad, bad, bad, bad, bad)
            @test_throws MethodError GasExchangeEnvironment{BT}(298.15, 800.0, 40.0, 21000.0, 101325.0, 0.6, 1200.0)
            @test_throws MethodError WaterVaporDiffusion(bad)
            @test_throws MethodError WaterVaporDiffusion{BT}(1.4)
            @test_throws MethodError EquilibriumSolver(SafeguardedNewton(4), bad, bad, bad, bad, 30)
            @test_throws MethodError EquilibriumSolution(
                (net_assimilation = 1.0,), ntuple(_ -> bad, 9)..., 1, EquilibriumConverged, true,
            )
            @test_throws MethodError EquilibriumSolution{BT,NamedTuple{(:net_assimilation,),Tuple{Float64}}}(
                (net_assimilation = 1.0,), ntuple(_ -> 1.0, 9)..., 1, EquilibriumConverged, true,
            )
            @test_throws MethodError CytochromeCapacity(bad, ConstantResponse(300.0))
            @test_throws MethodError CytochromeRates(bad, bad, bad, bad, bad, bad, bad, bad)
            @test_throws MethodError CytochromeRates{BT}(60.0, 150.0, 1.0, 40.0, 27000.0, 4.0, 1.28, 1.0)
            @test_throws MethodError CytochromeStoichiometry(bad, bad, bad, bad)
            @test_throws MethodError CytochromeStoichiometry{BT}(4.0, 8.0, 3.0, 3.5)
            @test_throws MethodError CytochromeElectronParameters(bad, bad, SerialLimitation())
            @test_throws MethodError PSIIRateParameters(bad, bad, bad, bad, bad)
            @test_throws MethodError PSIRateParameters(bad, bad, bad, bad)
            @test_throws MethodError CoupledPhotosystemInputs(bad, bad, bad, bad, bad)
            @test_throws MethodError CoupledPhotosystemInputs{BT}(100.0, 100.0, 40.0, 50.0, 0.6)
            @test_throws MethodError CytochromeFluorescenceResult(FluorescenceLit, false, ntuple(_ -> bad, 16)...)
            @test_throws MethodError CytochromeFluorescenceResult{BT}(FluorescenceLit, false, ntuple(_ -> 0.5, 16)...)
            @test_throws MethodError StateTransitionRule(bad, bad)
            @test_throws MethodError PAMFluorescenceObserver(bad, bad)
            @test_throws MethodError PAMFluorescenceResult(FluorescenceLit, ntuple(_ -> bad, 11)...)
            @test_throws MethodError PAMFluorescenceResult{BT}(FluorescenceLit, ntuple(_ -> 0.5, 11)...)
        end
        # Positive controls confirming exact arity/order of result constructors on Float64
        @test RateCompetitionResult(FluorescenceLit, false, ntuple(_ -> 0.5, 7)...) isa RateCompetitionResult{Float64}
        @test EquilibriumSolution((net_assimilation = 1.0,), ntuple(_ -> 1.0, 9)..., 1, EquilibriumConverged, true) isa EquilibriumSolution{Float64}
        @test CytochromeFluorescenceResult(FluorescenceLit, false, ntuple(_ -> 0.5, 16)...) isa CytochromeFluorescenceResult{Float64}
        @test PAMFluorescenceResult(FluorescenceLit, ntuple(_ -> 0.5, 11)...) isa PAMFluorescenceResult{Float64}

        # Required-keyword constructors still reject empty call with UndefKeywordError
        @test_throws UndefKeywordError StateTransitionRule()
        @test_throws UndefKeywordError PAMFluorescenceObserver()

        # Mixed Float32 / Float64 precision rejection on promoting constructors (exact arities)
        @test_throws MethodError Arrhenius(60.0f0, 65000.0, 298.15, 8.314)
        @test_throws MethodError LeafEnvironment(298.15f0, 800.0, 28.0, 21000.0)
        @test_throws MethodError GasExchangeEnvironment(298.15f0, 800.0, 40.0, 21000.0, 101325.0, 0.6, 1200.0)
        @test_throws MethodError PrescribedDualExcitation(300.0f0, 300.0)
        @test_throws MethodError Medlyn(0.01f0, 4.0, 1.6, 1e-6, 50.0)
        @test_throws MethodError SeriesDiffusion(1.4f0, 0.1, 1.6, 1e6)
        @test_throws MethodError SeriesDiffusion(; boundary_co2 = 1.4f0, mesophyll_co2 = 0.1, diffusivity_ratio = 1.6, mole_to_micromole = 1e6)
        @test_throws MethodError PSIIRateParameters(0.05f0, 0.55, 4.0, 2.0, 0.0)
        @test_throws MethodError PSIRateParameters(0.05f0, 0.55, 14.5, 5.0)

        # Packed isbitstype layout preserved
        @test isbitstype(LeafEnvironment{Float64}) && sizeof(LeafEnvironment{Float64}) == 32
        @test isbitstype(LeafEnvironment{Float32}) && sizeof(LeafEnvironment{Float32}) == 16
        @test isbitstype(GasExchangeEnvironment{Float64}) && sizeof(GasExchangeEnvironment{Float64}) == 56
        @test isbitstype(GasExchangeEnvironment{Float32}) && sizeof(GasExchangeEnvironment{Float32}) == 28

        # Direct scalar broadcasting without Ref(...) on standard Arrays & @view slices
        lp_root = dirname(@__DIR__)
        p64 = read_bonan_parameters(Float64, joinpath(lp_root, "examples", "parameters", "bonan_11_1.toml"))
        c3_model = bonan_leaf(Float64, p64)
        sites_std = [
            LeafEnvironment(293.15, 0.0, 30.0, 21000.0),
            LeafEnvironment(298.15, 600.0, 28.0, 21000.0),
            LeafEnvironment(303.15, 1200.0, 24.0, 21000.0),
        ]
        @test photosynthesize.(c3_model, sites_std) == photosynthesize.(Ref(c3_model), sites_std)
        @test photosynthesize.(c3_model, @view(sites_std[1:2:3])) == photosynthesize.(Ref(c3_model), sites_std[1:2:3])
    end

    @testset "Generated API docs sync" begin
        api_mod = Module(:_LeafPhotosynthesisAPIDocsCheck)
        Base.include(api_mod, joinpath(dirname(@__DIR__), "docs", "api_catalog.jl"))
        @test length(Base.invokelatest(getfield, api_mod, :API_PAGES)) == 6
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(@__DIR__)) $(joinpath(dirname(@__DIR__), "docs", "build_api_docs.jl")) --check`
        @test success(cmd)
    end

    if _HAS_EXTENDED_QA_DEPS
    @testset "OffsetArrays, @view slices, and Base.broadcastable scalar models" begin
        lp_root = dirname(@__DIR__)
        p64 = read_bonan_parameters(Float64, joinpath(lp_root, "examples", "parameters", "bonan_11_1.toml"))
        s64 = read_bonan_parameters(Float64, joinpath(lp_root, "examples", "parameters", "core_exchange.toml"))
        c3_prog, diff64 = core_leaf(Float64, p64, s64)
        c3_steady = LeafModel(;
            photosynthesis = c3_prog.photosynthesis,
            stomata = SteadyStateStomata(c3_prog.stomata.target),
            oxygen = c3_prog.oxygen,
        )
        solver64 = core_solver(Float64, s64)
        env_base = core_environment(Float64, p64, s64)
        vdt_case = setup_fluorescence(Float64)
        if !isdefined(@__MODULE__, :setup_johnson_berry_2021)
            include(joinpath(lp_root, "examples", "setup", "johnson_berry_2021.jl"))
        end
        jb = setup_johnson_berry_2021(Float64)
        cyto_steady = LeafModel(;
            photosynthesis = jb.photo,
            stomata = SteadyStateStomata(Medlyn(0.01, 4.0, 1.6, 1e-6, 50.0)),
            oxygen = CoreOxygenExchange(4.0, 1.5, 0.5, 1.0),
        )

        sites_vec = [
            LeafEnvironment(293.15, 0.0, 30.0, 21000.0),
            LeafEnvironment(298.15, 250.0, 28.0, 21000.0),
            LeafEnvironment(303.15, 900.0, 24.0, 21000.0),
            LeafEnvironment(308.15, 1400.0, 18.0, 21000.0),
            LeafEnvironment(288.15, 50.0, 36.0, 21000.0),
            LeafEnvironment(300.15, 600.0, 32.0, 21000.0),
        ]
        sites_off1 = OffsetArray(sites_vec, -3)
        sites_off2 = OffsetArray(reshape(sites_vec, 2, 3), 0:1, 10:12)
        sites_view = @view sites_vec[1:2:end]

        # Direct broadcasting WITHOUT Ref(...) verifies Base.broadcastable on scalar leaf components
        res_ref = photosynthesize.(c3_steady, sites_vec)
        res_off1 = photosynthesize.(c3_steady, sites_off1)
        res_off2 = photosynthesize.(c3_steady, sites_off2)
        res_view = photosynthesize.(c3_steady, sites_view)
        @test axes(res_off1) == axes(sites_off1)
        @test axes(res_off2) == axes(sites_off2)
        @test collect(res_off1) == res_ref
        @test vec(collect(res_off2)) == res_ref
        @test res_view == res_ref[1:2:end]

        # In-place broadcast! onto OffsetArray destination
        dest_off2 = similar(res_off2)
        broadcast!(photosynthesize, dest_off2, c3_steady, sites_off2)
        @test dest_off2 == res_off2

        # Turnover and oxygen exchange on OffsetArray & @view
        turn_off1 = biochemical_turnover.(c3_steady, sites_off1)
        @test axes(turn_off1) == axes(sites_off1)
        @test collect(turn_off1) == biochemical_turnover.(c3_steady, sites_vec)
        ox_off2 = oxygen_exchange.(c3_steady, sites_off2)
        @test axes(ox_off2) == axes(sites_off2)
        @test vec(collect(ox_off2)) == oxygen_exchange.(c3_steady, sites_vec)

        # Gas exchange & equilibrium solver on OffsetArray & @view
        gas_vec = [
            GasExchangeEnvironment(s.temperature, s.absorbed_par, 40.0, s.oxygen, env_base.pressure, 0.65, 1200.0)
            for s in sites_vec
        ]
        gs_vec = fill(StomatalState(0.18), length(sites_vec))
        gas_off1 = OffsetArray(gas_vec, 5)
        gs_off1 = OffsetArray(gs_vec, 5)
        gas_view = @view gas_vec[2:2:end]
        gs_view = @view gs_vec[2:2:end]

        ge_off1 = gas_exchange.(c3_steady, gas_off1, gs_off1, diff64)
        @test axes(ge_off1) == axes(gas_off1)
        @test collect(ge_off1) == gas_exchange.(c3_steady, gas_vec, gs_vec, diff64)
        @test gas_exchange.(c3_steady, gas_view, gs_view, diff64) ==
              gas_exchange.(c3_steady, gas_vec[2:2:end], gs_vec[2:2:end], diff64)

        eq_c3_off = equilibrium_gas_exchange.(c3_steady, gas_off1, diff64, solver64)
        eq_cyto_off = equilibrium_gas_exchange.(cyto_steady, gas_off1, diff64, solver64)
        @test axes(eq_c3_off) == axes(gas_off1)
        @test collect(eq_c3_off) == equilibrium_gas_exchange.(c3_steady, gas_vec, diff64, solver64)
        @test collect(eq_cyto_off) == equilibrium_gas_exchange.(cyto_steady, gas_vec, diff64, solver64)
        @test equilibrium_gas_exchange.(cyto_steady, gas_view, diff64, solver64) ==
              equilibrium_gas_exchange.(cyto_steady, gas_vec[2:2:end], diff64, solver64)

        # Fluorescence, state transitions, and PAM observer on OffsetArray & @view
        vdt_off = fluorescence.(vdt_case.regulated, psii_flux_inputs.(vdt_case.leaf.photosynthesis, sites_off1))
        @test axes(vdt_off) == axes(sites_off1)
        @test collect(vdt_off) == fluorescence.(vdt_case.regulated, psii_flux_inputs.(vdt_case.leaf.photosynthesis, sites_vec))

        exc_off = state_transition_excitation.(jb.rule, jb.fl_ku2, jb.photo, sites_off2)
        @test axes(exc_off) == axes(sites_off2)
        jb_fl_off = fluorescence.(jb.fl_ku2, photosystem_flux_inputs.(jb.photo, sites_off2, exc_off))
        pam_off = pam_fluorescence.(jb.pam_leaked, jb.fl_ku2, jb_fl_off, 0.425, 0.425)
        @test axes(pam_off) == axes(sites_off2)
        exc_ref = state_transition_excitation.(jb.rule, jb.fl_ku2, jb.photo, sites_vec)
        jb_fl_ref = fluorescence.(jb.fl_ku2, photosystem_flux_inputs.(jb.photo, sites_vec, exc_ref))
        @test vec(collect(pam_off)) == pam_fluorescence.(jb.pam_leaked, jb.fl_ku2, jb_fl_ref, 0.425, 0.425)
    end

    @testset "ForwardDiff gradients, parameters, & boundary smoothness vs 256-bit BigFloat" begin
        lp_root = dirname(@__DIR__)
        p64 = read_bonan_parameters(Float64, joinpath(lp_root, "examples", "parameters", "bonan_11_1.toml"))
        s64 = read_bonan_parameters(Float64, joinpath(lp_root, "examples", "parameters", "core_exchange.toml"))
        p_bf = read_bonan_parameters(BigFloat, joinpath(lp_root, "examples", "parameters", "bonan_11_1.toml"))
        s_bf = read_bonan_parameters(BigFloat, joinpath(lp_root, "examples", "parameters", "core_exchange.toml"))

        c3_prog, diff64 = core_leaf(Float64, p64, s64)
        c3_prog_bf, diff_bf = core_leaf(BigFloat, p_bf, s_bf)
        vdt64 = setup_fluorescence(Float64)
        vdt_bf = setup_fluorescence(BigFloat)
        if !isdefined(@__MODULE__, :setup_johnson_berry_2021)
            include(joinpath(lp_root, "examples", "setup", "johnson_berry_2021.jl"))
        end
        jb64 = setup_johnson_berry_2021(Float64)
        jb_bf = setup_johnson_berry_2021(BigFloat)

        # 1. Environment gradients across C3 & Cytochrome photosynthesis + oxygen + turnover
        x_site_rubisco = [298.15, 1200.0, 22.0, 21000.0]
        x_site_light   = [296.15, 180.0, 38.0, 21000.0]
        for x_site in (x_site_rubisco, x_site_light)
            g_c3 = ForwardDiff.gradient(
                x -> photosynthesize(c3_prog.photosynthesis, LeafEnvironment(x[1], x[2], x[3], x[4])).net_assimilation,
                x_site,
            )
            g_c3_ref = _bigfloat_grad(
                x -> photosynthesize(c3_prog_bf.photosynthesis, LeafEnvironment(x[1], x[2], x[3], x[4])).net_assimilation,
                x_site,
            )
            @test all(isfinite, g_c3)
            @test isapprox(g_c3, g_c3_ref; rtol = 1e-11, atol = 1e-13)

            g_ox = ForwardDiff.gradient(
                x -> oxygen_exchange(c3_prog, LeafEnvironment(x[1], x[2], x[3], x[4])).net_release,
                x_site,
            )
            g_ox_ref = _bigfloat_grad(
                x -> oxygen_exchange(c3_prog_bf, LeafEnvironment(x[1], x[2], x[3], x[4])).net_release,
                x_site,
            )
            @test isapprox(g_ox, g_ox_ref; rtol = 1e-11, atol = 1e-13)

            g_cyto = ForwardDiff.gradient(
                x -> photosynthesize(jb64.photo, LeafEnvironment(x[1], x[2], x[3], x[4])).net_assimilation,
                x_site,
            )
            g_cyto_ref = _bigfloat_grad(
                x -> photosynthesize(jb_bf.photo, LeafEnvironment(x[1], x[2], x[3], x[4])).net_assimilation,
                x_site,
            )
            @test all(isfinite, g_cyto)
            @test isapprox(g_cyto, g_cyto_ref; rtol = 1e-11, atol = 1e-13)
        end

        # 2. Equilibrium gas exchange AD (Medlyn & BallBerry; SafeguardedNewton even with coarse tol & Bisection; warm-start)
        med64 = SteadyStateStomata(Medlyn(0.01, 4.0, 1.6, 1e-6, 50.0))
        med_bf = SteadyStateStomata(Medlyn(BigFloat("0.01"), BigFloat(4), BigFloat("1.6"), BigFloat("1e-6"), BigFloat(50)))
        bb64 = SteadyStateStomata(c3_prog.stomata.target)
        bb_bf = SteadyStateStomata(c3_prog_bf.stomata.target)

        coarse_newton = EquilibriumSolver(;
            method = SafeguardedNewton(4),
            co2_atol = 1e-4, co2_rtol = 0.0, flux_atol = 1e-4, flux_rtol = 0.0,
            max_iterations = 40,
        )
        bisect_solver = EquilibriumSolver(;
            method = Bisection(),
            co2_atol = 1e-7, co2_rtol = 0.0, flux_atol = 1e-7, flux_rtol = 0.0,
            max_iterations = 80,
        )
        tight_bf_solver = EquilibriumSolver(;
            method = SafeguardedNewton(4),
            co2_atol = BigFloat(10)^(-32),
            co2_rtol = BigFloat(0),
            flux_atol = BigFloat(10)^(-32),
            flux_rtol = BigFloat(0),
            max_iterations = 80,
        )

        x_eq = [298.15, 450.0, 40.0, 21000.0, 101325.0, 0.62, 1350.0]
        for (stom64, stom_bf) in ((med64, med_bf), (bb64, bb_bf))
            g_ref_c3 = _bigfloat_grad(
                x -> equilibrium_gas_exchange(
                    c3_prog_bf.photosynthesis,
                    stom_bf,
                    GasExchangeEnvironment(x[1], x[2], x[3], x[4], x[5], x[6], x[7]),
                    diff_bf,
                    tight_bf_solver,
                ).carbon.net_assimilation,
                x_eq,
            )
            for slv in (coarse_newton, bisect_solver)
                g_ad = ForwardDiff.gradient(
                    x -> equilibrium_gas_exchange(
                        c3_prog.photosynthesis,
                        stom64,
                        GasExchangeEnvironment(x[1], x[2], x[3], x[4], x[5], x[6], x[7]),
                        diff64,
                        slv,
                    ).carbon.net_assimilation,
                    x_eq,
                )
                @test all(isfinite, g_ad)
                @test isapprox(g_ad, g_ref_c3; rtol = 1e-10, atol = 1e-12)
            end

            # Warm-start initial_co2 must carry zero spurious sensitivity to initial_co2 itself
            # and must not pollute environment AD tangents even when warm-started at the converged root.
            sol_primal = equilibrium_gas_exchange(
                c3_prog.photosynthesis,
                stom64,
                GasExchangeEnvironment(x_eq...),
                diff64,
                coarse_newton,
            )
            d_init = ForwardDiff.derivative(
                ci0 -> equilibrium_gas_exchange(
                    c3_prog.photosynthesis,
                    stom64,
                    GasExchangeEnvironment(x_eq...),
                    diff64,
                    coarse_newton;
                    initial_co2 = ci0,
                ).intercellular_co2,
                25.0,
            )
            @test d_init == 0.0
            g_warm = ForwardDiff.gradient(
                x -> equilibrium_gas_exchange(
                    c3_prog.photosynthesis,
                    stom64,
                    GasExchangeEnvironment(x[1], x[2], x[3], x[4], x[5], x[6], x[7]),
                    diff64,
                    coarse_newton;
                    initial_co2 = sol_primal.biochemical_co2 + 999.0 * (x[1] - x_eq[1]),
                ).carbon.net_assimilation,
                x_eq,
            )
            @test isapprox(g_warm, g_ref_c3; rtol = 1e-10, atol = 1e-12)

            g_ref_cyto = _bigfloat_grad(
                x -> equilibrium_gas_exchange(
                    jb_bf.photo,
                    stom_bf,
                    GasExchangeEnvironment(x[1], x[2], x[3], x[4], x[5], x[6], x[7]),
                    diff_bf,
                    tight_bf_solver,
                ).carbon.net_assimilation,
                x_eq,
            )
            for slv in (coarse_newton, bisect_solver)
                g_ad_cyto = ForwardDiff.gradient(
                    x -> equilibrium_gas_exchange(
                        jb64.photo,
                        stom64,
                        GasExchangeEnvironment(x[1], x[2], x[3], x[4], x[5], x[6], x[7]),
                        diff64,
                        slv,
                    ).carbon.net_assimilation,
                    x_eq,
                )
                @test all(isfinite, g_ad_cyto)
                @test isapprox(g_ad_cyto, g_ref_cyto; rtol = 1e-10, atol = 1e-12)
            end
        end

        # 3. Fluorescence, StateTransitionRule (high light KN > 0 AND low light KN == 0), & PAM observer AD
        g_vdt = ForwardDiff.gradient(
            x -> fluorescence(vdt64.regulated, psii_flux_inputs(vdt64.leaf.photosynthesis, LeafEnvironment(x[1], x[2], x[3], x[4]))).fluorescence_yield,
            x_site_light,
        )
        g_vdt_ref = _bigfloat_grad(
            x -> fluorescence(vdt_bf.regulated, psii_flux_inputs(vdt_bf.leaf.photosynthesis, LeafEnvironment(x[1], x[2], x[3], x[4]))).fluorescence_yield,
            x_site_light,
        )
        @test isapprox(g_vdt, g_vdt_ref; rtol = 1e-11, atol = 1e-13)

        for x_fl in ([298.15, 800.0, 28.0, 21000.0], [298.15, 20.0, 28.0, 21000.0])
            g_qt = ForwardDiff.gradient(
                x -> begin
                    site = LeafEnvironment(x[1], x[2], x[3], x[4])
                    exc = state_transition_excitation(jb64.rule, jb64.fl_ku2, jb64.photo, site)
                    fl = fluorescence(jb64.fl_ku2, photosystem_flux_inputs(jb64.photo, site, exc))
                    pam_fluorescence(jb64.pam_leaked, jb64.fl_ku2, fl, 0.425, 0.425).apparent_photochemical_yield
                end,
                x_fl,
            )
            g_qt_ref = _bigfloat_grad(
                x -> begin
                    site = LeafEnvironment(x[1], x[2], x[3], x[4])
                    exc = state_transition_excitation(jb_bf.rule, jb_bf.fl_ku2, jb_bf.photo, site)
                    fl = fluorescence(jb_bf.fl_ku2, photosystem_flux_inputs(jb_bf.photo, site, exc))
                    pam_fluorescence(jb_bf.pam_leaked, jb_bf.fl_ku2, fl, BigFloat("0.425"), BigFloat("0.425")).apparent_photochemical_yield
                end,
                x_fl,
            )
            @test all(isfinite, g_qt)
            @test isapprox(g_qt, g_qt_ref; rtol = 1e-11, atol = 1e-13)
        end

        # 4. Single-parameter AD without re-wrapping untouched Float64 fields
        site_fixed = LeafEnvironment(298.15, 800.0, 24.0, 21000.0)
        site_fixed_bf = LeafEnvironment(BigFloat("298.15"), BigFloat(800), BigFloat(24), BigFloat(21000))

        fn_vcmax = (vc, base_leaf, site_arg) -> begin
            t0 = base_leaf.photosynthesis.temperature
            a0 = t0.vcmax.arrhenius
            vc_law = PeakedArrhenius(;
                arrhenius = Arrhenius(vc, a0.activation_energy, a0.temperature_ref, a0.gas_constant),
                deactivation_energy = t0.vcmax.deactivation_energy,
                entropy = t0.vcmax.entropy,
            )
            t_new = C3TemperatureResponse(;
                vcmax = vc_law, jmax = t0.jmax, respiration = t0.respiration,
                kc = t0.kc, ko = t0.ko, gamma_star = t0.gamma_star,
            )
            m_new = C3Photosynthesis(;
                temperature = t_new,
                electron_transport = base_leaf.photosynthesis.electron_transport,
                stoichiometry = base_leaf.photosynthesis.stoichiometry,
                limitation = base_leaf.photosynthesis.limitation,
                compensation = base_leaf.photosynthesis.compensation,
            )
            photosynthesize(m_new, site_arg).net_assimilation
        end
        d_vc = ForwardDiff.derivative(vc -> fn_vcmax(vc, c3_prog, site_fixed), 60.0)
        d_vc_ref = _bigfloat_deriv(vc -> fn_vcmax(vc, c3_prog_bf, site_fixed_bf), 60.0)
        @test isapprox(d_vc, d_vc_ref; rtol = 1e-11)

        fn_cyto_b6f = (b6f, base_photo, site_arg) -> begin
            t0 = base_photo.temperature
            vq0 = t0.vqmax
            vq_new = CytochromeCapacity(b6f, vq0.turnover_rate)
            t_new = CytochromeTemperatureResponse(;
                vcmax = t0.vcmax, vqmax = vq_new, respiration = t0.respiration,
                kc = t0.kc, ko = t0.ko, gamma_star = t0.gamma_star,
                linear_atp_per_electron = t0.linear_atp_per_electron,
                cyclic_atp_per_electron = t0.cyclic_atp_per_electron,
            )
            m_new = CytochromePhotosynthesis(;
                temperature = t_new,
                electron_transport = base_photo.electron_transport,
                stoichiometry = base_photo.stoichiometry,
                limitation = base_photo.limitation,
                compensation = base_photo.compensation,
            )
            photosynthesize(m_new, LeafEnvironment(site_arg.temperature, eltype(site_arg.temperature)(200), site_arg.co2, site_arg.oxygen)).net_assimilation
        end
        d_b6f = ForwardDiff.derivative(b -> fn_cyto_b6f(b, jb64.photo, site_fixed), 1.4)
        d_b6f_ref = _bigfloat_deriv(b -> fn_cyto_b6f(b, jb_bf.photo, site_fixed_bf), 1.4)
        @test d_b6f > 0 && isapprox(d_b6f, d_b6f_ref; rtol = 1e-11)

        d_med_g1 = ForwardDiff.derivative(
            g1 -> equilibrium_gas_exchange(
                c3_prog.photosynthesis,
                SteadyStateStomata(Medlyn(0.01, g1, 1.6, 1e-6, 50.0)),
                GasExchangeEnvironment(x_eq...),
                diff64,
                coarse_newton,
            ).stomatal_conductance,
            4.0,
        )
        d_med_g1_ref = _bigfloat_deriv(
            g1 -> equilibrium_gas_exchange(
                c3_prog_bf.photosynthesis,
                SteadyStateStomata(Medlyn(BigFloat("0.01"), g1, BigFloat("1.6"), BigFloat("1e-6"), BigFloat(50))),
                GasExchangeEnvironment(BigFloat.(x_eq)...),
                diff_bf,
                tight_bf_solver,
            ).stomatal_conductance,
            4.0,
        )
        @test d_med_g1 > 0 && isapprox(d_med_g1, d_med_g1_ref; rtol = 1e-10)

        # PrognosticStomata tendency AD at dark/zero uptake
        drivers_zero = StomatalDrivers(0.0, 40.0, 101325.0, 0.65, 1200.0)
        d_tau = ForwardDiff.derivative(
            tau -> stomatal_tendency(
                PrognosticStomata(Medlyn(0.01, 4.0, 1.6, 1e-6, 50.0), 300.0, tau),
                StomatalState(0.12),
                drivers_zero,
            ),
            600.0,
        )
        @test isfinite(d_tau) && isapprox(d_tau, (0.12 - 0.01) / 600.0^2; rtol = 1e-12)

        # 5. Boundary & second-order (Hessian) AD smoothness checks
        @test ForwardDiff.gradient(x -> colimit(x[1], x[2], QuadraticLimitation(1.0)), [10.0, 10.0]) == [0.5, 0.5]
        @test LeafPhotosynthesis.colimit_slope(10.0, 10.0, 1.0, 0.0, QuadraticLimitation(1.0)) == (10.0, 1.0)

        # QuadraticLimitation at exact θ == 1.0 with distinct branches (a != b):
        # ∂r/∂θ and ∂(colimit_slope)/∂θ must match the smooth quadratic limit (not zero)
        d_theta_colimit = ForwardDiff.derivative(θ -> colimit(8.0, 12.0, QuadraticLimitation(θ)), 1.0)
        d_theta_colimit_ref = _bigfloat_left_deriv(θ -> colimit(BigFloat(8), BigFloat(12), QuadraticLimitation(θ)), 1.0)
        @test isapprox(d_theta_colimit, 16.0; rtol = 1e-13)
        @test isapprox(d_theta_colimit, d_theta_colimit_ref; rtol = 1e-11)
        d_theta_slope = ForwardDiff.derivative(
            θ -> LeafPhotosynthesis.colimit_slope(8.0, 12.0, 1.5, -0.5, QuadraticLimitation(θ))[2],
            1.0,
        )
        d_theta_slope_ref = _bigfloat_left_deriv(
            θ -> LeafPhotosynthesis.colimit_slope(
                BigFloat(8), BigFloat(12), BigFloat("1.5"), BigFloat("-0.5"), QuadraticLimitation(θ),
            )[2],
            1.0,
        )
        @test isapprox(d_theta_slope, 8 * 1.5 - 4 * (-0.5); rtol = 1e-13)
        @test isapprox(d_theta_slope, d_theta_slope_ref; rtol = 1e-11)

        # Compensation-point boundary assimilation == 0.0 in stomatal target_conductance
        for (st64, st_bf) in ((med64.target, med_bf.target), (bb64.target, bb_bf.target))
            dg_da = ForwardDiff.derivative(
                a -> target_conductance(st64, StomatalDrivers(a, 38.0, 101325.0, 0.65, 1200.0)),
                0.0,
            )
            dg_da_ref = _bigfloat_right_deriv(
                a -> target_conductance(
                    st_bf,
                    StomatalDrivers(a, BigFloat(38), BigFloat(101325), BigFloat("0.65"), BigFloat(1200)),
                ),
                0.0,
            )
            @test dg_da > 0
            @test isapprox(dg_da, dg_da_ref; rtol = 1e-12)
        end

        # Exact darkness (absorbed_par == 0.0): physical initial light-response slopes
        # for BOTH gas_exchange and equilibrium_gas_exchange across C3 and Cytochrome
        gs_dark64 = StomatalState(0.08)
        gs_dark_bf = StomatalState(BigFloat("0.08"))
        for (p_m64, p_mbf) in ((c3_prog.photosynthesis, c3_prog_bf.photosynthesis), (jb64.photo, jb_bf.photo))
            ge_fn64 = par -> gas_exchange(
                p_m64,
                GasExchangeEnvironment(298.15, par, 40.0, 21000.0, 101325.0, 0.65, 1200.0),
                gs_dark64,
                diff64,
            )
            ge_fnbf = par -> gas_exchange(
                p_mbf,
                GasExchangeEnvironment(
                    BigFloat("298.15"), par, BigFloat(40), BigFloat(21000),
                    BigFloat(101325), BigFloat("0.65"), BigFloat(1200),
                ),
                gs_dark_bf,
                diff_bf,
            )
            dAn_ge = ForwardDiff.derivative(par -> ge_fn64(par).carbon.net_assimilation, 0.0)
            dCc_ge = ForwardDiff.derivative(par -> ge_fn64(par).biochemical_co2, 0.0)
            dCi_ge = ForwardDiff.derivative(par -> ge_fn64(par).intercellular_co2, 0.0)
            dRes_ge = ForwardDiff.derivative(
                par -> begin
                    ex = ge_fn64(par)
                    supply = ex.conductance * (diff64.mole_to_micromole / 101325.0) * (40.0 - ex.biochemical_co2)
                    ex.carbon.net_assimilation - supply
                end,
                0.0,
            )
            @test dAn_ge > 0 && isapprox(dAn_ge, _bigfloat_right_deriv(par -> ge_fnbf(par).carbon.net_assimilation, 0.0); rtol = 1e-10)
            @test dCc_ge < 0 && isapprox(dCc_ge, _bigfloat_right_deriv(par -> ge_fnbf(par).biochemical_co2, 0.0); rtol = 1e-10)
            @test dCi_ge < 0 && isapprox(dCi_ge, _bigfloat_right_deriv(par -> ge_fnbf(par).intercellular_co2, 0.0); rtol = 1e-10)
            @test abs(dRes_ge) <= 1e-12

            for slv in (coarse_newton, bisect_solver)
                eq_fn64 = par -> equilibrium_gas_exchange(
                    p_m64,
                    med64,
                    GasExchangeEnvironment(298.15, par, 40.0, 21000.0, 101325.0, 0.65, 1200.0),
                    diff64,
                    slv,
                )
                eq_fnbf = par -> equilibrium_gas_exchange(
                    p_mbf,
                    med_bf,
                    GasExchangeEnvironment(
                        BigFloat("298.15"), par, BigFloat(40), BigFloat(21000),
                        BigFloat(101325), BigFloat("0.65"), BigFloat(1200),
                    ),
                    diff_bf,
                    tight_bf_solver,
                )
                dAn_eq = ForwardDiff.derivative(par -> eq_fn64(par).carbon.net_assimilation, 0.0)
                dCc_eq = ForwardDiff.derivative(par -> eq_fn64(par).biochemical_co2, 0.0)
                dCi_eq = ForwardDiff.derivative(par -> eq_fn64(par).intercellular_co2, 0.0)
                dRes_eq = ForwardDiff.derivative(par -> eq_fn64(par).diffusion_residual, 0.0)
                @test dAn_eq > 0 && isapprox(dAn_eq, _bigfloat_right_deriv(par -> eq_fnbf(par).carbon.net_assimilation, 0.0); rtol = 1e-10)
                @test dCc_eq < 0 && isapprox(dCc_eq, _bigfloat_right_deriv(par -> eq_fnbf(par).biochemical_co2, 0.0); rtol = 1e-10)
                @test dCi_eq < 0 && isapprox(dCi_eq, _bigfloat_right_deriv(par -> eq_fnbf(par).intercellular_co2, 0.0); rtol = 1e-10)
                @test abs(dRes_eq) <= 1e-12
            end
        end

        # Zero biochemical capacity (Vcmax == 0.0 and Vqmax == 0.0) right-sided sensitivities
        vc0_eq64 = vc -> begin
            t0 = c3_prog.photosynthesis.temperature
            a0 = t0.vcmax.arrhenius
            vc_law = PeakedArrhenius(
                Arrhenius(vc, a0.activation_energy, a0.temperature_ref, a0.gas_constant),
                t0.vcmax.deactivation_energy, t0.vcmax.entropy,
            )
            m = C3Photosynthesis(;
                temperature = C3TemperatureResponse(vc_law, t0.jmax, t0.respiration, t0.kc, t0.ko, t0.gamma_star),
                electron_transport = c3_prog.photosynthesis.electron_transport,
                stoichiometry = c3_prog.photosynthesis.stoichiometry,
                limitation = c3_prog.photosynthesis.limitation,
                compensation = c3_prog.photosynthesis.compensation,
            )
            equilibrium_gas_exchange(
                m, med64,
                GasExchangeEnvironment(298.15, 800.0, 40.0, 21000.0, 101325.0, 0.65, 1200.0),
                diff64, coarse_newton,
            )
        end
        vc0_eqbf = vc -> begin
            t0 = c3_prog_bf.photosynthesis.temperature
            a0 = t0.vcmax.arrhenius
            vc_law = PeakedArrhenius(
                Arrhenius(vc, a0.activation_energy, a0.temperature_ref, a0.gas_constant),
                t0.vcmax.deactivation_energy, t0.vcmax.entropy,
            )
            m = C3Photosynthesis(;
                temperature = C3TemperatureResponse(vc_law, t0.jmax, t0.respiration, t0.kc, t0.ko, t0.gamma_star),
                electron_transport = c3_prog_bf.photosynthesis.electron_transport,
                stoichiometry = c3_prog_bf.photosynthesis.stoichiometry,
                limitation = c3_prog_bf.photosynthesis.limitation,
                compensation = c3_prog_bf.photosynthesis.compensation,
            )
            equilibrium_gas_exchange(
                m, med_bf,
                GasExchangeEnvironment(
                    BigFloat("298.15"), BigFloat(800), BigFloat(40), BigFloat(21000),
                    BigFloat(101325), BigFloat("0.65"), BigFloat(1200),
                ),
                diff_bf, tight_bf_solver,
            )
        end
        dAn_vc0 = ForwardDiff.derivative(vc -> vc0_eq64(vc).carbon.net_assimilation, 0.0)
        dCc_vc0 = ForwardDiff.derivative(vc -> vc0_eq64(vc).biochemical_co2, 0.0)
        @test dAn_vc0 > 0 && isapprox(dAn_vc0, _bigfloat_right_deriv(vc -> vc0_eqbf(vc).carbon.net_assimilation, 0.0); rtol = 1e-10)
        @test dCc_vc0 < 0 && isapprox(dCc_vc0, _bigfloat_right_deriv(vc -> vc0_eqbf(vc).biochemical_co2, 0.0); rtol = 1e-10)

        dq2_b6f0 = ForwardDiff.derivative(
            b -> cytochrome_electron_supply(
                jb64.photo.electron_transport,
                400.0,
                temperature_response(CytochromeCapacity(b, jb64.photo.temperature.vqmax.turnover_rate), 298.15),
            ).psii_open_fraction,
            0.0,
        )
        dq2_b6f0_ref = _bigfloat_right_deriv(
            b -> cytochrome_electron_supply(
                jb_bf.photo.electron_transport,
                BigFloat(400),
                temperature_response(CytochromeCapacity(b, jb_bf.photo.temperature.vqmax.turnover_rate), BigFloat("298.15")),
            ).psii_open_fraction,
            0.0,
        )
        @test dq2_b6f0 > 0 && isapprox(dq2_b6f0, dq2_b6f0_ref; rtol = 1e-11)

        dq2_par0 = ForwardDiff.derivative(
            par -> cytochrome_electron_supply(jb64.photo.electron_transport, par, 120.0).psii_open_fraction,
            0.0,
        )
        dq2_par0_ref = _bigfloat_right_deriv(
            par -> cytochrome_electron_supply(jb_bf.photo.electron_transport, par, BigFloat(120)).psii_open_fraction,
            0.0,
        )
        @test dq2_par0 < 0 && isapprox(dq2_par0, dq2_par0_ref; rtol = 1e-11)

        # Second-order ForwardDiff.hessian accuracy on equilibrium_gas_exchange under coarse solver tolerances
        very_coarse_newton = EquilibriumSolver(;
            method = SafeguardedNewton(4),
            co2_atol = 1e-3, co2_rtol = 0.0, flux_atol = 1e-3, flux_rtol = 0.0,
            max_iterations = 40,
        )
        tight_f64_newton = EquilibriumSolver(;
            method = SafeguardedNewton(4),
            co2_atol = 1e-12, co2_rtol = 0.0, flux_atol = 1e-12, flux_rtol = 0.0,
            max_iterations = 80,
        )
        for p_hess in (c3_prog.photosynthesis, jb64.photo)
            h_coarse = ForwardDiff.hessian(
                v -> equilibrium_gas_exchange(
                    p_hess, med64,
                    GasExchangeEnvironment(v[1], v[2], v[3], 21000.0, 101325.0, 0.62, 1350.0),
                    diff64, very_coarse_newton,
                ).carbon.net_assimilation,
                [298.15, 450.0, 40.0],
            )
            h_tight = ForwardDiff.hessian(
                v -> equilibrium_gas_exchange(
                    p_hess, med64,
                    GasExchangeEnvironment(v[1], v[2], v[3], 21000.0, 101325.0, 0.62, 1350.0),
                    diff64, tight_f64_newton,
                ).carbon.net_assimilation,
                [298.15, 450.0, 40.0],
            )
            @test all(isfinite, h_coarse)
            @test isapprox(h_coarse, h_tight; rtol = 1e-10, atol = 1e-12)
        end

        # PeakedArrhenius at exact z == 0 (T == hd / ds): opposite +/- directional partials & smooth 2nd derivative
        pa = PeakedArrhenius(Arrhenius(60.0, 65000.0, 298.15, 8.314), 200000.0, 650.0)
        T_crit = 200000.0 / 650.0
        y_plus  = temperature_response(pa, ForwardDiff.Dual(T_crit, +1.0))
        y_minus = temperature_response(pa, ForwardDiff.Dual(T_crit, -1.0))
        @test ForwardDiff.value(y_plus) == ForwardDiff.value(y_minus)
        @test isapprox(ForwardDiff.partials(y_plus)[1], -ForwardDiff.partials(y_minus)[1]; rtol = 1e-14)
        h_pa = ForwardDiff.hessian(v -> temperature_response(pa, v[1]), [T_crit])
        @test isfinite(h_pa[1, 1])

        # SaturationQuenching at saturation == 0 across exponents < 1, == 1, > 1
        for expo in (1.0, 2.83)
            sq = SaturationQuenching(2.5, expo, 1.0)
            dq0 = ForwardDiff.derivative(s -> quenching_rate(sq, s), 0.0)
            @test isfinite(dq0) && isapprox(dq0, expo == 1.0 ? 5.0 : 0.0; rtol = 1e-14)
        end
        # Sublinear exponent < 1 at dark zero-partial Dual(0.0, 0.0) & Dual(-0.0, 0.0) avoids 0 * Inf = NaN
        @test ForwardDiff.derivative(kmax -> quenching_rate(SaturationQuenching(kmax, 0.5, 1.0), 0.0 * kmax), 2.5) === 0.0
        @test quenching_rate(SaturationQuenching(2.5, 0.5, 1.0), ForwardDiff.Dual(-0.0, 0.0)) == 0.0

        # diffusion_co2_root(m, k, gamma, rd, ca, g) at gamma = rd = ca = 0 (d == 0)
        g_diff0 = ForwardDiff.gradient(
            v -> LeafPhotosynthesis.diffusion_co2_root(v[1], v[2], v[3], v[4], v[5], v[6]),
            [20.0, 40.0, 0.0, 0.0, 0.0, 0.5],
        )
        @test all(isfinite, g_diff0)
        @test isapprox(g_diff0, [0.0, 0.0, 0.5, 1.0, 0.5, 0.0]; rtol = 1e-13)
    end
    end
end
