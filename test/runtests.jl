using Test
using LeafPhotosynthesis
include("../examples/core_exchange.jl")
include("reference_bonan.jl")
include("device_checks.jl")

parameter_path = joinpath(@__DIR__, "..", "examples", "parameters", "bonan_11_1.toml")
reference_parameters = read_bonan_parameters(BigFloat, parameter_path)

# Keep allocation measurement inside a specialization, outside global scope.
function leaf_allocations(model, environment)
    photosynthesize(model, environment)
    return @allocated photosynthesize(model, environment)
end

include("core_exchange.jl")
include("turnover.jl")
include("equilibrium.jl")
include("../examples/setup/oxygen_response.jl")
include("oxygen_compensation.jl")
include("oxygen_batches.jl")
include("../examples/setup/fluorescence.jl")
include("fluorescence.jl")
include("fluorescence_batches.jl")
@testset "Fluorescence CPU batches" begin
    check_fluorescence_batches(identity, () -> nothing)
end
@testset "Oxygen-aware CPU batches" begin
    check_oxygen_batches(identity, () -> nothing)
end

@testset "LeafPhotosynthesis" begin
    for FT in (Float32, Float64)
        @testset "$FT" begin
            p = read_bonan_parameters(FT, parameter_path)
            model = bonan_leaf(FT, p)
            e, c = p["experiment"], p["constants"]
            env = bonan_environment(FT, p; incident_par = e["incident_par"],
                internal_co2 = e["internal_co2_fraction"] * e["atmospheric_co2"],
                temperature_celsius = e["temperature_celsius"])
            @test isbitstype(typeof(model))
            @test isbitstype(typeof(env))
            result = @inferred photosynthesize(model, env)
            @test all(x -> x isa FT, values(result))
            @test leaf_allocations(model, env) == 0

            @testset "Reference normalization" begin
                for law in (model.photosynthesis.temperature.vcmax,
                    model.photosynthesis.temperature.jmax,
                    model.photosynthesis.temperature.respiration,
                    model.photosynthesis.temperature.kc,
                    model.photosynthesis.temperature.ko,
                    model.photosynthesis.temperature.gamma_star)
                    base = law isa PeakedArrhenius ? law.arrhenius : law
                    @test (@inferred temperature_response(law, base.temperature_ref)) ≈ base.value_ref
                    @test temperature_response(law, base.temperature_ref) isa FT
                end
            end

            @testset "Darkness and compensation point" begin
                r = biochemical_rates(model.photosynthesis.temperature, env.temperature)
                dark = LeafEnvironment{FT}(env.temperature, zero(FT), env.co2, env.oxygen)
                d = @inferred photosynthesize(model, dark)
                @test d.electron_transport === zero(FT)
                @test d.net_assimilation == -r.respiration
                at_compensation = LeafEnvironment{FT}(env.temperature, env.absorbed_par,
                    r.gamma_star, env.oxygen)
                @test photosynthesize(model, at_compensation).net_assimilation == -r.respiration
                low = LeafEnvironment{FT}(env.temperature, env.absorbed_par, zero(FT), env.oxygen)
                @test photosynthesize(model, low).net_assimilation == -r.respiration
                s = model.photosynthesis
                signed_model = LeafModel(; photosynthesis = C3Photosynthesis(;
                    temperature = s.temperature, electron_transport = s.electron_transport,
                    stoichiometry = s.stoichiometry, limitation = MinimumLimitation(),
                    compensation = SignedAssimilation()))
                @test photosynthesize(signed_model, low).net_assimilation < -r.respiration
                @test all(x -> x isa FT, values(@inferred photosynthesize(signed_model, low)))
                # Dark endpoints do not suppress interior illuminated leaves.
                population = [dark, env, dark]
                batch = photosynthesize.(Ref(model), population)
                @test batch[2] == result
                @test batch[[2, 3, 1]] == photosynthesize.(Ref(model), population[[2, 3, 1]])
                @test isempty(photosynthesize.(Ref(model), LeafEnvironment{FT}[]))
            end

            @testset "Stable quadratic and physical limits" begin
                theta = p["electron_transport"]["curvature"]
                quadratic = QuadraticLimitation(theta)
                @test colimit(zero(FT), zero(FT), quadratic) === zero(FT)
                for a in FT.((0, 1e-20, 1, 1e20)), b in FT.((0, 1e-20, 1, 1e20))
                    value = @inferred colimit(a, b, quadratic)
                    @test value isa FT
                    @test isfinite(value)
                    @test zero(FT) <= value <= min(a, b) * (one(FT) + 8eps(FT))
                    @test value ≈ colimit(b, a, quadratic)
                    # Polynomial residual evaluated in high precision and normalized by scale.
                    scale = max(BigFloat(a), BigFloat(b), BigFloat(1))
                    residual = theta * (BigFloat(value) / scale)^2 -
                        ((BigFloat(a) + b) / scale) * (BigFloat(value) / scale) +
                        (BigFloat(a) / scale) * (BigFloat(b) / scale)
                    @test abs(residual) <= 16eps(FT)
                end
                for a in FT.((-2, 0, 2)), b in FT.((-2, 0, 2))
                    @test colimit(a, b, QuadraticLimitation(one(FT))) ≈ min(a, b)
                end
                @test_throws ArgumentError QuadraticLimitation(zero(FT))
                @test_throws ArgumentError QuadraticLimitation(FT(NaN))
                # Ratio underflows in FP32, but the smaller rate is representable.
                @test colimit(FT(1e-30), FT(1e30), quadratic) ≈ FT(1e-30)
            end

            @testset "All Bonan 11.1 response points" begin
                curves = bonan_curves(FT, p)
                @test length(curves["light"]) == 201
                @test length(curves["co2"]) == 97
                @test length(curves["temperature"]) == 26
                for rows in values(curves), row in rows
                    @test all(x -> x isa FT, values(row))
                    oracle = reference_bonan(reference_parameters, row.incident_par,
                        row.internal_co2, row.temperature_celsius)
                    for key in keys(oracle)
                        @test isapprox(row[key], oracle[key]; rtol = 128eps(FT),
                            atol = 32eps(FT) * p["biochemistry"]["vcmax_ref"])
                    end
                end
                mktempdir() do directory
                    write_tables(directory, curves)
                    @test length(readdir(directory)) == length(curves)
                    for name in keys(curves)
                        @test countlines(joinpath(directory, name * ".csv")) == length(curves[name]) + 1
                    end
                end
            end
        end
    end
    check_leaf_batches(identity, () -> nothing, parameter_path)
    check_core_batches(identity, () -> nothing, parameter_path)
    check_turnover_batches(identity, () -> nothing, parameter_path)
    check_equilibrium_batches(identity, () -> nothing, parameter_path)
end
