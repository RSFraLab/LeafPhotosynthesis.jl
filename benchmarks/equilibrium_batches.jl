# CUDA_VISIBLE_DEVICES=1 julia --threads=8 --project=test/gpu benchmarks/equilibrium_batches.jl OUTPUT.csv
using CUDA
using TOML
include("../examples/core_exchange.jl")

function threaded_equilibrium!(output, leaf, input, diffusion, solver)
    Threads.@threads :static for i in eachindex(input)
        @inbounds output[i] = equilibrium_gas_exchange(leaf, input[i], diffusion, solver)
    end
    return output
end

function measure_equilibrium(leaf, pattern, d, solver, settings)
    n = settings["batch_size"]
    input = [pattern[mod1(i, length(pattern))] for i in 1:n]
    output = [equilibrium_gas_exchange(leaf, first(pattern), d, solver) for _ in 1:n]
    threaded_equilibrium!(output, leaf, input, d, solver)
    all(r -> r.converged, output) || error("CPU equilibrium failed")
    cpu = minimum([@elapsed threaded_equilibrium!(output, leaf, input, d, solver)
        for _ in 1:settings["cpu_samples"]])
    cpu_iterations = sum(r.iterations for r in output) / n
    cpu_max_iterations = maximum(r.iterations for r in output)
    indices = unique(round.(Int, range(1, n; length = length(pattern))))
    expected = output[indices]
    device_input = CuArray{eltype(input)}(undef, n)
    device_output = CuArray{eltype(output)}(undef, n)
    # Warm transfer dispatch as well as the compute kernel; first-call compilation
    # would otherwise contaminate the first precision/method's transfer timings.
    copyto!(device_input, input)
    CUDA.synchronize()
    upload = @elapsed begin
        copyto!(device_input, input)
        CUDA.synchronize()
    end
    broadcast!(equilibrium_gas_exchange, device_output, Ref(leaf), device_input, Ref(d), Ref(solver))
    CUDA.synchronize()
    copyto!(output, device_output)
    CUDA.synchronize()
    gpu = minimum([CUDA.@elapsed broadcast!(equilibrium_gas_exchange, device_output,
        Ref(leaf), device_input, Ref(d), Ref(solver)) for _ in 1:settings["gpu_samples"]])
    download = @elapsed begin
        copyto!(output, device_output)
        CUDA.synchronize()
    end
    all(r -> r.converged, output) || error("GPU equilibrium failed")
    for (actual, wanted) in zip(output[indices], expected)
        isapprox(actual.carbon.net_assimilation, wanted.carbon.net_assimilation;
            rtol = 8solver.flux_rtol, atol = 8solver.flux_atol) || error("CPU/GPU mismatch")
    end
    gpu_iterations = sum(r.iterations for r in output) / n
    gpu_max_iterations = maximum(r.iterations for r in output)
    CUDA.unsafe_free!(device_input)
    CUDA.unsafe_free!(device_output)
    return (; precision = string(typeof(first(pattern).temperature)), method = string(nameof(typeof(solver.method))),
        leaves = n, cpu_threads = Threads.nthreads(), cpu_seconds = cpu, gpu_seconds = gpu,
        upload_seconds = upload, download_seconds = download, cpu_iterations, cpu_max_iterations,
        gpu_iterations, gpu_max_iterations, threaded_speedup = cpu / gpu)
end

function main(path)
    CUDA.functional(true) || error("A functional CUDA device is required")
    CUDA.allowscalar(false)
    settings = TOML.parsefile(joinpath(@__DIR__, "workload.toml"))["equilibrium"]
    println("GPU: ", CUDA.name(CUDA.device()), "; Julia ", VERSION, "; CPU threads: ", Threads.nthreads())
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        header_written = false
        for FT in (Float32, Float64)
            p = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "examples", "parameters", "bonan_11_1.toml"))
            s = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "examples", "parameters", "core_exchange.toml"))
            base, d = core_leaf(FT, p, s)
            leaf = LeafModel(; photosynthesis = base.photosynthesis, stomata = SteadyStateStomata(base.stomata.target))
            env = core_environment(FT, p, s)
            pattern = [GasExchangeEnvironment(FT(t), FT(light), FT(ca), env.oxygen,
                env.pressure, env.relative_humidity, FT(vpd)) for t in settings["temperature"]
                for light in settings["absorbed_par"] for ca in settings["ambient_co2"] for vpd in settings["vpd"]]
            for method in (SafeguardedNewton(; bisection_interval = Int(s["solver"]["bisection_interval"])), Bisection())
                solver = core_solver(FT, s; method)
                result = measure_equilibrium(leaf, pattern, d, solver, settings)
                if !header_written
                    println(io, join(keys(result), ','))
                    header_written = true
                end
                println(io, join(values(result), ','))
                flush(io)
                println(result)
                flush(stdout)
                GC.gc()
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(only(ARGS))
end
