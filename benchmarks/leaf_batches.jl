# Run on an idle GPU, using the optional CUDA environment:
# CUDA_VISIBLE_DEVICES=1 julia --threads=8 --project=test/gpu benchmarks/leaf_batches.jl OUTPUT.csv
using CUDA
using TOML
include("../examples/bonan_11_1.jl")

function threaded_leaf!(output, leaf, input)
    Threads.@threads :static for i in eachindex(input)
        @inbounds output[i] = photosynthesize(leaf, input[i])
    end
    return output
end

function measure_batch(leaf, pattern, n, timing)
    # Repeated heterogeneous reference environments exercise temperature/light/CO₂
    # branches. This is a kernel workload, not an Earth-system simulation.
    output_type = typeof(photosynthesize(leaf, first(pattern)))
    input = Vector{eltype(pattern)}(undef, n)
    output = Vector{output_type}(undef, n)
    for i in eachindex(input)
        @inbounds input[i] = pattern[mod1(i, length(pattern))]
    end
    warm_input, warm_output = copy(pattern), Vector{output_type}(undef, length(pattern))
    broadcast!(photosynthesize, warm_output, Ref(leaf), warm_input)
    threaded_leaf!(warm_output, leaf, warm_input)
    serial = minimum([@elapsed broadcast!(photosynthesize, output, Ref(leaf), input)
        for _ in 1:timing["cpu_samples"]])
    threaded = minimum([@elapsed threaded_leaf!(output, leaf, input)
        for _ in 1:timing["cpu_samples"]])
    indices = unique(round.(Int, range(firstindex(input), lastindex(input); length = timing["check_samples"])))
    expected = output[indices]
    device_input = CuArray{eltype(input)}(undef, n)
    device_output = CuArray{output_type}(undef, n)
    upload = @elapsed begin
        copyto!(device_input, input)
        CUDA.synchronize()
    end
    # Warm-up excludes compilation from event timing.
    broadcast!(photosynthesize, device_output, Ref(leaf), device_input)
    CUDA.synchronize()
    gpu = minimum([CUDA.@elapsed broadcast!(photosynthesize, device_output, Ref(leaf), device_input)
        for _ in 1:timing["gpu_samples"]])
    download = @elapsed begin
        copyto!(output, device_output)
        CUDA.synchronize()
    end
    FT = typeof(first(pattern).temperature)
    tolerance = FT(timing["check_tolerance_eps"]) * eps(FT)
    for (got, want) in zip(output[indices], expected), key in keys(want)
        isapprox(got[key], want[key]; atol = tolerance, rtol = tolerance) ||
            error("CPU/GPU benchmark sample mismatch at $key")
    end
    payload_bytes = n * (sizeof(eltype(input)) + sizeof(output_type))
    CUDA.unsafe_free!(device_input)
    CUDA.unsafe_free!(device_output)
    return (; precision = string(FT), leaves = n, payload_bytes, cpu_threads = Threads.nthreads(),
        cpu_serial_seconds = serial, cpu_threaded_seconds = threaded,
        gpu_seconds = gpu, upload_seconds = upload, download_seconds = download,
        serial_speedup = serial / gpu, threaded_speedup = threaded / gpu)
end

function main(output_path)
    CUDA.functional(true) || error("A functional GPU is required")
    CUDA.allowscalar(false)
    settings = TOML.parsefile(joinpath(@__DIR__, "workload.toml"))
    w = settings["workload"]
    full = round(Int, w["latitudes"] * w["longitudes"] * w["land_fraction"] * w["leaves_per_column"])
    sizes = [w["smaller_batches"]; full]
    println("GPU: ", CUDA.name(CUDA.device()), "; Julia ", VERSION, "; CPU threads: ", Threads.nthreads())
    mkpath(dirname(abspath(output_path)))
    open(output_path, "w") do io
        header_written = false
        for FT in (Float32, Float64)
            p = read_bonan_parameters(FT, joinpath(@__DIR__, "..", "examples", "parameters", "bonan_11_1.toml"))
            leaf = bonan_leaf(FT, p)
            curves = bonan_curves(FT, p)
            rows = vcat((curves[key] for key in ("light", "co2", "temperature"))...)
            pattern = [bonan_environment(FT, p; incident_par = row.incident_par,
                internal_co2 = row.internal_co2, temperature_celsius = row.temperature_celsius) for row in rows]
            for n in sizes
                result = measure_batch(leaf, pattern, n, settings["timing"])
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
