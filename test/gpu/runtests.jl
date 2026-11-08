using Test
using CUDA
using LeafPhotosynthesis
include("../../examples/core_exchange.jl")
include("../device_checks.jl")
include("../../examples/setup/oxygen_response.jl")
include("../oxygen_batches.jl")
include("../../examples/setup/fluorescence.jl")
include("../fluorescence_batches.jl")

# Do not silently skip this explicitly requested GPU acceptance test.
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
@testset "CUDA leaf physics" begin
    path = joinpath(@__DIR__, "..", "..", "examples", "parameters", "bonan_11_1.toml")
    check_leaf_batches(CuArray, CUDA.synchronize, path)
    check_core_batches(CuArray, CUDA.synchronize, path)
    check_turnover_batches(CuArray, CUDA.synchronize, path)
    check_equilibrium_batches(CuArray, CUDA.synchronize, path)
    check_oxygen_batches(CuArray, CUDA.synchronize)
    check_fluorescence_batches(CuArray, CUDA.synchronize)
end
include("../excitation_batches.jl")
@testset "Explicit excitation CUDA batches" begin
    check_excitation_batches(CuArray,CUDA.synchronize)
end
