using Test
using CUDA
using LeafPhotosynthesis
include("../../examples/bonan_11_1.jl")
include("../device_checks.jl")

# Do not silently skip this explicitly requested GPU acceptance test.
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
@testset "CUDA leaf physics" begin
    path = joinpath(@__DIR__, "..", "..", "examples", "parameters", "bonan_11_1.toml")
    check_leaf_batches(CuArray, CUDA.synchronize, path)
end
