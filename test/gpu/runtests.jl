# CUDA acceptance test: the CPU batch checks of the main suite, run on CuArrays. Not run in CI.
#   julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
#   CUDA_VISIBLE_DEVICES=<idle device> julia --project=test/gpu test/gpu/runtests.jl
using CUDA
include("../common.jl")   # Test, LeafPhotosynthesis and the example fixtures, as in every group
include("../device_checks.jl")
include("../oxygen_batches.jl")
include("../fluorescence_batches.jl")
include("../excitation_fixture.jl")
include("../excitation_batches.jl")
include("../cytochrome_batches.jl")
include("../closed_stomata_batches.jl")

# Do not silently skip this explicitly requested GPU acceptance test.
CUDA.functional(true) || error("A functional CUDA device is required")
CUDA.allowscalar(false)
@testset "CUDA leaf physics" begin
    check_leaf_batches(CuArray, CUDA.synchronize, parameter_path)
    check_core_batches(CuArray, CUDA.synchronize, parameter_path)
    check_turnover_batches(CuArray, CUDA.synchronize, parameter_path)
    check_equilibrium_batches(CuArray, CUDA.synchronize, parameter_path)
    check_oxygen_batches(CuArray, CUDA.synchronize)
    check_fluorescence_batches(CuArray, CUDA.synchronize)
end
@testset "Explicit excitation CUDA batches" begin
    check_excitation_batches(CuArray, CUDA.synchronize)
end
@testset "Cytochrome CUDA batches" begin
    check_cytochrome_batches(CuArray, CUDA.synchronize)
end
@testset "Closed-stomata CUDA batches" begin
    check_closed_stomata_batches(CuArray, CUDA.synchronize)
end
