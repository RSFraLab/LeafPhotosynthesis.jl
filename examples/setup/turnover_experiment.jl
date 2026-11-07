# Keep loading and experiment assembly out of the biophysical workflow.
isdefined(@__MODULE__, :core_leaf) || include(joinpath(@__DIR__, "..", "core_exchange.jl"))

function turnover_experiment(::Type{FT}) where {FT}
    directory = joinpath(@__DIR__, "..", "parameters")
    biochemistry = read_bonan_parameters(FT, joinpath(directory, "bonan_11_1.toml"))
    settings = read_bonan_parameters(FT, joinpath(directory, "core_exchange.toml"))
    sweep = read_bonan_parameters(FT, joinpath(directory, "turnover.toml"))["experiment"]
    base, _ = core_leaf(FT, biochemistry, settings)
    # Turnover does not require the optional oxygen or stomatal calculations.
    leaf = LeafModel(; photosynthesis = base.photosynthesis,
        oxygen = NoOxygenExchange(), stomata = NoStomatalDynamics())
    forcing = core_environment(FT, biochemistry, settings)
    rates = biochemical_rates(leaf.photosynthesis.temperature, forcing.temperature)
    ratios = collect(sweep["co2_gamma_ratio_start"]:sweep["co2_gamma_ratio_step"]:sweep["co2_gamma_ratio_stop"])
    # Include the compensation point exactly, even if a changed grid misses it.
    sort!(unique!(push!(ratios, one(FT))))
    (;leaf,forcing,gamma_star=rates.gamma_star,ratios)
end

function write_turnover_example(rows,::Type{FT},directory) where {FT}
    write_tables(joinpath(directory,string(FT)),Dict("turnover"=>rows))
end
