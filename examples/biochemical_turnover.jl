# Available electrons, realized flow and carbon cancellation near compensation.
include("setup/turnover_experiment.jl")

function turnover_response(::Type{FT}) where {FT}
    experiment = turnover_experiment(FT)
    (;leaf,forcing,gamma_star,ratios) = experiment
    map(ratios) do ratio
        environment = LeafEnvironment(forcing.temperature,forcing.absorbed_par,
            ratio*gamma_star,forcing.oxygen)
        fluxes = biochemical_turnover(leaf,environment)
        carbon = photosynthesize(leaf,environment)
        (;co2_gamma_ratio=ratio,co2=environment.co2,gamma_star,
            fluxes...,net_assimilation=carbon.net_assimilation)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__,"output","turnover") : only(ARGS)
    for FT in (Float32,Float64)
        write_turnover_example(turnover_response(FT),FT,directory)
    end
end
