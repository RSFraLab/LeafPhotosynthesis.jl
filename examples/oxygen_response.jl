# Change O2 while keeping its reference compensation-point calibration explicit.
include("setup/oxygen_response.jl")

function oxygen_response(::Type{FT}) where {FT}
    case = setup_oxygen_response(FT)
    map(oxygen_response_environments(FT, case)) do environment
        responses = map((case.prescribed, case.scaled)) do leaf
            steady = LeafModel(;
                leaf.photosynthesis,
                stomata = SteadyStateStomata(leaf.stomata.target),
                leaf.oxygen,
            )
            exchange = equilibrium_gas_exchange(
                steady,
                environment,
                case.diffusion,
                case.solver,
            )
            exchange.converged || error("Leaf equilibrium failed: $(exchange.status)")
            rates = biochemical_rates(
                leaf.photosynthesis.temperature,
                environment.temperature,
                environment.oxygen,
            )
            site = LeafEnvironment(
                environment.temperature,
                environment.absorbed_par,
                exchange.biochemical_co2,
                environment.oxygen,
            )
            oxygen = oxygen_exchange(leaf, site)
            (;
                gamma_star = rates.gamma_star,
                assimilation = exchange.carbon.net_assimilation,
                conductance = exchange.stomatal_conductance,
                net_oxygen = oxygen.net_release,
            )
        end
        prescribed, scaled = responses
        (;
            temperature = environment.temperature,
            oxygen = environment.oxygen,
            prescribed_gamma = prescribed.gamma_star,
            scaled_gamma = scaled.gamma_star,
            prescribed_assimilation = prescribed.assimilation,
            scaled_assimilation = scaled.assimilation,
            prescribed_conductance = prescribed.conductance,
            scaled_conductance = scaled.conductance,
            prescribed_net_oxygen = prescribed.net_oxygen,
            scaled_net_oxygen = scaled.net_oxygen,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "oxygen_response") : only(ARGS)
    for FT in (Float32, Float64)
        write_oxygen_response(oxygen_response(FT), FT, directory)
    end
end
