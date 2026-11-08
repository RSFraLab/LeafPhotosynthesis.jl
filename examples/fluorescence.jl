# Shared PSII excitation drives carbon uptake and a replaceable fluorescence source.
include("setup/fluorescence.jl")

function fluorescence_response(::Type{FT}) where {FT}
    case = setup_fluorescence(FT)
    map(case.environments) do site
        carbon = photosynthesize(case.leaf, site)
        inputs = psii_flux_inputs(case.leaf.photosynthesis, site)
        unquenched = fluorescence(case.unquenched, inputs)
        regulated = fluorescence(case.regulated, inputs)
        fluorescence_row(site, carbon, inputs, unquenched, regulated)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    directory = isempty(ARGS) ? joinpath(@__DIR__, "output", "fluorescence") : only(ARGS)
    for FT in (Float32, Float64)
        write_fluorescence(fluorescence_response(FT), FT, directory)
    end
end
