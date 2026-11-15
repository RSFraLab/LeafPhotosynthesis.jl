# Test groups. Each group runs in a fresh module that first includes common.jl (imports and
# the example code used as fixtures) and then the group's files in the order listed, so a
# group never depends on another group. With no arguments every group runs, in this order.
#
#   julia --project=. -e 'using Pkg; Pkg.test()'                                # all groups
#   julia --project=. -e 'using Pkg; Pkg.test(test_args = ["cytochrome", "qa"])' # some groups
#   julia --project=test test/runtests.jl cytochrome qa                         # same
#
# The CUDA tests are not a group: they run from gpu/runtests.jl in their own environment.
using Test

const TEST_GROUPS = [
    "core_exchange" => ["core_exchange.jl", "closed_stomata_and_domains.jl"],
    "turnover" => ["turnover.jl"],
    "equilibrium" => ["equilibrium.jl", "closed_stomata_batches.jl", "equilibrium_closed_stomata.jl",
        "stomatal_humidity_and_ratio.jl"],
    "oxygen" => ["oxygen_batches.jl", "oxygen_compensation.jl", "oxygen_budget_stoichiometry.jl",
        "oxygen_inputs_form.jl"],
    "fluorescence" => ["fluorescence_batches.jl", "fluorescence.jl", "quenching_domain.jl"],
    "leaf" => ["reference_bonan.jl", "device_checks.jl", "leaf.jl", "temperature_domains.jl",
        "result_shapes_and_component_types.jl"],
    "excitation" => ["excitation_fixture.jl", "excitation_reference.jl",
        "excitation_batches.jl", "excitation.jl"],
    "cytochrome" => ["cytochrome_batches.jl", "cytochrome.jl", "cytochrome_additions.jl",
        "npq_and_c3_domain.jl"],
    "flexas" => ["flexas_2002.jl"],
    "qa" => ["qa_and_ad.jl", "precision_and_ties.jl"],
    "jet" => ["jet.jl"],
]

function group_module(name, files)
    paths = [joinpath(@__DIR__, file) for file in ["common.jl"; files]]
    body = Expr(:block, (:(include($path)) for path in paths)...)
    return Expr(:module, true, Symbol("TestGroup_", name), body)
end

function run_test_groups(requested)
    available = first.(TEST_GROUPS)
    selected = isempty(requested) || requested == ["all"] ? available : unique(requested)
    unknown = setdiff(selected, available)
    isempty(unknown) ||
        error("Unknown test group(s) $(join(unknown, ", ")); available: $(join(available, ", "))")
    failed = String[]
    for (name, files) in TEST_GROUPS
        name in selected || continue
        @info "Test group \"$name\""
        try
            Core.eval(Main, group_module(name, files))
        catch err
            # A failed top-level testset has printed its report already; show anything else.
            cause = err
            while cause isa LoadError
                cause = cause.error
            end
            if !(cause isa Test.TestSetException)
                showerror(stderr, err, catch_backtrace())
                println(stderr)
            end
            push!(failed, name)
        end
    end
    isempty(failed) || error("Failed test group(s): ", join(failed, ", "))
    return nothing
end

run_test_groups(ARGS)
