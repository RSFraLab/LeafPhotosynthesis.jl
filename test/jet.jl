# Static analysis with JET of the main entry points of both photosynthesis models, in Float32
# and Float64: `test_call` reports possible errors (no method, undefined names), `test_opt`
# reports runtime dispatch. Reports are limited to code in LeafPhotosynthesis.
using JET

function jet_entry_points(::Type{FT}) where {FT}
    bio = read_bonan_parameters(FT, parameter_path)
    core = read_bonan_parameters(FT, joinpath(EXAMPLES, "parameters", "core_exchange.toml"))
    c3_leaf, diffusion = core_leaf(FT, bio, core)
    c3_steady = LeafModel(; photosynthesis = c3_leaf.photosynthesis,
        stomata = SteadyStateStomata(c3_leaf.stomata.target), oxygen = c3_leaf.oxygen)
    gas = core_environment(FT, bio, core)
    site = LeafEnvironment(gas.temperature, gas.absorbed_par, FT(28), gas.oxygen)
    state = StomatalState(FT(0.2))
    solver = core_solver(FT, core)
    vdt = setup_fluorescence(FT)
    psii = fluorescence_inputs(vdt.leaf.photosynthesis, site)

    cyto = setup_cytochrome_stress(FT)
    jb = setup_johnson_berry_2021(FT)
    dual = PrescribedDualExcitation(FT(400), FT(380))
    photosystems = fluorescence_inputs(jb.photo, site, dual)
    jb_result = fluorescence(jb.fl_ku2, photosystems)
    a2, a1 = FT(0.425), FT(0.425)

    c3 = [
        "photosynthesize" => (photosynthesize, (c3_leaf, site)),
        "biochemical_turnover" => (biochemical_turnover, (c3_leaf, site)),
        "oxygen_exchange" => (oxygen_exchange, (c3_leaf, site)),
        "gas_exchange" => (gas_exchange, (c3_leaf, gas, state, diffusion)),
        "leaf_tendencies" => (leaf_tendencies, (c3_leaf, gas, state, diffusion)),
        "equilibrium_gas_exchange" => (equilibrium_gas_exchange, (c3_steady, gas, diffusion, solver)),
        "assimilation_sensitivity" =>
            (assimilation_sensitivity, (c3_leaf.photosynthesis, gas, state, diffusion)),
        "fluorescence_inputs" => (fluorescence_inputs, (vdt.leaf.photosynthesis, site)),
        "fluorescence" => (fluorescence, (vdt.regulated, psii)),
        "gas_exchange, closed stomata" => (gas_exchange, (c3_leaf, gas, StomatalState(zero(FT)), diffusion)),
        "equilibrium_gas_exchange, zero minimum conductance" => (equilibrium_gas_exchange,
            (c3_leaf.photosynthesis, SteadyStateStomata(BallBerry{FT}(0.0, 9.0, 1e-6)), gas, diffusion, solver)),
        "transpiration_sensitivity" => (transpiration_sensitivity,
            (WaterVaporDiffusion(FT(1.4)), state, FT(1400), gas.pressure)),
        "target_conductance" => (target_conductance,
            (c3_leaf.stomata.target, StomatalDrivers(FT(10), FT(30), gas.pressure, FT(0.6), FT(1200)), diffusion)),
        "stomatal_tendency" => (stomatal_tendency,
            (c3_leaf, state, StomatalDrivers(FT(10), FT(30), gas.pressure, FT(0.6), FT(NaN)), diffusion)),
        "oxygen_exchange, inputs" => (oxygen_exchange, (c3_leaf.oxygen, OxygenFluxInputs(FT(16), FT(4), FT(1)))),
    ]
    cytochrome = [
        "photosynthesize" => (photosynthesize, (cyto.cyto_leaf, site)),
        "photosynthesize, dual excitation" => (photosynthesize, (jb.photo, site, dual)),
        "biochemical_turnover" => (biochemical_turnover, (cyto.cyto_leaf, site)),
        "oxygen_exchange" => (oxygen_exchange, (cyto.cyto_leaf, site)),
        "gas_exchange" => (gas_exchange, (cyto.cyto_leaf, gas, state, cyto.diffusion)),
        "equilibrium_gas_exchange" =>
            (equilibrium_gas_exchange, (cyto.cyto_leaf, gas, cyto.diffusion, cyto.solver)),
        "target_conductance, Medlyn" => (target_conductance,
            (cyto.cyto_leaf, StomatalDrivers(FT(10), FT(30), gas.pressure, FT(NaN), FT(1200)), cyto.diffusion)),
        "assimilation_sensitivity" => (assimilation_sensitivity,
            (cyto.cyto_leaf.photosynthesis, gas, state, cyto.diffusion, dual)),
        "fluorescence_inputs" => (fluorescence_inputs, (jb.photo, site, dual)),
        "fluorescence" => (fluorescence, (jb.fl_ku2, photosystems)),
        "state_transition_excitation" =>
            (state_transition_excitation, (jb.rule, jb.fl_ku2, jb.photo, site)),
        "pam_fluorescence" => (pam_fluorescence, (jb.pam_leaked, jb.fl_ku2, jb_result, a2, a1)),
    ]
    return ("C3" => c3, "Cytochrome" => cytochrome)
end

@testset "JET static analysis of the entry points" begin
    for FT in (Float32, Float64), (model, calls) in jet_entry_points(FT)
        @testset "$model $FT" begin
            for (name, (f, args)) in calls
                @testset "$name" begin
                    types = Tuple{map(typeof, args)...}
                    test_call(f, types; target_modules = (LeafPhotosynthesis,))
                    test_opt(f, types; target_modules = (LeafPhotosynthesis,))
                end
            end
        end
    end
end
