function setup_excitation_test(::Type{FT}) where {FT}
    fluorescence=setup_fluorescence(FT)
    directory=joinpath(@__DIR__, "../examples/parameters")
    bio=read_bonan_parameters(FT, joinpath(directory, "bonan_11_1.toml"))
    core=read_bonan_parameters(FT, joinpath(directory, "core_exchange.toml"))
    base, diffusion=core_leaf(FT, bio, core)
    leaf=LeafModel(;
        photosynthesis = fluorescence.leaf.photosynthesis,
        stomata = base.stomata,
        oxygen = base.oxygen,
    )
    environment=core_environment(FT, bio, core)
    settings=read_bonan_parameters(FT, joinpath(@__DIR__, "parameters/excitation.toml"))
    for key in ("absorbed_par", "psii_fractions", "co2_gamma_factors")
        settings["inputs"][key]=FT.(settings["inputs"][key])
    end
    state=StomatalState(core["experiment"]["initial_conductance"])
    solver=core_solver(FT, core)
    steady=LeafModel(;
        leaf.photosynthesis,
        leaf.oxygen,
        stomata = SteadyStateStomata(leaf.stomata.target),
    )
    s=core["stomata"]
    medlyn=Medlyn{FT}(;
        minimum_conductance = s["minimum_conductance"],
        slope = s["medlyn_slope"],
        diffusivity_ratio = diffusion.diffusivity_ratio,
        micromole_to_mole = s["micromole_to_mole"],
        vpd_floor = s["vpd_floor"],
    )
    bisection=core_solver(FT, core; method = Bisection())
    return (;
        leaf,
        steady,
        diffusion,
        environment,
        state,
        solver,
        fluorescence,
        settings,
        medlyn,
        bisection,
    )
end

# A new excitation family declares its own input contract. This test-only family
# permits already allocated PSII photons outside the environment's PAR interval.
struct ExtraBandTestExcitation{FT} <: AbstractExcitation
    photons::FT
end
LeafPhotosynthesis.valid_excitation(x::ExtraBandTestExcitation, par) =
    isfinite(par)&&par>=zero(par)&&isfinite(x.photons)&&x.photons>=zero(x.photons)
LeafPhotosynthesis.psii_excitation(p, par, x::ExtraBandTestExcitation) =
    valid_excitation(x, par) ? x.photons : oftype(par, NaN)
function LeafPhotosynthesis.electron_transport(
    p::ElectronTransportParameters,
    par,
    jmax,
    x::ExtraBandTestExcitation,
)
    valid_excitation(x, par) || return oftype(par, NaN)
    colimit(p.psii_yield*x.photons, jmax, p.limitation)
end

# Components that support only broadband light implement the excitation argument for
# `BroadbandExcitation` alone; the default argument reaches them, other light fails dispatch.
struct BroadbandOnlyElectrons{FT} <: LeafPhotosynthesis.AbstractLeafScalar
    gain::FT
    psii_fraction::FT
end
LeafPhotosynthesis.electron_transport(p::BroadbandOnlyElectrons, par, jmax, ::BroadbandExcitation) =
    min(p.gain*par, jmax)
struct BroadbandOnlyPhotosynthesis{FT} <: AbstractPhotosynthesis
    gain::FT
end
LeafPhotosynthesis.photosynthesize(
    p::BroadbandOnlyPhotosynthesis,
    env::LeafEnvironment,
    ::BroadbandExcitation,
) = (; net_assimilation = p.gain*env.absorbed_par)
struct BroadbandOnlyOxygen <: AbstractOxygenExchange end
LeafPhotosynthesis.oxygen_exchange(
    ::BroadbandOnlyOxygen,
    p,
    env::LeafEnvironment,
    ::BroadbandExcitation,
) = (; test_value = env.oxygen)

function with_test_electrons(p, electrons)
    C3Photosynthesis(;
        p.temperature,
        p.stoichiometry,
        p.limitation,
        p.assimilation_sign,
        electron_transport = electrons,
    )
end

function excitation_allocations(case, site, input)
    p=case.leaf.photosynthesis
    photosynthesize(p, site, input)
    biochemical_turnover(p, site, input)
    fluorescence_inputs(p, site, input)
    gas_exchange(case.leaf, case.environment, case.state, case.diffusion, input)
    equilibrium_gas_exchange(
        case.steady,
        case.environment,
        case.diffusion,
        case.solver,
        input,
    )
    return @allocated begin
        photosynthesize(p, site, input)
        biochemical_turnover(p, site, input)
        fluorescence_inputs(p, site, input)
        gas_exchange(case.leaf, case.environment, case.state, case.diffusion, input)
        equilibrium_gas_exchange(
            case.steady,
            case.environment,
            case.diffusion,
            case.solver,
            input,
        )
    end
end
