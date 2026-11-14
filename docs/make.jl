# Build the LeafPhotosynthesis.jl manual (Documenter + DocumenterVitepress).
#
#   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
#   julia --project=docs docs/make.jl
#   python3 -m http.server --bind 127.0.0.1 --directory docs/build/1 8000
#
# The build runs every example in Float32 and Float64, checks the results, writes the chart
# data into docs/src/public/generated/ (git-ignored), runs the doctests and builds the site.
# Set LEAF_DOCS_FORMAT=html for a plain Documenter HTML site (no Node.js needed), written to
# docs/build/ and readable from the file system.
using Documenter
using DocumenterVitepress
using LeafPhotosynthesis

include(joinpath(@__DIR__, "gallery_data.jl"))

const DOCS = @__DIR__
const REPO = "github.com/RSFraLab/LeafPhotosynthesis.jl"
const GENERATED = joinpath(DOCS, "src", "public", "generated")
const IN_CI = get(ENV, "CI", "false") == "true"
const HTML_FALLBACK = get(ENV, "LEAF_DOCS_FORMAT", "vitepress") == "html"

rm(GENERATED; force = true, recursive = true)
mkpath(GENERATED)
gallery = GalleryData.write_gallery_data(GENERATED)
GalleryData.check_gallery_data(gallery, GENERATED,
    joinpath(DOCS, "src", "public", "leaf"), joinpath(DOCS, "src"))

DocMeta.setdocmeta!(LeafPhotosynthesis, :DocTestSetup, :(using LeafPhotosynthesis); recursive = true)

const PAGES = Any[
    "Home" => "index.md",
    "Example gallery" => "examples.md",
    "Flexas (2002) explorer" => "flexas_2002_explorer.md",
    "Design" => [
        "Design decisions" => "design.md",
        "References and provenance" => "references.md",
    ],
    "Walkthroughs" => [
        "Bonan 11.1 exercises" => "bonan_exercises.md",
        "Core gas exchange" => "core_exchange.md",
        "Equilibrium and sensitivities" => "equilibrium_and_sensitivities.md",
        "Biochemical turnover" => "biochemical_turnover.md",
        "Oxygen response" => "oxygen_response.md",
        "Photosystem excitation" => "excitation.md",
        "PSII fluorescence" => "fluorescence.md",
        "Cytochrome b6f stress" => "cytochrome_fluorescence.md",
        "Flexas (2002) water-stress SIF" => "flexas_2002_water_stress.md",
    ],
    "API reference" => [
        "Overview" => "api/index.md",
        "Leaf composition and environments" => "api/leaf_composition_and_environments.md",
        "Temperature and compensation" => "api/temperature_and_compensation.md",
        "Excitation and electron transport" => "api/excitation_and_electron_transport.md",
        "Carbon assimilation and turnover" => "api/carbon_assimilation_and_turnover.md",
        "Fluorescence and PAM" => "api/fluorescence_and_pam.md",
        "Stomata, diffusion and equilibrium" => "api/stomata_diffusion_and_equilibrium.md",
    ],
]

format = if HTML_FALLBACK
    Documenter.HTML(;
        mathengine = Documenter.KaTeX(),
        prettyurls = IN_CI,
        edit_link = nothing,
        assets = ["public/leaf/loader.js"],
        size_threshold = nothing,
        size_threshold_warn = nothing,
    )
else
    MarkdownVitepress(;
        repo = REPO,
        devbranch = "main",
        devurl = "dev",
        # Locally (not CI) build one site with base "/" in docs/build/1.
        deploy_decision = IN_CI ? nothing : Documenter.DeployDecision(; all_ok = false),
        description = "Composable leaf photosynthesis, gas exchange and fluorescence in Julia.",
    )
end

makedocs(;
    sitename = "LeafPhotosynthesis.jl",
    modules = [LeafPhotosynthesis],
    remotes = Dict(dirname(DOCS) => (Documenter.Remotes.GitHub("RSFraLab", "LeafPhotosynthesis.jl"), "main")),
    format,
    pages = PAGES,
    clean = true,
    doctest = true,
    checkdocs = :public,
    warnonly = false,
)

# Deployment happens only in CI (GitHub Actions); a local build stops here.
if IN_CI
    if HTML_FALLBACK
        deploydocs(; repo = REPO * ".git", devbranch = "main", push_preview = true)
    else
        DocumenterVitepress.deploydocs(;
            repo = REPO * ".git",
            target = joinpath(DOCS, "build"),
            devbranch = "main",
            devurl = "dev",
            push_preview = true,
        )
    end
end
