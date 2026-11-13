using Documenter
using LeafPhotosynthesis

include(joinpath(@__DIR__, "api_catalog.jl"))

function documenter_inline_math(txt::AbstractString)
    return replace(
        txt,
        r"```[\s\S]*?```|`[^`\n]+`|\$\$[\s\S]+?\$\$|\$(?!\s)([^\$]+?)(?<!\s)\$" =>
            function (m)
                if startswith(m, '$') && !startswith(m, "\$\$")
                    inner = SubString(m, nextind(m, firstindex(m)), prevind(m, lastindex(m)))
                    return "``" * replace(inner, '\n' => ' ') * "``"
                end
                return m
            end,
    )
end

function stage_walkthrough(src_path::String, dst_path::String)
    txt = read(src_path, String)
    # Rewrite repo-relative source/test links and README links for staged Documenter pages.
    txt = replace(
        txt,
        r"\]\(\.\./((?:examples|test|src)/[^)]+)\)" =>
            s"](https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/\1)",
    )
    txt = replace(txt, "](api/README.md)" => "](api/index.md)")
    txt = replace(txt, "](README.md#" => "](index.md#")
    txt = replace(txt, "](README.md)" => "](index.md)")
    write(dst_path, documenter_inline_math(txt))
end

function build_documenter_site()
    mktempdir() do stage_src
        mkpath(joinpath(stage_src, "api"))
        cp(joinpath(@__DIR__, "assets"), joinpath(stage_src, "assets"); force = true)
        cp(joinpath(@__DIR__, "examples.html"), joinpath(stage_src, "examples.html"); force = true)
        for jf in ("oxygen_response_validation.json", "fluorescence_validation.json", "excitation_validation.json")
            src_jf = joinpath(@__DIR__, jf)
            isfile(src_jf) && cp(src_jf, joinpath(stage_src, jf); force = true)
        end

        stage_walkthrough(joinpath(@__DIR__, "README.md"), joinpath(stage_src, "index.md"))
        walkthrough_files = [
            "design.md",
            "bonan_exercises.md",
            "core_exchange.md",
            "biochemical_turnover.md",
            "oxygen_response.md",
            "excitation.md",
            "fluorescence.md",
            "cytochrome_fluorescence.md",
            "flexas_2002_water_stress.md",
            "equilibrium_and_sensitivities.md",
            "references.md",
        ]
        for wf in walkthrough_files
            stage_walkthrough(joinpath(@__DIR__, wf), joinpath(stage_src, wf))
        end

        api_pages_nav = Pair{String,String}[
            "Overview & Index" => "api/index.md",
        ]

        idx_buf = IOBuffer()
        println(idx_buf, "# API Reference")
        println(idx_buf)
        println(idx_buf, "```@docs")
        println(idx_buf, "LeafPhotosynthesis")
        println(idx_buf, "```")
        println(idx_buf)
        println(idx_buf, "## Contents")
        println(idx_buf)
        println(idx_buf, "```@contents")
        println(idx_buf, "Pages = [", join(["\"$(p.slug).md\"" for p in API_PAGES], ", "), "]")
        println(idx_buf, "Depth = 2")
        println(idx_buf, "```")
        println(idx_buf)
        println(idx_buf, "## Symbol Index")
        println(idx_buf)
        println(idx_buf, "```@index")
        println(idx_buf, "Pages = [\"index.md\", ", join(["\"$(p.slug).md\"" for p in API_PAGES], ", "), "]")
        println(idx_buf, "```")
        write(joinpath(stage_src, "api", "index.md"), String(take!(idx_buf)))

        for page in API_PAGES
            rel = "api/$(page.slug).md"
            push!(api_pages_nav, page.title => rel)
            pbuf = IOBuffer()
            println(pbuf, "# $(page.title)")
            println(pbuf)
            summary_doc = documenter_inline_math(strip(page.summary))
            println(pbuf, summary_doc)
            println(pbuf)
            println(pbuf, "```@contents")
            println(pbuf, "Pages = [\"$(page.slug).md\"]")
            println(pbuf, "Depth = 2")
            println(pbuf, "```")
            for (sec_title, syms) in page.sections
                println(pbuf)
                println(pbuf, "## $sec_title")
                println(pbuf)
                println(pbuf, "```@docs")
                for s in syms
                    println(pbuf, s)
                end
                println(pbuf, "```")
            end
            write(joinpath(stage_src, "api", "$(page.slug).md"), String(take!(pbuf)))
        end

        DocMeta.setdocmeta!(
            LeafPhotosynthesis,
            :DocTestSetup,
            :(using LeafPhotosynthesis);
            recursive = true,
        )

        repo_root = dirname(@__DIR__)
        build_dir = joinpath(@__DIR__, "build")
        makedocs(
            sitename = "LeafPhotosynthesis.jl",
            source = stage_src,
            build = build_dir,
            modules = [LeafPhotosynthesis],
            remotes = Dict(
                repo_root => (Documenter.Remotes.GitHub("RSFraLab", "LeafPhotosynthesis.jl"), "main"),
            ),
            clean = true,
            doctest = true,
            checkdocs = :exports,
            format = Documenter.HTML(
                mathengine = Documenter.KaTeX(),
                prettyurls = get(ENV, "CI", nothing) == "true",
                edit_link = nothing,
                collapselevel = 1,
                size_threshold = nothing,
                size_threshold_warn = nothing,
            ),
            pages = Any[
                "Home" => "index.md",
                "Design & Concepts" => [
                    "Design decisions" => "design.md",
                    "References & provenance" => "references.md",
                ],
                "Biophysical Walkthroughs" => [
                    "Bonan 11.1 exercises" => "bonan_exercises.md",
                    "Core gas exchange" => "core_exchange.md",
                    "Equilibrium & sensitivities" => "equilibrium_and_sensitivities.md",
                    "Biochemical turnover" => "biochemical_turnover.md",
                    "Oxygen response" => "oxygen_response.md",
                    "Photosystem excitation" => "excitation.md",
                    "PSII fluorescence" => "fluorescence.md",
                    "Cytochrome b6f stress" => "cytochrome_fluorescence.md",
                    "Flexas (2002) water-stress SIF" => "flexas_2002_water_stress.md",
                ],
                "API Reference" => api_pages_nav,
            ],
        )
    end
end

build_documenter_site()
