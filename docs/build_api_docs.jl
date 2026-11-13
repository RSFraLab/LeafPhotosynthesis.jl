# Generate structured Markdown API reference pages in `docs/api/` directly from
# `LeafPhotosynthesis` docstrings and `docs/api_catalog.jl`.
#
# Usage:
#   julia --project=. docs/build_api_docs.jl          # write docs/api/*.md
#   julia --project=. docs/build_api_docs.jl --check  # verify zero drift and coverage

using LeafPhotosynthesis

include(joinpath(@__DIR__, "api_catalog.jl"))

const REPO_ROOT = dirname(@__DIR__)
const API_DIR = joinpath(@__DIR__, "api")

function symbol_anchor(sym::Symbol)
    return lowercase(replace(String(sym), r"[^A-Za-z0-9_-]" => ""))
end

function symbol_kind(mod::Module, sym::Symbol)
    val = getfield(mod, sym)
    if val isa DataType || val isa UnionAll
        if isabstracttype(val)
            return "Abstract Type"
        elseif val <: Enum
            return "Enum Type"
        else
            return "Struct"
        end
    elseif val isa Function
        return "Function"
    elseif val isa Enum
        return "Enum Value"
    else
        return string(typeof(val))
    end
end

function verify_catalog_coverage(mod::Module)
    expected = sort!(filter(!=(nameof(mod)), collect(names(mod))))
    seen = Symbol[]
    duplicates = Symbol[]
    seen_set = Set{Symbol}()
    for page in API_PAGES
        for (_, syms) in page.sections
            for s in syms
                if s in seen_set
                    push!(duplicates, s)
                else
                    push!(seen_set, s)
                    push!(seen, s)
                end
            end
        end
    end
    !isempty(duplicates) && error("Duplicate symbols in API_PAGES: $(duplicates)")
    missing_syms = setdiff(Set(expected), seen_set)
    extra_syms = setdiff(seen_set, Set(expected))
    !isempty(missing_syms) && error("Exported symbols missing from API_PAGES: $(sort!(collect(missing_syms)))")
    !isempty(extra_syms) && error("Non-exported symbols listed in API_PAGES: $(sort!(collect(extra_syms)))")
    return seen
end

function build_symbol_index()
    idx = Dict{String,Tuple{String,String}}()
    for page in API_PAGES
        fname = page.slug * ".md"
        for (_, syms) in page.sections
            for s in syms
                idx[String(s)] = (fname, symbol_anchor(s))
            end
        end
    end
    return idx
end

function rewrite_prose_line(line::AbstractString, current_file::String, sym_index::Dict{String,Tuple{String,String}})
    # Rewrite Documenter `[`Symbol`](@ref)` cross-references to relative Markdown links.
    out = replace(line, r"\[`([^`]+)`\]\(@ref\)" => function (m)
        inner = match(r"^\[`([^`]+)`\]\(@ref\)$", m).captures[1]
        key = replace(inner, r"^LeafPhotosynthesis\." => "")
        haskey(sym_index, key) || error("Unresolved @ref target `$inner` in docstring")
        target_file, anchor = sym_index[key]
        href = target_file == current_file ? "#$anchor" : "$target_file#$anchor"
        return "[`$inner`]($href)"
    end)
    # Convert Julia Markdown double-backtick inline math ``...`` to GitHub/KaTeX $...$.
    out = replace(out, r"``([^`]+)``" => s"$\1$")
    # Shift top-level docstring subheaders (`# Examples`, `# Arguments`) below `### Symbol`.
    if startswith(out, "# ")
        out = "#### " * SubString(out, 3)
    elseif startswith(out, "## ")
        out = "#### " * SubString(out, 4)
    end
    return out
end

function format_docstring_markdown(raw_text::AbstractString, current_file::String, sym_index::Dict{String,Tuple{String,String}})
    lines = split(chomp(raw_text), '\n'; keepempty = true)
    sig_lines = String[]
    i = 1
    while i <= length(lines)
        l = lines[i]
        if startswith(l, "    ")
            push!(sig_lines, l[5:end])
            i += 1
        elseif isempty(strip(l)) && !isempty(sig_lines) && i < length(lines) && startswith(lines[i + 1], "    ")
            push!(sig_lines, "")
            i += 1
        else
            break
        end
    end
    while i <= length(lines) && isempty(strip(lines[i]))
        i += 1
    end

    buf = IOBuffer()
    if !isempty(sig_lines)
        println(buf, "```julia")
        for sl in sig_lines
            println(buf, sl)
        end
        println(buf, "```")
        if i <= length(lines)
            println(buf)
        end
    end

    in_fence = false
    for j in i:length(lines)
        l = lines[j]
        stripped = lstrip(l)
        if startswith(stripped, "```")
            if !in_fence
                in_fence = true
                if startswith(stripped, "```jldoctest")
                    println(buf, replace(l, "```jldoctest" => "```julia"; count = 1))
                else
                    println(buf, l)
                end
            else
                in_fence = false
                println(buf, l)
            end
        elseif in_fence
            println(buf, l)
        else
            println(buf, rewrite_prose_line(l, current_file, sym_index))
        end
    end
    return String(take!(buf))
end

function one_line_summary(meta::IdDict, mod::Module, sym::Symbol, current_file::String, sym_index::Dict{String,Tuple{String,String}})
    multidoc = meta[Docs.Binding(mod, sym)]
    d = multidoc.docs[first(multidoc.order)]
    raw = join(d.text, "")
    lines = split(chomp(raw), '\n'; keepempty = true)
    i = 1
    while i <= length(lines) && (startswith(lines[i], "    ") || isempty(strip(lines[i])))
        i += 1
    end
    para = String[]
    while i <= length(lines) && !isempty(strip(lines[i]))
        s = strip(lines[i])
        startswith(s, "- ") && break
        push!(para, s)
        i += 1
    end
    joined = join(para, " ")
    if endswith(joined, ":")
        joined = chop(joined) * "."
    end
    # Keep up to first sentence period (ignoring common scholarly abbreviations).
    m = match(r"^(.+?(?<!\bEq)(?<!\bEqs)(?<!\bal)(?<!\bvs)(?<!\bcf)(?<!\bFig)\.)(?:\s|$)", joined)
    sentence = m === nothing ? joined : m.captures[1]
    return rewrite_prose_line(sentence, current_file, sym_index)
end

function render_symbol_entry(meta::IdDict, mod::Module, sym::Symbol, current_file::String, sym_index::Dict{String,Tuple{String,String}})
    b = Docs.Binding(mod, sym)
    haskey(meta, b) || error("Missing Docs.meta entry for $sym")
    multidoc = meta[b]
    buf = IOBuffer()
    println(buf, "### `$sym`")
    println(buf)
    for (k, sig) in enumerate(multidoc.order)
        d = multidoc.docs[sig]
        k > 1 && println(buf)
        raw = join(d.text, "")
        print(buf, format_docstring_markdown(raw, current_file, sym_index))
        src_path = String(d.data[:path])
        src_line = Int(d.data[:linenumber])
        display_src = replace(relpath(src_path, REPO_ROOT), '\\' => '/')
        gh_src = "https://github.com/RSFraLab/LeafPhotosynthesis.jl/blob/main/$display_src#L$src_line"
        println(buf)
        println(buf, "*Kind:* $(symbol_kind(mod, sym)) · *Source:* [`$display_src:$src_line`]($gh_src)")
    end
    return String(take!(buf))
end

function render_process_page(page, meta::IdDict, mod::Module, sym_index::Dict{String,Tuple{String,String}})
    fname = page.slug * ".md"
    buf = IOBuffer()
    println(buf, "# $(page.title)")
    println(buf)
    walk_links = join(["[$(label)]($(href))" for (label, href) in page.walkthroughs], " · ")
    println(buf, "[API Overview](README.md) · Walkthroughs: $walk_links · [Interactive Gallery](../examples.html)")
    println(buf)
    for line in split(strip(page.summary), '\n')
        println(buf, strip(line))
    end
    println(buf)
    println(buf, "## Contents")
    println(buf)
    for (sec_title, syms) in page.sections
        sym_links = join(["[`$s`](#$(symbol_anchor(s)))" for s in syms], ", ")
        println(buf, "- **$sec_title**: $sym_links")
    end
    for (sec_title, syms) in page.sections
        println(buf)
        println(buf, "---")
        println(buf)
        println(buf, "## $sec_title")
        for sym in syms
            println(buf)
            print(buf, render_symbol_entry(meta, mod, sym, fname, sym_index))
        end
    end
    return String(take!(buf))
end

function render_overview_page(meta::IdDict, mod::Module, sym_index::Dict{String,Tuple{String,String}})
    buf = IOBuffer()
    println(buf, "# API Reference")
    println(buf)
    println(buf, "[Repository README](https://github.com/RSFraLab/LeafPhotosynthesis.jl#readme) · [Worked Examples & Plots](../README.md) · [Interactive Gallery](../examples.html) · [Design Decisions](../design.md)")
    println(buf)

    mod_doc = meta[Docs.Binding(mod, nameof(mod))]
    raw_mod = join(mod_doc.docs[first(mod_doc.order)].text, "")
    print(buf, format_docstring_markdown(raw_mod, "README.md", sym_index))
    println(buf)

    println(buf, "## Biophysical Process Map")
    println(buf)
    println(buf, "All $(length(sym_index)) exported types, status codes, and functions in `LeafPhotosynthesis` are organized by biophysical process:")
    println(buf)
    println(buf, "| Biophysical Process | Key Types & Entry Points | Walkthroughs |")
    println(buf, "| --- | --- | --- |")
    for page in API_PAGES
        fname = page.slug * ".md"
        all_syms = reduce(vcat, [syms for (_, syms) in page.sections])
        preview = join(["[`$s`]($fname#$(symbol_anchor(s)))" for s in first(all_syms, min(5, length(all_syms)))], ", ")
        suffix = length(all_syms) > 5 ? " (+$(length(all_syms) - 5) more)" : ""
        walks = join(["[$(label)]($(href))" for (label, href) in page.walkthroughs], ", ")
        println(buf, "| [**$(page.title)**]($fname) | $preview$suffix | $walks |")
    end
    println(buf)
    println(buf, "## Core Conventions")
    println(buf)
    println(buf, "- **Units**: Leaf temperature in `K`; gas partial pressures and Michaelis–Menten constants (`co2`, `oxygen`, `kc`, `ko`, `gamma_star`) in `Pa`; photon, electron, carbon, and oxygen area fluxes in `μmol m⁻² leaf s⁻¹`; molar conductances (`g_sw`, `g_sc`, `g_bw`, `g_bc`) and transpiration in `mol m⁻² leaf s⁻¹`.")
    println(buf, "- **No hidden constants**: Every physical, stoichiometric, and reference parameter is passed explicitly by the caller (see `examples/parameters/*.toml`).")
    println(buf, "- **Type & device contract**: Parameter and environment structs are parametric on scalar type `FT`, allocate zero heap memory per leaf evaluation, propagate `ForwardDiff.Dual` tangents, and broadcast directly as scalars over CPU and CUDA arrays.")
    println(buf)
    println(buf, "## Complete Symbol Directory")
    for page in API_PAGES
        fname = page.slug * ".md"
        println(buf)
        println(buf, "### [$(page.title)]($fname)")
        println(buf)
        println(buf, "| Symbol | Kind | Summary |")
        println(buf, "| --- | --- | --- |")
        for (_, syms) in page.sections
            for s in syms
                kind = symbol_kind(mod, s)
                summ = one_line_summary(meta, mod, s, "README.md", sym_index)
                println(buf, "| [`$s`]($fname#$(symbol_anchor(s))) | $kind | $summ |")
            end
        end
    end
    return String(take!(buf))
end

function build_all(; check_only::Bool = false)
    verify_catalog_coverage(LeafPhotosynthesis)
    sym_index = build_symbol_index()
    meta = Docs.meta(LeafPhotosynthesis)

    outputs = Pair{String,String}[
        joinpath(API_DIR, "README.md") => render_overview_page(meta, LeafPhotosynthesis, sym_index),
    ]
    for page in API_PAGES
        push!(outputs, joinpath(API_DIR, page.slug * ".md") => render_process_page(page, meta, LeafPhotosynthesis, sym_index))
    end

    if check_only
        for (path, expected) in outputs
            isfile(path) || error("Missing generated API file: $path")
            actual = read(path, String)
            actual == expected || error("Stale generated API file: $path (re-run `julia --project=. docs/build_api_docs.jl`)")
        end
        println("Verified $(length(outputs)) API Markdown files covering $(length(sym_index)) exported symbols with zero drift.")
    else
        mkpath(API_DIR)
        for (path, content) in outputs
            write(path, content)
        end
        println("Generated $(length(outputs)) API Markdown files in docs/api/ covering $(length(sym_index)) exported symbols.")
    end
end

build_all(; check_only = ("--check" in ARGS))
