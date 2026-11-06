using TOML

"""Read all numerical settings as FT. The reusable physics never reads files."""
function read_bonan_parameters(::Type{FT}, path) where {FT}
    raw = TOML.parsefile(path)
    return Dict(section => Dict(key => value isa Real ? FT(value) : value
        for (key, value) in entries) for (section, entries) in raw)
end
