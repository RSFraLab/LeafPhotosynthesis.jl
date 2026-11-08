# Host-side output; kept outside the biophysics examples.
function write_tables(output_directory, curves)
    mkpath(output_directory)
    for (name, rows) in curves
        open(joinpath(output_directory, name * ".csv"), "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for row in rows
                println(io, join((value === nothing ? "" : value for value in values(row)), ','))
            end
        end
    end
    return output_directory
end
