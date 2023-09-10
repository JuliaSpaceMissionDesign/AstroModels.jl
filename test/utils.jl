using DelimitedFiles

function read_gdf_file(filename::AbstractString)
    file = open(filename, "r")
    current_line = 1 
    end_of_head_line = 0

    while !eof(file)
        current_line += 1
        line = split(readline(file))

        # empty lines 
        length(line) < 1 && continue    
        if line[1] == "end_of_head"
            end_of_head_line = current_line
            break
        end
    end 

    return readdlm(filename; skipstart = end_of_head_line-1)
end