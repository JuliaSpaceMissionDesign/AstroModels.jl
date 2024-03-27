export GravityHarmonicsICGEMData, GravityHarmonicsPDSData, parse_data

# ----------------------------------
# PDS gravity harmonics data (SHA)
# ----------------------------------

"""
    GravityHarmonicsPDSData{N, T} <: AbstractGravityHarmonicsData{N, T}

Data container for spherical harmonics data from the Geosciences Node of NASA's 
[Planetary Data System (PDS)](https://pds-geosciences.wustl.edu/default.htm).

### Fields 

- `μ` -- gravitational constant (km³/s²)
- `radius` -- reference radius of the spherical harmonic development (km)
- `max_degree` -- maximum degree of the model 
- `normalized` -- normalized coefficients
- `Clm, Slm` -- spherical harmonics constants

### References 
- https://pds-geosciences.wustl.edu/dataserv/gravity_model_desc.htm
"""
struct GravityHarmonicsPDSData{N, T} <: AbstractGravityHarmonicsData{N, T}
    μ::T # km³/s²
    radius::T # km
    max_degree::Int

    normalized::Bool

    # Coefficients
    Clm::Matrix{T}
    Slm::Matrix{T}
end

function Base.show(io::IO, gd::GravityHarmonicsPDSData{N, T}) where {N, T}
    println(io, "GravityHarmonicsPDSData{$N, $T}(norm=$(gd.normalized))")
end

"""
    parse_data(::Type{T}, ::Type{GravityHarmonicsPDSData}, filename::AbstractString;
        maxdegree::Int=-1) where T

Parse data associated to a `GravityHarmonics` model from NASA's PDS SHA files.

### Parameters

- `filename` -- name of the file to be parsed 
- `maxdegree` -- maximum degree of the harmonics to be saved in memory after reading. 
                 Optional, default is -1 indicating that all the coefficients shall be loaded.

### References
- [https://pds-geosciences.wustl.edu/dataserv/gravity_model_desc.htm](https://pds-geosciences.wustl.edu/dataserv/gravity_model_desc.htm)
- [https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/ggmro_095a_sha.lbl](https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/ggmro_095a_sha.lbl)
- [https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.lbl](https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.lbl)
- [https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.tab](https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.tab)
"""
function parse_data(::Type{T}, ::Type{GravityHarmonicsPDSData}, filename::AbstractString; 
    maxdegree::Int=-1) where T 

    fn = open(filename, "r")

    line = strip.(split(strip(readline(fn)), ","))
    radius = parse(T, line[1])
    μ = parse(T, line[2])
    deg = parse(Int, line[4])
    normalized = parse(Int, line[6]) == 1

    ref_lon = parse(T, line[7])
    ref_lon != zero(T) && @warn "Reference longitude for the a spherical harmonics model is $ref_lon (deg)"
    ref_lat = parse(T, line[8])
    ref_lat != zero(T) && @warn "Reference latitude for the a spherical harmonics model is $ref_lat (deg)"
    
    # Maximum degree check
    if maxdegree == -1 
        maxdegree = deg  
    end 

    # Load delimited file
    raw_sha = readdlm(filename, ','; skipstart=1)
    raw_sha = raw_sha[1:end, :]

    # Extract coefficients
    Clm = zeros(maxdegree+1, maxdegree+1)
    Slm = zeros(maxdegree+1, maxdegree+1)
    for row in eachrow(raw_sha) 
        i = Int(row[1]+1)
        j = Int(row[2]+1)
        if i <= maxdegree+1
            if j <= maxdegree+1
                Clm[i, j] = row[3]
                Slm[i, j] = row[4]
            end
        else 
            break
        end
    end

    return GravityHarmonicsPDSData{maxdegree, T}(μ, radius, maxdegree, normalized, Clm, Slm)
end

"""
    GravityHarmonicsICGEMData{N, T} <: AbstractGravityHarmonicsData{N, T}

Data container for ICGEM spherical harmonics gfc data.

### Fields 
- `μ` -- gravitational constant (km³/s²)
- `radius` -- reference radius of the spherical harmonic development (km)
- `max_degree` -- maximum degree of the model 
- `normalized` -- normalization 
- `tide_system` -- tides used in the model development
- `Clm, Slm` -- spherical harmonics constants

### References 
- [http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf](http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf)
"""
struct GravityHarmonicsICGEMData{N, T} <: AbstractGravityHarmonicsData{N, T}
    μ::T # km³/s²
    radius::T # km
    max_degree::Int

    normalized::Bool 
    tide_system::Symbol 

    # Coefficients
    Clm::Matrix{T}
    Slm::Matrix{T}
end

function Base.show(io::IO, gd::GravityHarmonicsICGEMData{N, T}) where {N, T}
    println(io, "GravityHarmonicsICGEMData{$N, $T}(norm=$(gd.normalized), tide=$(gd.tide_system))")
end

"""
    parse_data(::Type{T}, ::Type{GravityHarmonicsICGEMData}, filename::AbstractString;
        maxdegree::Int=-1) where T

Parse data associated to a `GravityHarmonics` model from International Centre for Global 
Earth Models `gfc` models.

### Parameters

- `filename` -- name of the file to be parsed 
- `maxdegree` -- maximum degree of the harmonics to be saved in memory after reading. 
                 Optional, default is -1 indicating that all the coefficients shall be loaded.

### References
- [http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf](http://icgem.gfz-potsdam.de/ICGEM-Format-2023.pdf)
"""
function parse_data(::Type{T}, ::Type{GravityHarmonicsICGEMData}, filename::AbstractString; 
    maxdegree::Int=-1) where T 
    file = open(filename, "r")

    # --------------------------------------------------------------------------------------
    # Parse header 
    header_line_start = 1 
    header_line_end = 0 
    current_line = 0 
    boh_found = false # begin of head found, optional
    eoh_found = false # end of head found, required
    header = String[]

    # Look for the header keywords begin and end of head
    while !eof(file)
        current_line += 1
        line = strip.(split(readline(file)))

        # empty lines 
        length(line) < 1 && continue

        if line[1] == "begin_of_head"
            if boh_found 
                throw(
                    LoadError(
                        filename, current_line, 
                        "Invalid ICGEM format: two `begin_of_head` keywords found!"
                    )
                )
            end 
            header_line_start = current_line
            boh_found = true 
        elseif line[1] == "end_of_head"
            eoh_found = true
            header_line_end   = current_line
            break
        else boh_found && (current_line != header_line_start)
            push!(header, join(line, " "))
        end
    end

    # `end_of_head` keyword is mandatory.
    if !eoh_found
        throw(
            LoadError(
                filename, current_line,
                "Invalid ICGEM format: the mandatory keyword `end_of_head` was not found!"
            )
        )
    end

    DATA =  Dict(
        :product_type => :unknown, :modelname => :unknown, :gravity_constant => zero(T),
        :radius => zero(T), :max_degree => 0, :tide_system => :zero_tide,
        :norm => :fully_normalized,
    )

    MANDATORY_KEYS = (
        :product_type, :modelname, :gravity_constant, :radius, :max_degree
    )

    # Mandatory
    jheader = join(header, "\n")
    for k in MANDATORY_KEYS
        m = match(Regex("$(k) .*[\\d\\w]"), jheader)

        if m === nothing 
            throw(ErrorException("Invalid ICGEM format: no keyword $k found!"))
        end

        matched = split(m.match)
        if k in (:gravity_constant, :radius, :max_degree)
            DATA[k] = parse(Float64, replace(matched[2], r"[D,d]" => "e"))
        else 
            DATA[k] = Symbol(split(matched[2], ".")[1])
        end
    end

    if DATA[:product_type] != :gravity_field
        throw(ErrorException("Only ICGEM files of type :gravity_field are supported."))
    end

    # Optional 
    for k in (:tide_system, :norm)
        m = match(Regex("$(k) .*[\\d\\w]"), jheader)

        if m !== nothing 
            matched = split(m.match)
            if matched[1] == "tide_system"
                DATA[k] = Symbol(matched[2])
                if !(DATA[k] in (:zero_tide, :tide_free, :unknown))
                    @warn "Invalid ICGEM format: the tide system is not valid! Assuming :unknown."
                    DATA[k] = :unknown
                end
            elseif matched[1] == "norm"
                DATA[k] = Symbol(matched[2])
                if !(DATA[k] in (:fully_normalized, :unnormalized))
                    @warn "Invalid ICGEM format: the norm is not valid! Assuming :fully_normalized."
                    DATA[k] = :fully_normalized
                end
            end
        end
    end

    # --------------------------------------------------------------------------------------
    # Read coefficients
    raw = readdlm(filename; skipstart=header_line_end+1)

    deg = Int(DATA[:max_degree])  
    # Maximum degree check
    if maxdegree == -1 
        maxdegree = deg
    elseif maxdegree > deg 
        throw(
            ArgumentError(
                "Maximum degree provided ($maxdegree) is greater than the one in the file ($deg)."
            )
        )
    end 

    # Check if there is any number to be converted. In this case, we assume that
    # the number is written in FORTRAN format.
    raw_coefs = @view(raw[:, 4:6])
    for i in eachindex(raw_coefs)
        if !(typeof(raw_coefs[i]) <: AbstractFloat)
            try
                data             = replace(raw_coefs[i], r"[D,d]" => "e")
                raw_coefs[i] = parse(T, data)
            catch e
                throw(
                    ErrorException(
                        "Invalid ICGME format: could not convert the coefficient \"$(raw_coefs[i])\" to a $(T)."
                    )
                )
            end
        end
    end

    # Parse the coefficients matrices
    # At the moment only STATIC coefficients (in time) are considered
    Clm = zeros(maxdegree+1, maxdegree+1)
    Slm = zeros(maxdegree+1, maxdegree+1)
    for row in eachrow(raw)
        # looping over gfc or gfct lines only
        if row[1] in ("gfc", "gfct")
            i = row[2]+1
            j = row[3]+1
            if i <= maxdegree+1
                if j <= maxdegree+1
                    Clm[i, j] = row[4]
                    Slm[i, j] = row[5]
                end
            end
        end
    end

    return GravityHarmonicsICGEMData{maxdegree, T}(
        DATA[:gravity_constant]/1e9, # to km^3/s^2
        DATA[:radius]/1e3, # to km 
        maxdegree, 
        DATA[:norm] == :fully_normalized, 
        DATA[:tide_system], 
        Clm, Slm
    )
end
