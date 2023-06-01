"""
    PDSGravityHarmonicsData{N, T} <: AbstractGravityModelData{T}

Data container for spherical harmonics data from the Geosciences Node of NASA's 
[Planetary Data System (PDS)](https://pds-geosciences.wustl.edu/default.html).

### Fields 

- `μ` -- gravitational constant (km³/s²)
- `radius` -- reference radius of the spherical harmonic development (km)
- `max_degree` -- maximum degree of the model 
- `normalized` -- normalized coefficients
- `uncertainty` -- if the coefficients has uncertainty
- `Clm, Slm` -- spherical harmonics constants
- `Clm_unc, Slm_unc` -- spherical harmonics constants uncertainty
"""
struct PDSGravityHarmonicsData{N, T} <: AbstractGravityModelData
    μ::T # km³/s²
    radius::T # km
    max_degree::Int
    normalized::Bool
    uncertainty::Bool

    Clm::Vector{T}
    Slm::Vector{T}
    Clm_unc::Vector{T}
    Slm_unc::Vector{T}
end

function Base.show(io::IO, gd::PDSGravityHarmonicsData{N, T}) where {N, T}
    println(io, "SphericalHarmonicsData{$N, $T}(norm=$(gd.normalized), uncertainty=$(gd.uncertainty))")
end

"""
    parse_data(::Type{T}, ::Type{PDSGravityHarmonicsData}, filename::AbstractString;
        maxdegree::Int=-1, uncertainty::Bool=false) where T

Parse data associated to a `GravityHarmonics` model from NASA's PDS SHA files.

### Fields

- `filename` -- name of the file to be parsed 
- `maxdegree` -- maximum degree of the harmonics to be saved in memory after reading. Optional, default is -1, i.e. all.
- `uncertainty` -- flag to save spherical harmonics constants

### References

- https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/ggmro_095a_sha.lbl
- https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200a_sha.lbl
- https://pds-geosciences.wustl.edu/messenger/mess-h-rss_mla-5-sdp-v1/messrs_1001/data/shadr/jgmess_160a_sha.tab
"""
function parse_data(::Type{T}, ::Type{PDSGravityHarmonicsData}, filename::AbstractString; 
    maxdegree::Int=-1, uncertainty::Bool=false) where T 

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
       
    # Number of rows 
    rows = Int(1/2*(maxdegree+1)*(maxdegree+2)-1)
    
    # Load delimited file
    raw_sha = readdlm(filename, ','; skipstart=1)
    raw_sha = raw_sha[1:rows, :]

    # Extract coefficients
    @views Clm = convert(Vector{T}, raw_sha[:, 3])
    @views Slm = convert(Vector{T}, raw_sha[:, 4])
    if uncertainty
        @views ClmUnc = convert(Vector{T}, raw_sha[:, 5])
        @views SlmUnc = convert(Vector{T}, raw_sha[:, 6])
    else 
        ClmUnc = [] 
        SlmUnc = []
    end
    return PDSGravityHarmonicsData{maxdegree, T}(
        μ, radius, maxdegree, normalized, uncertainty, 
        Clm, Slm, ClmUnc, SlmUnc
        )
end
