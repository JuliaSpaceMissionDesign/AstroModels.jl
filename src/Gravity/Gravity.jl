# module Gravity 

using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
import JSMDInterfaces.Models: parse_data, parse_model
using JSMDInterfaces.Errors
using Logging
using DelimitedFiles
using StaticArrays

# Abstract types/methods
include("abstract.jl")

# Gravity harmonics
include("harmonics/utils.jl")
include("harmonics/parse.jl")
include("harmonics/harmonics.jl")
include("harmonics/compute.jl")

# end

long = -120.0 
lat = -30.0
point = SA[ cosd(lat)*cosd(long), cosd(lat)*sind(long), sind(lat) ] * m.radius