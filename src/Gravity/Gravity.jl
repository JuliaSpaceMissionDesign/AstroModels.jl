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
include("parse_harmonics.jl")
include("harmonics.jl")
include("compute_harmonics.jl")

# end
