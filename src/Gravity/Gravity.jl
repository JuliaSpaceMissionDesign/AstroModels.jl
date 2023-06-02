using JSMDInterfaces.Models: AbstractJSMDModelData, AbstractJSMDModel
import JSMDInterfaces.Models: parse_data, parse_model
using JSMDInterfaces.Errors
using Logging
using DelimitedFiles

include("abstract.jl")
include("parse_harmonics.jl")
include("harmonics.jl")