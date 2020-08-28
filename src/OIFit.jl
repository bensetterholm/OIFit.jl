module OIFit

export +, ==, length, Vis2
export annulus, disk, gauss, pseudolorentzian

import Base: +, ==, length, merge

using CMPFit
using Measurements
using OIFITS
using Random
using SpecialFunctions
using Statistics

include("data.jl")
include("model.jl")
include("fit.jl")

end # module
