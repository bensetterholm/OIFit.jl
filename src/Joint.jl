module Joint

export +, ==, length, Vis2
export annulus, disk, gauss, pseudolorentzian

import Base: +, ==, length, merge

using CMPFit
using Measurements
using OIFITS
using Random
using Statistics

include("data.jl")
include("model.jl")
include("fit.jl")

end # module
