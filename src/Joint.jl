module Joint

export +, ==, length, Vis2
export annulus, disk, gauss, pseudolorentzian

import Base: +, ==, length, merge

using Measurements
using OIFITS
using Statistics

include("data.jl")
include("model.jl")

end # module
