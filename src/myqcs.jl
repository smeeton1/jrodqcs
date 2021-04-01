#
#
# Loads the modules used in myqcs.
#
#

module myqcs

using ITensors
using Random
using DelimitedFiles
include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")
include("noteb.jl")
include("wave.jl")
include("wave_comp.jl")
include("openQASM.jl")
include("measure.jl")

end
