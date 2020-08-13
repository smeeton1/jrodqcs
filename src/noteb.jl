module noteb

using ITensors

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")

struct qc_network

 qubit_N    :: Int
 agency_mat :: Int[]
 gates      :: ITensor[]
 init_state :: ITensor[]
 measure    :: ITensor[]
 indexs     :: Index[]
 
 qc_network(n,q) = (n = new(), n.qubit_N = q, n.indexs = parcing.Index_setup(q), n.agency_mat = Int[], n.gates = ITensor[], n.init_state = ITensor[], n.measure = ITensor[])
 
 set_init()
 
 add_gate()
 
 add_measure()
 
 contract()
 
 graph()
 
 wave_out()
 
 density_out()
 
 measure_out()
 
 

end

end
