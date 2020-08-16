module noteb

using ITensors

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")

export qc_network

struct qc_network

 qubit_N    :: Int
 agency_mat :: Int[]
 gates      :: ITensor[]
 init_state :: ITensor[]
 out
 measure     
 indexs     :: Index[]
 
 qc_network(n,q) = (n = new(), n.qubit_N = q, n.indexs = parcing.Index_setup(q))
 
 set_init(s) = (init_state = parcing.set_qinit(s,indexs))
 
 add_gate(g,n) = (push!(gates,parcing.gate_set(g,n,indexs)))
 
 add_measure(g,n) = (push!(measure,(g,n)))
 
 wave_in(n) = (parcing.wave_out(init_state[n]))
 
 density_in(n) = (parcing.density_out(init_state[n]))

 wave_out(n) = (parcing.wave_out(out[n]))
 
 density_out(n) = (parcing.density_out(out[n]))
 
 qmeasure_out(n) = (tensor_fun.Q_Meas(out[n]))
 
 () -> (set_init;add_gate;add_measure;wave_out;density_out;qmeasure_out)

end

function contract(n::qc_network)
n.gates = tensor_fun.Ten_split(n.gates)
n.out = tensor_fun.Contract_Lines(Q,T)

end

function graph(n::qc_network)

end

end
