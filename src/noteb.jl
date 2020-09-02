module noteb

using ITensors
using Random

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")

export qc_network, set_init, add_gate, add_measure, wave_in, wave_out, density_in, density_out, qmeasure_out, contract

mutable struct qc_network

 qubit_N    
 agency_mat 
 gates      
 init_state 
 out
 measure     
 indexs     
 
 qc_network(q) = (qubit_N = q, indexs = parcing.Index_setup(q), init_state = [], gates = [], measure = [],out = [])
end

function set_init(T,s)

 push!(T.init_state,parcing.set_qinit(s,T.indexs,T.qubit_N))


end
 
 
function add_gate(T,g,n) 
 if isa(g,Array)
  for i=1:length(g)
    push!(T.gates,parcing.gate_set(g[1],n[1],T.indexs))
  end
 else
  push!(T.gates,parcing.gate_set(g,n,T.indexs))
 end
end
 
function add_measure(T,g,n) 

 push!(T.measure,(g,n))

end
 
function wave_in(T,n) 
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.init_state[1][n[i]])
  end
 else
  parcing.write_wave_out(T.init_state[1][n])
 end
 
end
 
function density_in(T,n) 
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.init_state[1][n[i]])
  end
 else
   parcing.write_density_out(T.init_state[1][n])
 end
 
end

function wave_out(T,n) 
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.out[1][n[i]])
  end
 else
  parcing.write_wave_out(T.out[1][n])
 end
 
end
 
function density_out(T,n)
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.out[1][n[i]])
  end
 else
   parcing.write_density_out(T.out[1][n])
 end
 
end
 
function qmeasure_out(T,n) 

 return tensor_fun.Q_Meas(T.out[1][n])
 
end

function split(T)
 fin=[]
 N = size(T.gates,1)
 hold = ITensor(ComplexF64,T.indexs[1,1])
 for i=1:N
  if order(T.gates[i])>2
   index=inds(T.gates[i])
   U,S,V =svd(T.gates[i],(index[1],index[2])) 
   U=U*S
   push!(fin,U)
   push!(T.gates,hold)
   while order(V)>3
    index=inds(V)
    n=size(index,1)
    U,S,V =svd(V,(index[1],index[2],index(n)))
    U=U*S
    push!(fin,U)
    push!(T.gates,hold)
   end
   push!(fin,V)
  
  else
   push!(fin,T.gates[i])
  end
 
 end
 T.gates[:] = fin
end


function contract(T)
 
 push!(T.out,tensor_fun.Contract_Lines(T.init_state[1],T.gates))

end

function graph(n)

end

end
