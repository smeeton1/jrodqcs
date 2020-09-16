module noteb

using ITensors
using Random

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")

export qc_network, set_init, add_gate, add_measure, wave_in, wave_out, density_in, density_out, qmeasure_out, contract
export Split_N, Split_T, Split

mutable struct qc_network

 qubit_N    
 edge 
 gates      
 init_state 
 out
 measure     
 indexs 
 record
 
 qc_network(q) = (qubit_N = q, indexs = parcing.Index_setup(q), init_state = [], gates = [], measure = [],out = [], edge = [], record =[[string(q)]])
end

function Set_init(T,s)
 push!(T.init_state,parcing.set_qinit(s,T.indexs,T.qubit_N))
 if length(T.record)<2
   push!(T.record,[s])
 else
   s=[s," "]
   splice!(T.record,2:1,[s])
 end


end
 
function Search_edge(T,n)
  a=[]
  for i=1:length(T.edge)
   if n in T.edge[i]
    push!(a,i)
   end
  end
  return a
end
 
function Add_gate(T,g,n) 
 if isa(g,Array)
  for i=1:length(g)
    push!(T.gates,parcing.gate_set(g[i],n[i],T.indexs))
    push!(T.record,[g[i],string(n[i])])
    if length(T.edge)<1
     if isa(n[i],Array)
      for j=1:length(n[i])
       push!(T.edge,[n[i][j],T.qubit_N+1])
      end
     else
      push!(T.edge,[n[i],T.qubit_N+1])
     end
    else
     if isa(n[i],Array)
       for j=1:length(n[i]) 
         n1=n[i][j]
         for k=1:length(T.edge)
          if n1 == T.edge[k][1] && T.edge[i][2] != T.qubit_N+length(T.gates)
           n1 = T.edge[k][2]
          end
         end
         push!(T.edge,[n1,T.qubit_N+length(T.gates)])
       end
     else
      n1=n[i]
      for k=1:length(T.edge)
       if n1 == T.edge[k][1] && T.edge[i][2] != T.qubit_N+length(T.gates)
         n1 = T.edge[k][2]
       end
      end
      push!(T.edge,[n1,T.qubit_N+length(T.gates)])
     end
    end
  end
 else
  push!(T.gates,parcing.gate_set(g,n,T.indexs))
  push!(T.record,[g,string(n)])
  if length(T.edge)<1
    if isa(n,Array)
      for j=1:length(n)
       push!(T.edge,[n[j],T.qubit_N+1])
      end
     else
      push!(T.edge,[n,T.qubit_N+1])
     end
  else
   if isa(n,Array)
    for j=1:length(n)
      n1=n[j]
      for i=1:length(T.edge)
       if n1 == T.edge[i][1] && T.edge[i][2] != T.qubit_N+length(T.gates)
        n1 = T.edge[i][2]
       end
      end
     push!(T.edge,[n1,T.qubit_N+length(T.gates)])
    end
   else
    n1=n
    for i=1:length(T.edge)
     if n1 == T.edge[i][1] && T.edge[i][2] != T.qubit_N+length(T.gates)
       n1 = T.edge[i][2]
     end
    end
    push!(T.edge,[n1,T.qubit_N+length(T.gates)])
   end
  end
 end

end
 
function Add_measure(T,g,n) 

 push!(T.measure,(g,n))

end
 
function Wave_in(T,n) 
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.init_state[1][n[i]])
  end
 else
  parcing.write_wave_out(T.init_state[1][n])
 end
 
end
 
function Density_in(T,n) 
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.init_state[1][n[i]])
  end
 else
   parcing.write_density_out(T.init_state[1][n])
 end
 
end

function Wave_out(T,n) 
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.out[1][n[i]])
  end
 else
  parcing.write_wave_out(T.out[1][n])
 end
 
end
 
function Density_out(T,n)
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.out[1][n[i]])
  end
 else
   parcing.write_density_out(T.out[1][n])
 end
 
end
 
function Qmeasure_out(T,n) 

 return tensor_fun.Q_Meas(T.out[1][n])
 
end

function Split(T)#push V on to the end of gates so then just need to add edges and not move all nodes
 if solver == "simple"
   Split_L(T)
 else
   Split_N(T)
 end
end


function Split_L(T)#push V on to the end of gates so then just need to add edges and not move all nodes
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

function Split_N(T)#push V on to the end of gates so then just need to add edges and not move all nodes
 fin=[]
 N = size(T.gates,1)
 hold = ITensor(ComplexF64,T.indexs[1,1])
 for i=1:N
  if order(T.gates[i])>2
   index=inds(T.gates[i])
   U,S,V =svd(T.gates[i],(index[1],index[2])) 
   U=U*S
   T.gates[i]=U
   while order(V)>3
    index=inds(V)
    n=size(index,1)
    U,S,V =svd(V,(index[1],index[2],index(n)))
    U=U*S
    push!(T.gates,U)
   end
   push!(T.gates,V)
   a=Search_edge(T,i+T.T.qubit_N)
   l=length(a)/2.0
   for j=2:2*l
    if j>=l
      T.edge[a[j]][2]=length(T.gates)-l+j
    else #check indexs for other gate
     for k=2:l
      if inds(T.gates[T.edge[a[j]][2]])==inds(T.gates[length(T.gates)-l+k])
         T.edge[a[j]][1]=length(T.gates)-l+k
      end
     end
    end
   end
   
   push!(T.edge,[i,length(length(T.gates))-l+2])
   if l>2
    for k=2:l-1
     push!(T.edge,[length(length(T.gates))-l+k,length(length(T.gates))-l+k+1])
    end
   end
   
  end
 end

end



function Contract(T)
 if solver == "simple"
  push!(T.out,tensor_fun.Contract_Lines(T.init_state[1],T.gates))
 end
 if solver == "node"
  push!(T.out,tensor_fun.Contract_Node(T.init_state[1],T.gates,T.edge))
 end

end

function Load_circuit(T,filename)
 a = parcing.Read_InPutFile(filename)
 T = qc_network(parse(Int64,a[1][1]))
 Set_init(T,a[2])
 Add_gate(T,a[3:end][1],a[3:end][2:end])
 T.record = a

end

function Save_circuit(T,filename)

   parcing.Write_OutPutFile(T.record,filename)

end

function graph(n)

end

end
