#
#
# Contains Definition of structures and functions need for notebook functionality. 
#
#

module noteb

using ITensors
using Random
using DelimitedFiles

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")
include("wave.jl")
include("wave_comp.jl")

export qc_network, set_init, add_gate, add_measure, wave_in, wave_out, density_in, density_out, qmeasure_out, contract
export Split_N, Split_T, Split


##########
#
# Definition of qc_network structure used to store the circuit information for use in a notebook.
# The structure is initiated using qc_network(q) with q being the number of qubits.
#
##########

mutable struct qc_network

 qubit_N    
 edge 
 gates      
 init_state 
 cbits
 out
 measure     
 indexs 
 record
 verbose      
 solver   
 form   
 
 qc_network(q,f="Wave",v=false) = (qubit_N = q, indexs = parcing.Index_setup(q,f), init_state = [], cbits = [], gates = [], measure = [],out = [], edge = [], record =[[string(q)]],verbose = v,solver  = "Node",form    = f)
end


# Sets the initial condition for the network. Takes network and a string containing initial state.  
function Set_init(T,s)
 push!(T.init_state,parcing.set_qinit(s,T.indexs,T.qubit_N,T.form))
 if length(T.record)<2
   push!(T.record,[s])
 else
   s=[s," "]
   splice!(T.record,2:1,[s])
 end


end
 
function Set_cinit(T,s)
 for i=1:length(s)
    push!(T.cbits,parse(Int,s[i]))
 end
 if length(T.record)<3
   push!(T.record,[s])
 else
   s=[s," "]
   splice!(T.record,3:2,[s])
 end


end
 
 
# Finds the node proceeding current node being added. Takes qubit and network. 
function Find_pre_node(n,T)
   n1=n
    for k=1:length(T.edge)
       if n1 == T.edge[k][1] && T.edge[k][2] != T.qubit_N+length(T.gates)
         indexs=inds(T.gates[T.edge[k][2]-T.qubit_N])
         for l=1:length(indexs)
          if (indexs[l]) in inds(T.gates[end]) && T.indexs[n,1] in indexs
           n1 = T.edge[k][2]
          end
         end
       end
       
    end
  return n1

end

function Set_Cbits(T,n,s="0")

    if s == "0"
        for i=1:n
            push!(T.cbits,0)
        end
    else
        for i=1:n
            push!(T.cbits,eval(Meta.parse(s[i])))
        end
    end


end

# Add a set of gates to the network. 
# In the case of one gate the inputs are the network a string for the gate and qubits for the gate.
# If a gate effects more then one qubit then they are given in square brackets [i,j] 
# For multiple gates the are given in square brackets ["1","2',"3"]
# The qubits are given in a similar way. [[i,j],k,l]
function Add_gate(T,g,n) 
 if isa(g,Array) 
  for i=1:length(g)
    if g[i] == "m"
       push!(T.gates,ITensor(ComplexF64,T.indexs[n[i][1],1],T.indexs[n[i][1],2]))
       push!(T.measure,[g[i],n[i]])
    else
        if T.form == "Density" 
            push!(T.gates,parcing.gate_set(g[i],n[i],T.indexs))
            push!(T.measure,-1)
        else
            push!(T.gates,wave.gate_setW(g[i],n[i],T.indexs))
            push!(T.measure,-1)
        end
    end
    push!(T.record,[g[i],string(n[i])])
    if length(T.edge)<1
     if isa(n[i],Array)
      if g[i] == "CR" || g[i] == "CP" || g[i] == "Rx" || g[i] == "Ry" || g[i] == "Rz"
       for j=1:length(n[i])-1
        push!(T.edge,[Int64(n[i][j]),T.qubit_N+1])
       end
      else
       for j=1:length(n[i])
        push!(T.edge,[n[i][j],T.qubit_N+1])
       end
      end
     else
      push!(T.edge,[n[i],T.qubit_N+1])
     end
    else
     if isa(n[i],Array)
       if g[i] == "CR" || g[i] == "CP" || g[i] == "Rx" || g[i] == "Ry" || g[i] == "Rz"
        for j=1:length(n[i]) -1
         push!(T.edge,[Find_pre_node(Int64(n[i][j]),T),T.qubit_N+length(T.gates)])
        end
       else
        for j=1:length(n[i])
         push!(T.edge,[Find_pre_node(n[i][j],T),T.qubit_N+length(T.gates)])
        end
       end
     else
      push!(T.edge,[Find_pre_node(n[i],T),T.qubit_N+length(T.gates)])
     end
    end
  end
 else
  if g == "m"
       push!(T.gates,ITensor(ComplexF64,T.indexs[n,1],T.indexs[n,2]))
       push!(T.measure,[g,n])
  else
    if T.form == "Density" 
        push!(T.gates,parcing.gate_set(g,n,T.indexs))
        push!(T.measure,-1)
    else
        push!(T.gates,wave.gate_setW(g,n,T.indexs)) 
        push!(T.measure,-1)
    end
  end
  push!(T.record,[g,string(n)])
  if length(T.edge)<1
    if isa(n,Array)
      if g == "CR" || g == "CP" || g == "Rx" || g == "Ry" || g == "Rz"
       
       for j=1:length(n)-1
        push!(T.edge,[Int64(n[j]),T.qubit_N+1])
       end
      else
       for j=1:length(n)
        push!(T.edge,[n[j],T.qubit_N+1])
       end
      end
     else
      push!(T.edge,[n,T.qubit_N+1])
     end
  else
   if isa(n,Array)
    if g == "CR" || g == "CP" || g == "Rx" || g == "Ry" || g == "Rz"
     for j=1:length(n)-1
       push!(T.edge,[Find_pre_node(Int64(n[j]),T),T.qubit_N+length(T.gates)]);
     end
    else
     for j=1:length(n)
       push!(T.edge,[Find_pre_node(n[j],T),T.qubit_N+length(T.gates)]);
     end
    end
   else
    push!(T.edge,[Find_pre_node(n,T),T.qubit_N+length(T.gates)]);
   end
  end
 end

end
 
function Add_measure(T,g,n) 

 push!(T.measure,(g,n))

end

##########
#
# Functions to write out the initial conditions and the out waves.
#
##########
 
# Writing the initial state wavefunction for the qubits n where n can be an int or array of ints. 
function Wave_in(T,n=0) 
 if n==0
   parcing.write_wave_out(tensor_fun.Tproduct(T.init_state[1]),T.form)
 else
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.init_state[1][n[i]],T.form)
   println(" ")
  end
 else
  parcing.write_wave_out(T.init_state[1][n],T.form)
 end
 end
end

# Writing the initial state density matrix for the qubits n where n can be an int or array of ints.
function Density_in(T,n=0) 
 if n == 0
   parcing.write_density_out(tensor_fun.Tproduct(T.init_state[1]),T.form)
 else
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.init_state[1][n[i]],T.form)
     println(" ")
  end
 else
   parcing.write_density_out(T.init_state[1][n],T.form)
 end
 end
end

# Writing the output state wavefunction for the qubits n where n can be an int or array of ints.
function Wave_out(T,n) 
 if isa(n,Array)
  for i=1:length(n)
   parcing.write_wave_out(T.out[n[i]],T.form)
   println(" ")
  end
 else
  parcing.write_wave_out(T.out[n],T.form)
 end
 
end

# Writing the output state density matrix for the qubits n where n can be an int or array of ints.
function Density_out(T,n)
 if isa(n,Array)
  for i=1:length(n)
     parcing.write_density_out(T.out[n[i]],T.form)
     println(" ")
  end
 else
   parcing.write_density_out(T.out[n],T.form)
 end
 
end
 
# Performs a measurement on qubit n  
function Qmeasure_out(T,n) 

 return tensor_fun.Q_Meas(T.out[1][n])
 
end

##########
#
# Functions to split gates that work on multiple qubits.
#
##########

# Calls the splitting function depending on the value of solver
function Split(T)
 if T.verbose == true
  println(T.solver)
 end
 if T.solver == "Line"
   Split_L(T);
 else
   Split_N(T);
 end
end

function Split_L(T)
 fin=[]
 N = size(T.gates,1)
 hold = ITensor(ComplexF64,T.indexs[1,1])
 for i=1:N
  if order(T.gates[i])>3
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
 T.gates[:] = fin;
end

function Split_N(T)
 
 N = size(T.gates,1)
 for i=1:N
  if order(T.gates[i])>3
   index=inds(T.gates[i])
   l = Int64(length(index)/2)
   U,S,V =svd(T.gates[i],(index[1],index[2])) 
   U=U*S
   T.gates[i]=U
   while order(V)>3
    index=inds(V)
    n=size(index,1)
    U,S,V =svd(V,(index[1],index[2],index(n)))
    U=U*S
    push!(T.gates,U)
    push!(T.measure,T.measure[i])
   end
   push!(T.gates,V)
   push!(T.measure,T.measure[i])
   a=tensor_fun.Search_edge(T.edge,i+T.qubit_N)
   m=length(T.gates)+T.qubit_N
   for j=2:length(a)
    if j<=l 
      T.edge[a[j]][2]=m-l+j
    else 
     for k=2:l
      indexs=inds(T.gates[m-l+k-T.qubit_N])
      if indexs[1] in inds(T.gates[T.edge[a[j]][2]-T.qubit_N])
         T.edge[a[j]][1]=m-l+k
      end
     end
    end
   end
   
   push!(T.edge,[i+T.qubit_N,m-l+2])
   if l>2
    for k=2:l-1
     push!(T.edge,[m-l+k,m-l+k+1])
    end
   end
   
  end
 end

end

##########
#
# Functions to contract the network.
#
##########

function Contract(T,Depth::Int64=1000)
 if T.verbose == true
  println(T.solver)
 end
 if T.solver == "Line"
  push!(T.out,tensor_fun.Contract_Lines(T.init_state[1],T.gates,T.measure,T.cbits,T.verbose))
 end
 if T.solver == "Node"
  push!(T.out,tensor_fun.Contract_Node(T.init_state[1],T.gates,T.edge,T.measure,T.cbits,T.verbose,Depth))
 end

end

##########
#
# Functions to clear a network.
#
##########

function remove_gates(T)
 for i=1:length(T.gates)
   pop!(T.gates);
 end
 for i=1:length(T.edge)
   pop!(T.edge);
 end
end

function remove_init(T)
 for i=1:length(T.init_state)
   pop!(T.init_state);
 end
end

function remove_record(T)
 for i=1:length(T.record)
   pop!(T.record);
 end
end

function empty_network(T)
    remove_gates(T)
    remove_init(T)
    remove_record(T)
end


##########
#
# Functions to load and save circuits to and from file.
#
##########

function Load_circuit(filename)
 a = readdlm(filename)#(parcing.Read_InPutFile(filename)
 if component_def.verbose == true
  println(a)
 end
 T = qc_network(a[1][1]);
 s= a[2]
 Set_init(T,s)
 G=a[3:end,1]
 N=a[3:end,2:end]
 for i=1:length(N)
  if !isa(N[i],Int64)
  n=[]
   for j=1:length(N[i])
    if parcing.check_str(N[i][j])
      push!(n,parse(Int64,N[i][j]))
    end
   
   end
  N[i]=n
  else
  N[i]=N[i]
  end
 end

 Add_gate(T,G,N)
 
 if component_def.verbose == true
  println(T)
 end
 
 return T
end

function Save_circuit(T,filename)

   parcing.Write_OutPutFile(T.record,filename)

end

function graph(n)

end

##########
#
# Functions to change global variables when used in notebook.
#
##########


function change_solver(v)
 component_def.set_solver(v)
 if component_def.verbose == true
  println(component_def.solver)
 end
end

function change_verbose(v)
 component_def.set_verbose(v)
 if component_def.verbose == true
  println(component_def.verbose)
 end
end

function change_tol(v)
 component_def.set_tol(v)
 if component_def.verbose == true
  println(component_def.tol)
 end
end

function change_depth(v)
 component_def.set_depth(v)
 if component_def.verbose == true
  println(component_def.depth)
 end
end

function change_Form(v)
 component_def.set_form(v)
 if component_def.verbose == true
  println(component_def.form)
 end
end

end
