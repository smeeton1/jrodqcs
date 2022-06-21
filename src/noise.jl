#
#
# This module adds noise to the gates. The noise is added a gates that follows the original gate. Three types of noise are included.
#
#

module noise

using ITensors
using ITensors
using Random
using DelimitedFiles
include("component_def.jl")
include("parcing.jl")

# Sets a tensor for use in noise creation.
function set_Aeth(i::Index, j::Index,e)

    A = ITensor(ComplexF64,i,j)
    
    A[i(1),j(1)] = 1.0
    A[i(2),j(2)] = sqrt(1-e)
    A[i(3),j(3)] = sqrt(1-e)
    A[i(4),j(4)] = 1-e
    
    return A
    

end

function set_Aeth2(i::Index, j::Index,e)

    A = ITensor(ComplexF64,i,j)
    A[i(1),j(1)] = 1.0
    A[i(2),j(2)] = 1.0
    A[i(3),j(3)] = 1.0
    A[i(4),j(4)] = 1.0
    A[i(1),j(4)] = e
    
    return A
    

end

# Creates a dephasing noise tensor.
function DPh_noise_1qb(e,G)

    indexs=inds(G)
    g= (1-e)*component_def.IGate(indexs[1],indexs[2]) + e*component_def.ZGate(indexs[1],indexs[2])
    return g #Dag_mult(g,G,4)

end


# Creates a depolarizing noise tensor.
function DPo_noise_1qb1(e,G)

    indexs=inds(G)
    g=set_Aeth(indexs[1],indexs[2],e)
    return g #Mat_mult(G,g,4)

end

function DPo_noise_1qb2(e,G)

    indexs=inds(G)
    g =set_Aeth2(indexs[1],indexs[2],e)
    return g#Mat_mult(G,g,4)

end


# Adds amplitude damping noise to the gate.
function AD_noise_1qb(e,G)

    indexs=inds(G)
    G=(1-e)*G+e/2*component_def.IGate(indexs[1],indexs[2])
    return G
    
end

# Creates a dephasing noise tensor for a 2 qubit gate.
function DPh_noise_2qb(e,G)

    indexs=inds(G)
    g = (1-e)*component_def.IGate(indexs[1],indexs[2]) + e*component_def.ZGate(indexs[1],indexs[2])
    g1 = (1-e)*component_def.IGate(indexs[3],indexs[4]) + e*component_def.ZGate(indexs[3],indexs[4])
    G1 = g*g1
    return G1 #Dag_mult(g,G,4)

end

# Adds amplitude damping noise to a 2 qubit gate.
function AD_noise_2qb(e,G)

    indexs=inds(G)
    G=(1-e)*G+e/2*component_def.IGate(indexs[1],indexs[2])*component_def.IGate(indexs[3],indexs[4])
    return G
    
end

# Creates a depolarizing noise tensor for a 2 qubit gate.
function DPo_noise_2qb1(e,G)

    indexs=inds(G)
    g=set_Aeth(indexs[1],indexs[2],e)*set_Aeth(indexs[3],indexs[4],e)
    return g #Mat_mult(G,g,4)

end

function DPo_noise_2qb2(e,G)

    indexs=inds(G)
    g =set_Aeth2(indexs[1],indexs[2],e)*set_Aeth2(indexs[3],indexs[4],e)
    return g #Mat_mult(G,g,4)

end

#################################################################################
#
# The following functions choose between one or two qubit gates for the noise.
#
##################################################################################


function noise_DPO1(G,e)
    
    if order(G) == 2
        g = DPo_noise_1qb1(e,G)
    else
        g = DPo_noise_2qb1(e,G)
    end
 
    return g
    
end

function noise_DPO2(G,e)
    
    if order(G) == 2
        g = DPo_noise_1qb2(e,G)
    else
        g = DPo_noise_2qb2(e,G)
    end
 
    return g
    
end

function noise_AD(G,e)
    
    if order(G) == 2
        G = AD_noise_1qb(e,G)
    else
        G = AD_noise_2qb(e,G)
    end
    return G
    
end

function nosie_DPH(G,e)
    
    if order(G) == 2
        g = DPh_noise_1qb(e,G)
    else
        g = DPh_noise_2qb(e,G)
    end
    
    return g
    
end



# Sets the edge for the noise tensors.
function fix_edge(edge,n,m,order)
    
    add=true
    for i=1:length(edge)
        if edge[i][1]==n+3
            hold = edge[i][2]
            edge[i][2] = m+3
            push!(edge,[m+3,hold])
            add=false
        end
    end
    
    if add
        if order == 2
            push!(edge,[n+3,m+3])
        else
            push!(edge,[n+3,m+3])
            push!(edge,[n+3,m+3])
        end
    end

    return edge

end



#function to add type of noise to the gates.
function add_nosie_gate(T,e,n)
    
    if length(e)==3
        if e[1]!=0
          push!(T.gates,nosie_DPH(T.gates[n],e[1]))
          fix_edge(T.edge,n,length(T.gates),order(T.gates[n]))
        end
        if e[2]!=0
          push!(T.gates,noise_AD(T.gates[n],e[2])) 
        end
        if e[3]!=0
          push!(T.gates,noise_DPO1(T.gates[n],e[3]))
          fix_edge(T.edge,n,length(T.gates),order(T.gates[n]))
          push!(T.gates,noise_DPO2(T.gates[n],e[3]))
          fix_edge(T.edge,n,length(T.gates),order(T.gates[n]))
        end
        
    elseif length(e)==2
        if e[1]!=0
          push!(T.gates,nosie_DPH(T.gates[n],e[1]))
          fix_edge(T.edge,n,length(T.gates),order(T.gates[n]))
        end
        if e[2]!=0
          push!(T.gates,noise_AD(T.gates[n],e[2]))  
        end
        
    else 
      push!(T.gates,nosie_DPH(T.gates[n],e))
      fix_edge(T.edge,n,length(T.gates),order(T.gates[n]))
    end
   
    return T
    
end

# Adds noise for a constant e.
function add_nosie_eco(T,e)
   
   n=length(T.gates)
   for i=1:n
        T=add_nosie_gate(T,e,i)
   end
    
   return T 
    
end

# Adds noise for a variable e.
function add_nosie_ech(T,e)

   for i=1:length(T.gates) 
        T=add_nosie_gate(T,e[i],i)
   end
    
   return T 
    
end


#####################################################
#
# Master Noise adding function
#
######################################################
function add_nosie(T,e)
    
    if length(e)==length(T.gates)
        T=add_nosie_ech(T,e)
    else
        T=add_nosie_eco(T,e)
    end

    return T
    
end


end
