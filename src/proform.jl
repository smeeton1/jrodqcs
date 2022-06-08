#
#
# Functions for measuring and testing performance.
#
#




module proform

using ITensors
include("noteb.jl")
include("parcing.jl")



###########################################################################
#
# Function to estimate the number of flops used.
#
###########################################################################


function check_ind(ind1,ind2)
    fsum=0.0

    for i=1:length(ind1)
     for j=1:length(ind2)
        if id(ind1[i])==id(ind2[j])
            fsum+=2*dim(ind1[i])
        end
     end
    end
    
    return fsum

end


function flop_estimate(G1,G2)

    fsum=0.0

    ind1=inds(G1)
    ind2=inds(G2)

    fsun = check_ind(ind1,ind2)


    return fsum

end

################################################################
#
# calculates the fidelity of two tensors using eigenvaules
#
################################################################
function fidelity(A,B)
   
    
    AM=parcing.makematrix(A)
    BM=parcing.makematrix(B)
    la,va =eigen(AM, sortby = x -> -abs(x))
    lb,vb =eigen(BM, sortby = x -> -abs(x))
    
    println(la)
    println(lb)
    
    sum=0
    for i=1:length(lb)
        sum=sum+sqrt(sqrt(lb[i])*la[i]*sqrt(lb[i]))
    end
    
    return sum^2
end



###########################################################################
#
# Function to generate a random circuit with N qubits and of depth of M.
#
###########################################################################

# Creates N initial states set to the state 0
function init_set0(T)
    
    s="0"
    
    for i=2:T.qubit_N
        s*="0"
    end
    
    noteb.Set_init(T,s);
    
    return T
    
end

# Add a 2 qubit gate to each qubit
# if shift is true then the first qubit is skipped 
function add_2BG(T,shift)
    
    N=Int(floor(T.qubit_N/2))
    
    for i=1:N
        
        n=rand()
        if shift
            j = 2*i
        else
            j = 1+2*(i-1)
        end
        
        if j+1<=T.qubit_N
            if 0.0<=n<0.5
                noteb.Add_gate(T,"CN",[j,j+1]);
            else
                noteb.Add_gate(T,"CZ",[j,j+1]);
            end
        end
        
    end
    
    return T
    
end

# Adds a 1 qubit gate to each qubit
function add_1BG(T,shift)
    
    for i=1:T.qubit_N 
        n=rand()
        
        if 0.0<=n<0.2
            noteb.Add_gate(T,"H",i);
        end
        
        if 0.2<=n<0.4
            noteb.Add_gate(T,"Z",i);
        end
        
        if 0.4<=n<0.6
            noteb.Add_gate(T,"X",i);
        end
        
        if 0.6<=n<0.8
            noteb.Add_gate(T,"Y",i);
        end
        
        if 0.8<=n<=1.0
            noteb.Add_gate(T,"I",i);
        end
          
    end
    
    return T
    
end

# Creates a circuit of the form 
#
# Q-G-G-G- -G-G-G....
#     |       |
# Q-G-G-G-G-G-G-G....
#         |
# Q-G- -G-G-G- -G....
#
#Creates a circuit that has all the qubits connected.
function Random_Circuit(T,Depth)
    
    shift=false
    
    T = init_set0(T)
    
    for i=1:Depth
       
        if mod(i,2)==0
            T = add_2BG(T,shift)
            if shift
                shift = false
            else
                shift = true
            end
        else
            T = add_1BG(T,shift)
        end
        
        
    end
    
    return T
    
end

end
