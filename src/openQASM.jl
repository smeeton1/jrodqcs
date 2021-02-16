module openQASM

using ITensors
using Random
using DelimitedFiles

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")
include("wave.jl")
include("wave_comp.jl")
include("noteb.jl")


##################################################################
#                                                                #
# Structures to hold the data for gate function and registrars   #
#                                                                #
##################################################################

mutable struct registar
   name
   l
   qorc
   init
end

mutable struct Name_gate
    name 
    paramaters
    inputs
    gate
    Name_gate(v) = (name = v, paramaters = [], inputs = [], gate = [])
end


##################################################################
#                                                                #
# Basic function for sorting through input files                 #
#                                                                #
##################################################################

function findchar(a,s)
   i=1
    while (a[i]!=s) && i<length(a)
        i=i+1
    end
    if a[i] == s
        return i-1
    else
        return i
    end
end

function findChArray(a,s)
    go=true
    j=1
    while go && j<length(a)
           if findchar(a[j],s)<length(a[j]) 
                go=false
           else
                j=j+1
           end
    end
    return j
end


function splitequal(a)
    n=length(a)
    m=findchar(a[2],'=')
    s="0"
    if m<length(a[2])
        if m<length(a[2])+2
            s=a[2][m+2:end-1]
        else
            s=a[3][1:end-1]
        end
    end
    if n>2
        for i=3:n
           m=findchar(a[i],'=')
           if m<length(a[i])+2
               s=a[i][m+2:end-1] 
           end
        end
    end
    return s
    
end

function get_contained_string(a,n,target=',')
    n1=findchar(a[n+2:end],target)+n+1
    s=a[n+2:n1]
    
   return s,n1 
end

function split_onchar(a,s)
             n=findchar(a,s) 
             s=a[1:n]
             ch=a[n+2:end]
    return s, ch
end

function compress_string(a)
   b=""
   for i=1:length(a) 
        b=b*a[i]
   end
   return b
end


##################################################################
#                                                                #
# Functions to set up gate and registrars data                   #
#                                                                #
##################################################################

function qregis(a)
    n = findchar(a[2],'[')
    if n < length(a[2]) 
        n1 = findchar(a[2],']')
        m=parse(Int64,a[2][n+2:n1])
        A = registar(a[2][1:n],m,true,'0'^m)
    else
        n = findchar(a[1],'[')
        n1 = findchar(a[1],']')
        m=parse(Int64,a[1][n+2:n1])
        A = registar(a[2][1:end-1],m,true,'0'^m)
    end
    s=splitequal(a)
    if s != "0"
        A.init=s
    end
   return A 
end

function cregis(a)
    n = findchar(a[2],'[')
    if n < length(a[2]) 
        n1 = findchar(a[2],']')
        m=parse(Int64,a[2][n+2:n1])
        A = registar(a[2][1:n],m,false,'0'^m)
    else
        n = findchar(a[1],'[')
        n1 = findchar(a[1],']')
        m=parse(Int64,a[1][n+2:n1])
        A = registar(a[2][1:end-1],m,false,'0'^m)
    end
    s=splitequal(a)
    if s != "0"
        A.init=s
    end
   return A 
end

function add_gate(a)
    m=findchar(a[1][2],'(')
    A=Name_gate(a[1][2][1:m])
    k=2
    if m<length(a[1][2])
        l=findchar(a[1][2],':')+1
        j=true
        if l<length(a[1][2])
           push!(A.paramaters,a[1][2][l+1:end])
           k=k+1 
        else
           k=k+1 
        end
        while j==true
            if  findchar(a[1][k],')')<length(a[1][k])
                j=false
                push!(A.paramaters,a[1][k][1:end-1])
                k=k+1
            else
                push!(A.paramaters,a[1][k])
                k=k+1
            end
        end
    else
        k=3
    end
    
    if a[1][end]=="{"
        n=length(a[1])-1
    else
        n=length(a[1])
    end
    
    for i=k:n
       push!(A.inputs,a[1][i]) 
    end
    
    i=2
    j=true
    while j==true
        if a[i][end]=="}"
            j=false
            if length(a[i])>1
             push!(A.gate,a[i][1:end-1])
            end
        else
            push!(A.gate,a[i])
            i=i+1
        end
    end
    
    return A
end

end
