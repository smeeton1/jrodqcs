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


struct registar
   name
   l
   qorc
   init
end


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
   return A 
end



end
