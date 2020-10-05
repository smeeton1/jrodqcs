#
#
# Contains reading and writing function utility functions
#
#

module parcing

using ITensors
using DelimitedFiles
include("component_def.jl")


export Read_InPutFile, Write_OutPutFile, Index_setup, check_str, Tensor_Setup
export gate_set, set_qinit, measure_set, density_out, write_wave_out, isin, flattenA

##########
#
# Utility functions
#
##########

# Check if a is a string.
function check_str(a)
    try
        parse(Int64,a)
        true
    catch
        false
    end
end

# Checks if n is in A.
function isin(A,n)
 t = false
 for i=1:length(A)
   if isa(A[i],Array)
    for j=1:length(A[i])
     if isin(A[i][j],n)
      t=true
     end
    end
   else
     if A[i] == n
      t=true
     end
   end
 end
 return t

end

# Removes one layer of an array. [[i,j],k]--> [i,j,k]
function flattenA(A)
 b=[]
 for i=1:length(A)
   if length(A[i])>1
      for j=1:length(A[i])
       push!(b,A[i][j])
      end
    else
     push!(b,A[i])
    end
   end

 return b
end


##########
#
# Read in and set up functions
#
##########

# Sets up the indexs for the Tensors
function Index_setup(a)

    d=Index[]
    for i=1:a
        append!(d,[Index(4),Index(4)])
    end
    d=reshape(d,(a,2))
    return d

end

# Reads the input file
function Read_InPutFile(InPutFile)

 a=[]
 open(InPutFile) do f
   i=1
   while !eof(f)
    x=readline(f)
    append!(a,[split(x, " ")])
   end
 end
 
 
 return a

end

# Sets the initial states
function Inital_State(N,d,a)
 Q=ITensor[]
 for i=1:N
  push!(Q,component_def.Init_st(d[i,2],a[1+i]))
 end
 return Q
end

# Gets a tensor for a Gate
function gate_set(s,n,d)


   if s == "H"

       Ham=component_def.HGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "X"

       Ham=component_def.XGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "CN"

       Ham=component_def.CNotGate(d[n[1],1],d[n[1],2],d[n[2],1],d[n[2],2]) 
     
   end

 
   if s == "Y"

       Ham=component_def.YGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "Z"

       Ham=component_def.ZGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "Rx"

       Ham=component_def.RxGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "Ry"

       Ham=component_def.RyGate(d[n,1],d[n,2]) 
     
   end
   
   if s == "Rz"

       Ham=component_def.RzGate(d[n,1],d[n,2]) 
     
   end
   

 return Ham
end

function measure_set(s,n,d)

   
   if s == "X"

       Ham=component_def.XMeasure(d[n,1]) 
     
   end
   
   if s == "Y"
     
       Ham=component_def.YMeasure(d[n,1]) 
     
   end
   
   if s == "Z"
     
       Ham=component_def.ZMeasure(d[n,1]) 
     
   end

 return Ham
end

function set_qinit(S,d,N)
 Q=ITensor[]
 for i=1:N
  push!(Q,component_def.Init_st(d[i,2],S[i]))
 end
 return Q
end


##########
#
# Output functions
#
##########

# Writes an output file
function Write_OutPutFile(a,OutFile)

 io=open(OutFile, "w")
#   for i=1:length(a)
#     println(a[i])
    writedlm(io,a," ")
#  end
 close(io)

end

# Prints a wavefunction to the screen
function write_wave_out(T)
 if sign(T[2]) == 0
  println(sqrt(T[1]),' ', sqrt(T[4]))
 else
  println(sqrt(T[1]),' ', sign(T[2])*sqrt(T[4]))
 end
end

# Prints a density matrix to the screen
function write_density_out(T)

  println(T[1],' ',T[2])
  println(T[3],' ',T[4])

end



end
