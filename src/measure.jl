module measure

using ITensors
using Random
using DelimitedFiles

include("component_def.jl")
include("wave_comp.jl")
include("parcing.jl")


# Performs a partial trace for index I of tensor T.  
function Par_Trac(T::ITensor, I::Index)
 if hasinds(T,I)

  index=inds(T)

  A=T
  for j=1:length(index)
    if index[j]!=I
     temp= component_def.Trace(index[j])
     A= A*temp
    end

  end
  return A
 else
  #println(I, "is not in T")
  return T
 end


end


# Performs a measurement on qubit Q returning 1 or 0 for state of the qubit
function Q_Meas(Q,t)
 i=0
 while i<t
  x=bitrand()
  if x[1]
   if rand() < abs(Q[1])
    return 1
   end
  else
   if rand() < abs(Q[2])
    return 0 
   end
  end
  i=i+1 
 end
 x =bitrand()
 if x[1]
    return 1
 else
    return 0
 end

end 

function Get_Measure_stat(T,qn,t,v,form)

  A=Par_Trac(T, qn)
  Q=[]
  if form == "Wave"
    push!(Q,A[1]/sqrt(A[1]*A[1]+A[2]*A[2])*A[1]/sqrt(A[1]*A[1]+A[2]*A[2]))
    push!(Q,A[2]/sqrt(A[1]*A[1]+A[2]*A[2])*A[2]/sqrt(A[1]*A[1]+A[2]*A[2]))
  else
    push!(Q,A[1]/(A[1]+A[4]))
    push!(Q,A[4]/(A[1]+A[4]))
  end
  if v == true
   println(Q)
   println(Q[1]+Q[2])
   
  end
  n=Q_Meas(Q,t)
  
  if n==1 && form=="Wave"
    A[1] = 1
    A[2] = 0
  elseif n==0 && form=="Wave"
    A[1] = 0
    A[2] = 1
  elseif n==1 && form!="Wave"
    A[1] = 1
    A[2] = 0
    A[3] = 0
    A[4] = 0
  elseif n==0 && form!="Wave"  
    A[1] = 0
    A[2] = 0
    A[3] = 0
    A[4] = 1
  end
  
  return A
end

function create_newstate(T,qn,M,v)

  A = T
  temp = component_def.Trace(qn) 
  A = A*temp
  sum=0
  if dim(qn)==2
  for i =1:length(A.store)
    sum = sum +A.store[i]^2
  end
  for i =1:length(A.store)
    A.store[i] = A.store[i]/sum
  end 
  else
  N=isqrt(length(A.store))
  O=order(A)
  for i =1:N
    sum = sum +A.store[parcing.Fun_A(i,i,N,O)]
  end
  for i =1:length(A.store)
    A.store[i] = A.store[i]/sum
  end 
  
  end
  A = A*M

  return A

end

function measure_out(T,qn,t=100,form="Wave",v=false)

  M = Get_Measure_stat(T,qn,t,v,form)
  if v == true
    println("M= ",M)
    println("inds T = ",length(inds(T)) )
  end
  if length(inds(T)) > 1
    A = create_newstate(T,qn,M,v)
  else
    A=M
  end
  if v == true
    println("A= ",A)
  end
 
 
    if M[1] == 1
        c=0
    else
        c=1
    end
  if v == true
    println("C= ",c)
  end
 
 

  return A, c
end


end
