module measure

using ITensors

include("component_def.jl")
include("wave_comp.jl")
include("parcing.jl")
include("tensor_fun.jl")


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

function Get_Measure_stat(T,qn,form="Wave",t)

  A=tensor_fun.Par_Trac(T, qn)
  Q=[]
  if form == "Wave"
    push!(Q,A[1]/sqrt(A[1]*A[1]+A[2]*A[2]))
    push!(Q,A[2]/sqrt(A[1]*A[1]+A[2]*A[2]))
  else
    push!(Q,A[1]/sqrt(A[1]*A[1]+A[4]*A[4]))
    push!(Q,A[4]/sqrt(A[1]*A[1]+A[4]*A[4]))
  end
  
  n=Q_Meas(Q,t)
  
  if n==1 && form=="Wave"
    A[1] = 1
    A[2] = 0
  elseif n==0 && form=="Wave"
    A[1] = 0
    A[2] = 1
  else if n==1 && form!="Wave"
    A[1] = 1
    A[2] = 0
    A[3] = 0
    A[4] = 0
  else if n==0 && form!="Wave"  
    A[1] = 0
    A[2] = 0
    A[3] = 0
    A[4] = 1
  end
  
  return A
end

function create_newstate(T,qn,M)

  A = T
  temp = component_def.Trace(qn)
  A = A*temp
  A = A*M

  return A

end

function measure_out(T,qn,form="Wave",t)

  M = Get_Measure_stat(T,qn,form="Wave",t)
  A = create_newstate(T,qn,M)
  if form == "Wave"
    if M[1] == 1
        c=0
    else
        c=1
    end
  else
      if M[1] == 1
        c=1
    else
        c=0
    end
  
  end

  return A, c
end
