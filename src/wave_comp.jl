

module wave_comp

using ITensors

include("component_def.jl")



##########
#
# Sets up the initial conditions for one qubit.
# Can take an array for the density matrix of the qubit 
# or the character 1,0,+,- to identify the state of the qubit. 
#
##########

function Init_stW(i::Index,a)
  A = ITensor(ComplexF64,i)
  if isa(a,Array)
  A[i(1)]=parse(ComplexF64,a[1])
  if size(a,1)>=2
    A[i(2)]=parse(ComplexF64,a[2])
  else
    A[i(2)]= 0.0 + 0.0im
  end
  else 
   if a == '1'
    A[i(1)]= 0.0 + 0.0im
    A[i(2)]= 1.0 + 0.0im
   end
   if a == '0'
    A[i(1)]= 1.0 + 0.0im
    A[i(2)]= 0.0 + 0.0im
   end
   if a == '+'
    A[i(1)]= 1/sqrt(2) + 0.0im
    A[i(2)]= 1/sqrt(2) + 0.0im
   end
   if a == '-'
    A[i(1)]= 1/sqrt(2) + 0.0im
    A[i(2)]= -1/sqrt(2) + 0.0im
   end
  end 

  return A
  
end

##########
#
# Following are function to set the tensors for the gates. 
#
##########


function UGateW(i::Index, j::Index,n1,n2,n3)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=exp(1.0im*(n1+n2)/2.0)*cos(n3/2.0)
  A[i(1),j(2)]=-exp(-1.0im*(n1-n2)/2.0)*sin(n3/2.0)
  A[i(2),j(1)]=exp(1.0im*(n1-n2)/2.0)*sin(n3/2.0)
  A[i(2),j(2)]=exp(1.0im*(n1+n2)/2.0)*cos(n3/2.0)
  
  return A

end 

function CUGateW(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(1),k(2),l(2)]=1.0

  A[i(2),j(2),k(1),l(1)]=exp(1.0im*(n1+n2)/2.0)*cos(n3/2.0)
  A[i(2),j(2),k(1),l(2)]=-exp(-1.0im*(n1-n2)/2.0)*sin(n3/2.0)
  A[i(2),j(2),k(2),l(1)]=exp(1.0im*(n1-n2)/2.0)*sin(n3/2.0)
  A[i(2),j(2),k(2),l(2)]=exp(1.0im*(n1+n2)/2.0)*cos(n3/2.0)
  
  return A

end 


function CCXGateW(i::Index, j::Index, k::Index, l::Index, m::Index, n::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1),m(1),n(1)]=1.0
  A[i(1),j(1),k(1),l(1),m(2),n(2)]=1.0

  A[i(1),j(1),k(2),l(2),m(1),n(1)]=1.0
  A[i(1),j(1),k(2),l(2),m(2),n(2)]=1.0
  
  A[i(2),j(2),k(1),l(1),m(1),n(1)]=1.0
  A[i(2),j(2),k(1),l(1),m(2),n(2)]=1.0
  
  A[i(2),j(2),k(2),l(2),m(2),n(1)]=1.0
  A[i(2),j(2),k(2),l(2),m(1),n(2)]=1.0
  
  return A

end 

function HGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0/sqrt(2.0)
  A[i(1),j(2)]=1.0/sqrt(2.0)
  A[i(2),j(1)]=1.0/sqrt(2.0)
  A[i(2),j(2)]=-1.0/sqrt(2.0)
  
  return A

end 

function XGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(2)]=1.0
  A[i(2),j(1)]=1.0

  return A

end 

function YGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(2)]=1.0im
  A[i(2),j(1)]=-1.0im

  return A

end 

function ZGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=-1.0

  return A

end 

function IGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=1.0

  return A

end 

function SGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=1.0im

  return A

end 

function TGateW(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = 1.0
  A[i(2), j(2)] = cos(pi/4.0) -1.0 * sin(pi/4.0)im

  return A

end

function CNotGateW(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(1),k(2),l(2)]=1.0
  A[i(2),j(2),k(1),l(2)]=1.0
  A[i(2),j(2),k(2),l(1)]=1.0

  return A

end 

function RxGateW(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = cos(Phase / 2.0)^2.0
  A[i(1), j(2)] = (sin(Phase) / 2.0)im
  A[i(2), j(1)] = -1.0 * (sin(Phase) / 2.0)im
  A[i(2), j(2)] = (sin(Phase / 2.0)^2.0)

  return A

end 

function RyGateW(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = (cos(Phase / 2.0)^2.0)
  A[i(1), j(2)] = sin(Phase) / 2.0
  A[i(2), j(1)] = sin(Phase) / 2.0
  A[i(2), j(2)] = (sin(Phase / 2.0)^2.0)
  return A

end 

function RzGateW(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = 1.0
  A[i(2), j(2)] = cos(Phase) +1.0 * (sin(Phase))im

  return A

end 

function CRGateW(i::Index, j::Index, k::Index, l::Index,C)
  A = ITensor(ComplexF64,i,j,k,l)
  
  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(1),k(2),l(2)]=1.0
  A[i(2),j(2),k(1),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=exp(2.0 * pi *1.0im / (2.0^(C+1.0)))

  return A

end 

function CPhaseGateW(i::Index, j::Index, k::Index, l::Index,Phase)
  A = ITensor(ComplexF64,i,j,k,l)
  
  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(1),k(2),l(2)]=1.0
  A[i(2),j(2),k(1),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=(cos(Phase) +1.0 * sin(Phase)im)


  return A

end 

function CZGateW(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(1),k(2),l(2)]=1.0
  A[i(2),j(2),k(1),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=-1.0

  return A

end 

function SWGateW(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(2),l(1)]=1.0
  A[i(2),j(1),k(1),l(2)]=1.0
  A[i(2),j(2),k(2),l(2)]=1.0

  return A

end 






end
