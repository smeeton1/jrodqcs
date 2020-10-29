#
#
# Contains some basic global variable used through out the code and some basic set up functions
#
#

module component_def

using ITensors

export set_verbose, set_solver, set_tol, set_depth

export Init_st, HGate, XGate, CNotGate, XMeasure, ZMeasure, RxGate, RyGate, RzGate, YGate, ZGate, YMeasure, Project1, Project0, Trace

##########
#
# Global Variables
#
##########


verbose = false       # Turns on or off extra output information.
solver  = "Node"      # Chooses the type of solver. (Line, Node)
tol     = 1E-6        # Set tolerance for SVD 
depth   = 1000        # Set maximum depth to be calculated to. If depth is set to higher then circuit depth, 
                      # the whole circuit will be calculated
form    = "Density"   # Choose between working with density matrices or wave function. (Density, Wave)

set_verbose(v) = (global verbose = v) # Allows for the value of verbose to be changed (true, false)
set_solver(v)  = (global solver = v)  # Allows for the value of solver to be changed  (Line, Node)
set_tol(v)     = (global tol = v)     # Allows for the value of tol to be changed     (float)
set_depth(v)   = (global depth = v)   # Allows for the value of depth to be changed   (int)
set_form(v)    = (global form = v)    # Allows for the value of depth to be changed   (Density, Wave)

##########
#
# Sets up the initial conditions for one qubit.
# Can take an array for the density matrix of the qubit 
# or the character 1,0,+,- to identify the state of the qubit. 
#
##########

function Init_st(i::Index,a)
  A = ITensor(ComplexF64,i)
  if isa(a,Array)
  A[i(1)]=parse(ComplexF64,a[1])
  if size(a,1)>=2
    A[i(2)]=parse(ComplexF64,a[2])
        if size(a,1)>=3
            A[i(3)]=parse(ComplexF64,a[3])
                if size(a,1)>=4
                    A[i(4)]=parse(ComplexF64,a[4])
                else
                    A[i(4)]= 0.0 + 0.0im
                end
        else
            A[i(3)]= 0.0 + 0.0im
            A[i(4)]= 0.0 + 0.0im
        end
  else
    A[i(2)]= 0.0 + 0.0im
    A[i(3)]= 0.0 + 0.0im
    A[i(4)]= 0.0 + 0.0im
  end
  else 
   if a == '1'
    A[i(1)]= 1.0 + 0.0im
    A[i(2)]= 0.0 + 0.0im
    A[i(3)]= 0.0 + 0.0im
    A[i(4)]= 0.0 + 0.0im
   end
   if a == '0'
    A[i(1)]= 0.0 + 0.0im
    A[i(2)]= 0.0 + 0.0im
    A[i(3)]= 0.0 + 0.0im
    A[i(4)]= 1.0 + 0.0im
   end
   if a == '+'
    A[i(1)]= 0.5 + 0.0im
    A[i(2)]= 0.5 + 0.0im
    A[i(3)]= 0.5 + 0.0im
    A[i(4)]= 0.5 + 0.0im
   end
   if a == '-'
    A[i(1)]= 0.5 + 0.0im
    A[i(2)]= -0.5 + 0.0im
    A[i(3)]= -0.5 + 0.0im
    A[i(4)]= 0.5 + 0.0im
   end
  end 

  return A
  
end

##########
#
# Following are function to set the tensors for the gates. 
#
##########


function HGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0/2.0
  A[i(1),j(2)]=1.0/2.0
  A[i(1),j(3)]=1.0/2.0
  A[i(1),j(4)]=1.0/2.0
  A[i(2),j(1)]=1.0/2.0
  A[i(2),j(2)]=-1.0/2.0
  A[i(2),j(3)]=1.0/2.0
  A[i(2),j(4)]=-1.0/2.0
  A[i(3),j(1)]=1.0/2.0
  A[i(3),j(2)]=1.0/2.0
  A[i(3),j(3)]=-1.0/2.0
  A[i(3),j(4)]=-1.0/2.0
  A[i(4),j(1)]=1.0/2.0
  A[i(4),j(2)]=-1.0/2.0
  A[i(4),j(3)]=-1.0/2.0
  A[i(4),j(4)]=1.0/2.0
  
  return A

end 

function XGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(4)]=1.0
  A[i(2),j(3)]=1.0
  A[i(3),j(2)]=1.0
  A[i(4),j(1)]=1.0

  return A

end 

function YGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(4)]=1.0
  A[i(2),j(3)]=-1.0
  A[i(3),j(2)]=-1.0
  A[i(4),j(1)]=1.0

  return A

end 

function ZGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=-1.0
  A[i(3),j(3)]=-1.0
  A[i(4),j(4)]=1.0

  return A

end 

function IGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=1.0
  A[i(3),j(3)]=1.0
  A[i(4),j(4)]=1.0

  return A

end 

function SGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1),j(1)]=1.0
  A[i(2),j(2)]=1.0im
  A[i(3),j(3)]=1.0im
  A[i(4),j(4)]=-1.0

  return A

end 

function TGate(i::Index, j::Index)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = 1.0
  A[i(2), j(2)] = cos(pi/4.0) -1.0 * sin(pi/4.0)im
  A[i(3), j(3)] = cos(pi/4.0) + sin(pi/4.0)im
  A[i(4), j(4)] = 1.0

  return A

end

function CNotGate(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(1),l(2)]=1.0
  A[i(1),j(3),k(1),l(3)]=1.0
  A[i(1),j(4),k(1),l(4)]=1.0
  A[i(2),j(1),k(2),l(2)]=1.0
  A[i(2),j(2),k(2),l(1)]=1.0
  A[i(2),j(3),k(2),l(4)]=1.0
  A[i(2),j(4),k(2),l(3)]=1.0
  A[i(3),j(1),k(3),l(3)]=1.0
  A[i(3),j(2),k(3),l(4)]=1.0
  A[i(3),j(3),k(3),l(1)]=1.0
  A[i(3),j(4),k(3),l(2)]=1.0
  A[i(4),j(1),k(4),l(4)]=1.0
  A[i(4),j(2),k(4),l(3)]=1.0
  A[i(4),j(3),k(4),l(2)]=1.0
  A[i(4),j(4),k(4),l(1)]=1.0

  return A

end 

function RxGate(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = cos(Phase / 2.0)^2.0
  A[i(1), j(2)] = (sin(Phase) / 2.0)im
  A[i(1), j(3)] = -1.0 * (sin(Phase) / 2.0)im
  A[i(1), j(4)] = (sin(Phase / 2.0)^2.0)
  A[i(2), j(1)] = (sin(Phase) / 2.0)im
  A[i(2), j(2)] = (cos(Phase / 2.0)^2.0)
  A[i(2), j(3)] = (sin(Phase / 2.0)^2.0)
  A[i(2), j(4)] = -1.0 * (sin(Phase) / 2.0)im
  A[i(3), j(1)] = -1.0 * (sin(Phase) / 2.0)im
  A[i(3), j(2)] = (sin(Phase / 2.0)^ 2.0)
  A[i(3), j(3)] = (cos(Phase / 2.0)^2.0)
  A[i(3), j(4)] = (sin(tPhase) / 2.0)im
  A[i(4), j(1)] = (sin(Phase / 2.0)^2.0)
  A[i(4), j(2)] = -1.0 * (sin(Phase) / 2.0)im
  A[i(4), j(3)] = (sin(Phase) / 2.0)im
  A[i(4), j(4)] = (cos(Phase / 2.0)^2.0)

  return A

end 

function RyGate(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = (cos(Phase / 2.0)^2.0)
  A[i(1), j(2)] = sin(Phase) / 2.0
  A[i(1), j(3)] = sin(Phase) / 2.0
  A[i(1), j(4)] = (sin(Phase / 2.0)^2.0)
  A[i(2), j(1)] = -1.0 * sin(Phase) / 2.0
  A[i(2), j(2)] = (cos(Phase / 2.0)^2)
  A[i(2), j(3)] = -1.0 * (sin(Phase / 2)^2.0)
  A[i(2), j(4)] = sin(Phase) / 2.0
  A[i(3), j(1)] = -1.0 * sin(Phase) / 2.0
  A[i(3), j(2)] = -1.0 * (sin(Phase / 2.0)^2.0)
  A[i(3), j(3)] = (cos(Phase / 2.0)^2.0)
  A[i(3), j(4)] = sin(Phase) / 2.0
  A[i(4), j(1)] = (sin(Phase / 2.0)^2.0)
  A[i(4), j(2)] = -1.0 * sin(Phase) / 2.0
  A[i(4), j(3)] = -1.0 * sin(Phase) / 2.0
  A[i(4), j(4)] = (cos(Phase / 2.0)^2.0)

  return A

end 

function RzGate(i::Index, j::Index, Phase)
  A = ITensor(ComplexF64,i,j)

  A[i(1), j(1)] = 1.0
  A[i(2), j(2)] = cos(Phase) -1.0 * (sin(Phase))im
  A[i(3), j(3)] = cos(Phase) + (sin(Phase))im
  A[i(4), j(4)] = 1.0

  return A

end 

function CRGate(i::Index, j::Index, k::Index, l::Index,C)
  A = ITensor(ComplexF64,i,j,k,l)
  
  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(1),l(2)]=1.0
  A[i(1),j(3),k(1),l(3)]=1.0
  A[i(1),j(4),k(1),l(4)]=1.0
  A[i(2),j(1),k(2),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=exp(2.0 * pi * -1.0im / (2.0^(C+1.0)))
  A[i(2),j(3),k(2),l(3)]=1.0
  A[i(2),j(4),k(2),l(4)]=exp(2.0 * pi * -1.0im / (2.0^(C+1.0)))
  A[i(3),j(1),k(3),l(1)]=1.0
  A[i(3),j(2),k(3),l(2)]=1.0
  A[i(3),j(3),k(3),l(3)]=exp(2.0 * pi * 1.0im / (2.0^(C+1.0)))
  A[i(3),j(4),k(3),l(4)]=exp(2.0 * pi * 1.0im / (2.0^(C+1.0)))
  A[i(4),j(1),k(4),l(1)]=1.0
  A[i(4),j(2),k(4),l(2)]=exp(2.0 * pi * -1.0im / (2.0^(C+1.0)))
  A[i(4),j(3),k(4),l(3)]=exp(2.0 * pi * 1.0im / (2.0^(C+1.0)))
  A[i(4),j(4),k(4),l(4)]=1.0

  return A

end 

function CPhaseGate(i::Index, j::Index, k::Index, l::Index,Phase)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(1),l(2)]=1.0
  A[i(1),j(3),k(1),l(3)]=1.0
  A[i(1),j(4),k(1),l(4)]=1.0
  A[i(2),j(1),k(2),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=(cos(Phase) -1.0 * sin(Phase)im)
  A[i(2),j(3),k(2),l(3)]=1.0
  A[i(2),j(4),k(2),l(4)]=(cos(Phase) -1.0 * sin(Phase)im)
  A[i(3),j(1),k(3),l(1)]=1.0
  A[i(3),j(2),k(3),l(2)]=1.0
  A[i(3),j(3),k(3),l(3)]=(cos(Phase) +1.0 * sin(Phase)im)
  A[i(3),j(4),k(3),l(4)]=(cos(Phase) +1.0 * sin(Phase)im)
  A[i(4),j(1),k(4),l(1)]=1.0
  A[i(4),j(2),k(4),l(2)]=(cos(Phase) -1.0 * sin(Phase)im)
  A[i(4),j(3),k(4),l(3)]=(cos(Phase) +1.0 * sin(Phase)im)
  A[i(4),j(4),k(4),l(4)]=1.0

  return A

end 

function CZGate(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(1),l(2)]=1.0
  A[i(1),j(3),k(1),l(3)]=1.0
  A[i(1),j(4),k(1),l(4)]=1.0
  A[i(2),j(1),k(2),l(1)]=1.0
  A[i(2),j(2),k(2),l(2)]=-1.0
  A[i(2),j(3),k(2),l(3)]=1.0
  A[i(2),j(4),k(2),l(4)]=-1.0
  A[i(3),j(1),k(3),l(1)]=1.0
  A[i(3),j(2),k(3),l(2)]=1.0
  A[i(3),j(3),k(3),l(3)]=-1.0
  A[i(3),j(4),k(3),l(4)]=-1.0
  A[i(4),j(1),k(4),l(1)]=1.0
  A[i(4),j(2),k(4),l(2)]=-1.0
  A[i(4),j(3),k(4),l(3)]=-1.0
  A[i(4),j(4),k(4),l(4)]=1.0

  return A

end 

function SWGate(i::Index, j::Index, k::Index, l::Index)
  A = ITensor(ComplexF64,i,j,k,l)

  A[i(1),j(1),k(1),l(1)]=1.0
  A[i(1),j(2),k(2),l(1)]=1.0
  A[i(2),j(1),k(1),l(2)]=1.0
  A[i(2),j(2),k(2),l(2)]=1.0
  A[i(1),j(3),k(3),l(1)]=1.0
  A[i(1),j(4),k(4),l(1)]=1.0
  A[i(2),j(3),k(3),l(2)]=1.0
  A[i(2),j(4),k(4),l(2)]=1.0
  A[i(3),j(1),k(1),l(3)]=1.0
  A[i(3),j(2),k(2),l(3)]=1.0
  A[i(4),j(1),k(1),l(4)]=1.0
  A[i(4),j(2),k(2),l(4)]=1.0
  A[i(3),j(3),k(3),l(3)]=1.0
  A[i(3),j(4),k(4),l(3)]=1.0
  A[i(4),j(3),k(3),l(4)]=1.0
  A[i(4),j(4),k(4),l(4)]=1.0

  return A

end 

##
#
# Measurements
#
##

function XMeasure(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(2)]=1.0
  A[i(3)]=1.0

  return A

end 

function ZMeasure(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(1)]=1.0
  A[i(4)]=-1.0

  return A

end 

function YMeasure(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(2)]=1.0im
  A[i(3)]=-1.0im

  return A

end 

function Project1(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(4)]=1.0

  return A

end 

function Project0(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(1)]=1.0

  return A

end 

function Trace(i::Index)
  A = ITensor(ComplexF64,i)

  A[i(1)]=1.0
  A[i(4)]=1.0

  return A

end 

end
