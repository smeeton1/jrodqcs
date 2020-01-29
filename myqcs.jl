

include("jlsrc/component_def.jl")
include("jlsrc/parcing.jl")

using ITensors


##
#
# getting command line arguments and reading input file
#
##

fname = []
oname = []
verbous = false

 
for i=1:size(ARGS,1)
 if ARGS[i] == "-f"
   i=i+1
   push!(fname,ARGS[i])
 elseif ARGS[i] == "-h"
   println("This is a pregame to perform a tensor simulation of a  quantum circuit.")
   println("-h displays this help message.")
   println("-f is followed by input file name.")
   println("-o is followed by output file name old files will be over written.")
   println("-v turns on verbose running.")
   
   exit()
 elseif ARGS[i] == "-o"
   i=i+1
   push!(oname,ARGS[i])
 elseif ARGS[i] == "-v"
   verbous = true
 end 
end 


if !isempty(fname)
  a=Read_InPutFile(fname)
end


##
#
# Setting up the indexes for the tensors. Each qubit has two indexes associated with it. All indexes are of length 4.
#
##

if check_str(a[1])
    I=Index_setup(parse(Int64,a[1]))
else
    println("First line is not number of qubits.")
    exit()
end


##
#
# Setting up the tensors for the qubits. Each State tensor has one index.  Each qubit is set to start in |0>
#
##

Q1=ITensor(j)
Q2=ITensor(l)

Q1[j(1)]=1; 
Q2[l(1)]=1;

##
#
# Setting up the gates. These are called from another module I have made and is brought in by include("jlsrc/component_def.jl")
# To create the gate the appropriate indexes need to be supplied.  
#
##

# H and X gate acting on |Q1> with indexes i and j
H1=component_def.HGate(i, j)
X1=component_def.XGate(i, j)

# H  gate acting on |Q2> with indexes k and l
H2=component_def.HGate(k, l)

# Cnot gate acting on both qubits. It has the indexes for both with them order so that the top qubits index come first. 
Cnot=component_def.CNotGate(i, j, k, l)

# Measurements have the index chosen depending on how the indexes balance out through the contractions. 
# Sorry about not being able to give a better answer.
XM1=component_def.XMeasure(i)
XM2=component_def.XMeasure(l)


##
#
# Contraction of the tensors
#
##

# If all indexes are the same for neighboring tensors they are added.

A=(H1+X1)

# When the tensors have different indexes they are contracted.

O=Q1*A*Cnot*X1*H2*XM1*XM2
O=Q2*O

# The second line is there because it was the only way I could get it to stop crashing. I am guessing a bug in itensor.


# Printing answer to screen, which give an x measure of 0. If you remove H2 and change the index of XM2 to k then the output is 1.
println(O)

