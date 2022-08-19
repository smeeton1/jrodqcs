

include("src/jrodqcs.jl")

using ITensors


##
#
# getting command line arguments and reading input file
#
##

fname   = [];
oname   = [];
verbous = false;
QASM    = true;
stype   = "Wave";


i=1; 
while i<size(ARGS,1)+1
 if ARGS[i] == "-f"
   i=i+1;
   push!(fname,ARGS[i]);
 elseif ARGS[i] == "-h"
   println("This is a pregame to perform a tensor simulation of a  quantum circuit.")
   println("-h displays this help message.")
   println("-f is followed by input file name.")
   println("-o is followed by output file name old files will be over written.")
   println("-v turns on verbose running.")
   println("-D run simulation using density matrices.")
   println("-I file is not in openQASM format.")
   
   exit()
 elseif ARGS[i] == "-o"
   i=i+1;
   push!(oname,ARGS[i]);
 elseif ARGS[i] == "-v"
   verbous = true;
 elseif ARGS[i] == "-D"
   stype = "Density";
 elseif ARGS[i] == "-I"
   QASM    = false;
 end 
end 


if !isempty(fname)
  if QASM 
      QT = myqcs.openQASM.ReadopenQasm(fname[1],verbous);
  else
    a=Read_InPutFile(fname);
  end
else
   println("No input file given.")
end


myqcs.tensor_fun.simple_solver(QT);

if !isempty(oname)
   oname = "out"
end


if QT.form == "Density"
    myqcs.parcing.Write_OutPutFile_dens(QT.out,QT.form,oname)
else
    myqcs.parcing.Write_OutPutFile_wave(QT.out,QT.form,oname)
end




