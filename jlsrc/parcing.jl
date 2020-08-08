
module parcing

using ITensors
include("component_def.jl")


export Read_InPutFile, Read_InPutLine, Index_setup, check_str, Tensor_Setup, testing

function testing()
 println("hello")
end

function check_str(a)
    try
        parse(Int64,a)
        true
    catch
        false
    end
end

function Index_setup(a)

    d=Index[]
    for i=1:a
        append!(d,[Index(4),Index(4)])
    end
    d=reshape(d,(a,2))
    return d

end

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

function Read_InPutLine()



end

function Inital_State(N,d,a)
 Q=ITensor[]
 for i=1:N
  push!(Q,component_def.Init_st(d[i,2],a[1+i]))
 end
 return Q
end

function Tensor_Setup(N,d,a)
 Ham=[]
 for i=N+2:size(a,1)
   if a[i][1] == "H"

       push!(Ham,component_def.HGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "X"

       push!(Ham,component_def.XGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "CN"

       push!(Ham,component_def.CNotGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2]) )
     
    end

 
   if a[i][1] == "Y"

       push!(Ham,component_def.YGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Z"

       push!(Ham,component_def.ZGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Rx"

       push!(Ham,component_def.RxGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Ry"

       push!(Ham,component_def.RyGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Rz"

       push!(Ham,component_def.RzGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "XM"

       push!(Ham,component_def.XMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
   if a[i][1] == "YM"
     
       push!(Ham,component_def.YMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
   if a[i][1] == "ZM"
     
       push!(Ham,component_def.ZMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
 end
 return Ham
end

end
