
module parcing

using ITensors
include("component_def.jl")


export Read_InPutFile, Read_InPutLine, Index_setup, check_str, Ten_Add, Tensor_Setup

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
  push!(Q,Init_st(d[i,2],a[1+i]))
 end
 return Q
end

function Tensor_Setup(N,d,a)
 Ham=[]
 for i=N+2:size(a,1)
   if a[i][1] == "H"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+HGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,HGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,HGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "X"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+XGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,XGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,XGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "CN"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2])&&(order(Ham[end])==4)
         Ham[end]=Ham[end]+CNotGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2])
       else
         push!(Ham,CNotGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2]) )
       end
     else
       push!(Ham,CNotGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2]) )
     end
   end
 end

 return Ham
end

function Ten_Add(Tens)
 fin=[]
 N = size(Tens,1)
 println(N)
 re=[]
 append!(re,0)
 for i=1:N
   println(Tens[i])
   if !(i in re)
   push!(fin,Tens[i])
   ad = true
   j=i+1
   while ad
    println(j)
    if j>N
      ad = false
    elseif hasinds(Tens[j],Tens[i])
     if order(Tens[i])==order(Tens[j])
        fin[i]=fin[i]+Tens[j]
        append!(re,j)
        j=j+1
     else
        ad = false
     end
    else
      j=j+1
    end
   
   end
   end
 
 end

 return fin
 
end

end
