
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

 
   if a[i][1] == "Y"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+YGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,YGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,YGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "Z"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+ZGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,ZGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,ZGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "Rx"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+RxGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,RxGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,RxGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "Ry"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+RyGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,RyGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,RyGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "Rz"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])&&(order(Ham[end])==2)
         Ham[end]=Ham[end]+RzGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2])
       else
         push!(Ham,RzGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
       end
     else
       push!(Ham,RzGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     end
   end
   
   if a[i][1] == "XM"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1])&&(order(Ham[end])==1)
         Ham[end]=Ham[end]+XMeasure(d[parse(Int64,a[i][2]),1])
       else
         push!(Ham,XMeasure(d[parse(Int64,a[i][2]),1] )
       end
     else
       push!(Ham,XMeasure(d[parse(Int64,a[i][2]),1]) )
     end
   end
   
   if a[i][1] == "YM"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1])&&(order(Ham[end])==1)
         Ham[end]=Ham[end]+YMeasure(d[parse(Int64,a[i][2]),1])
       else
         push!(Ham,YMeasure(d[parse(Int64,a[i][2]),1] )
       end
     else
       push!(Ham,YMeasure(d[parse(Int64,a[i][2]),1]) )
     end
   end
   
   if a[i][1] == "ZM"
     if !isempty(Ham)
       if hasinds(Ham[end],d[parse(Int64,a[i][2]),1])&&(order(Ham[end])==1)
         Ham[end]=Ham[end]+ZMeasure(d[parse(Int64,a[i][2]),1])
       else
         push!(Ham,ZMeasure(d[parse(Int64,a[i][2]),1] )
       end
     else
       push!(Ham,ZMeasure(d[parse(Int64,a[i][2]),1]) )
     end
   end
   
 end
 return Ham
end

function Ten_Add(Tens)
 fin=[]
 N = size(Tens,1)
 println("Add N=",N)
 re=[]
 #append!(fin,Tens[1])
 append!(re,0)
 k=0
 for i=1:N
  # println(Tens[i])
   if !(i in re)
   push!(fin,Tens[i])

   k=k+1
   index=inds(fin[k])
   ad = true
   j=i+1
   while ad
    println("j=",j)
    if j>N
      ad = false
      
    elseif hasinds(Tens[j],index)
     if order(fin[k])==order(Tens[j])
        fin[k]=fin[k]+Tens[j]
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

function Ten_split(T)
 fin=[]
 N = size(T,1)
 println("Split N=",N)
 for i=1:N
 println(order(T[i]))
  if order(T[i])>2
   index=inds(T[i])
   U,S,V =svd(T[i],(index[1],index[2])) 
   U=U*S
   push!(fin,U)
   println(order(V))
   while order(V)>3
    index=inds(V)
    n=size(index,1)
    U,S,V =svd(V,(index[1],index[2],index(n)))
    U=U*S
    push!(fin,U)
   end
   push!(fin,V)
  
  else
   push!(fin,T[i])
  end
 
 end

 return fin
 
end

end
