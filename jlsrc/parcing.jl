
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

       push!(Ham,HGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "X"

       push!(Ham,XGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "CN"

       push!(Ham,CNotGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2],d[parse(Int64,a[i][3]),1],d[parse(Int64,a[i][3]),2]) )
     
    end

 
   if a[i][1] == "Y"

       push!(Ham,YGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Z"

       push!(Ham,ZGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Rx"

       push!(Ham,RxGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Ry"

       push!(Ham,RyGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "Rz"

       push!(Ham,RzGate(d[parse(Int64,a[i][2]),1],d[parse(Int64,a[i][2]),2]) )
     
   end
   
   if a[i][1] == "XM"

       push!(Ham,XMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
   if a[i][1] == "YM"
     
       push!(Ham,YMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
   if a[i][1] == "ZM"
     
       push!(Ham,ZMeasure(d[parse(Int64,a[i][2]),1]) )
     
   end
   
 end
 return Ham
end

function Ten_Add(Tens)
 fin=[]
 N = size(Tens,1)
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
 for i=1:N
  if order(T[i])>2
   index=inds(T[i])
   U,S,V =svd(T[i],(index[1],index[2])) 
   U=U*S
   push!(fin,U)
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
