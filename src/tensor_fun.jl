module tensor_fun

using ITensors
include("component_def.jl")


export Ten_Add, Ten_split, line_mps, Q_Meas
 
 
 
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


function line_mps(Q,T)
 N=size(T,1)
 index=inds(Q)
 for i=1:N
   if hasinds(T[i],index)
     Q=Q*T[i]
  
   end
  
 end

 return Q

end


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


function Contract_Lines(Q,T)
 N=size(Q,1)
 A=[]
 for i=1:N
   push!(A,line_mps(Q[i],T))
 end
 B=A[1]*A[2]
 for i=3:N
   B=B*A[i]
 end
 index=inds(B)
 #println(index)
 #println(length(index))
 for j=1:N
  A[j]=Par_Trac(B, index[j])
  #println("o= ",order(A[j]))
 end
 return A

end

function Q_Meas(Q)
 i=0
 while i<30
  x=bitrand()
  if x[1]
   if rand() < abs(Q[1])
    return 1
   end
  else
   if rand() < abs(Q[4])
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

end
