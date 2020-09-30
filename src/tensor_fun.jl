module tensor_fun

using ITensors
include("component_def.jl")


export Ten_Add, Ten_split, line_mps, Contract_Lines, Contract_Node, Q_Meas 
 
 
 
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

function Search_edge(T,n)
  a=[]
  for i=1:length(T)
   if n in T[i]
    push!(a,i)
   end
  end
  return a
end

function Search_edge_1(T,n)
  a=[]
  for i=1:length(T)
   if n == T[i][1]
    push!(a,i)
   end
  end
  return a
end

function Edge_contract(B,H,E,n)
 i=E[2]-n
 B=B*H[i]
end

function Contract_Node(Q,H,E)
#need to contract by column/ depth
#need to skip rows when no gate in (check multi row gates to see if skip)
 N=length(Q)+length(H)
 n=length(Q)
 Et=copy(E)
 Bh=copy(Q)
 println(Et)
 D=0
 N=[]
 for i=1:n
  push!(N,i)
 end
 D=1
 while length(Et)>0 && D < 5
  println("D=",D)
  println("N=",N)
  println("E_L=",length(Et))
  println(Et)
  D=D+1
  a=[]
  n1=[]
  for i=1:n
   push!(a,Search_edge_1(Et,N[i]))
  end
  println(a)
  
  for i=1:length(a)
   for j=1:length(a[i])
     push!(n1,[i,j,Et[a[i][j]][2]])
   end
  end
  println(n1)
  for i=1:length(n1)
    b= Search_edge(Et,n1[i][3])
    n2=0
    if length(b)>2
     for j=2:length(b)
      println(Et[b[j]])
      for k=1:length(n1)
       println(k!=i && !(n1[k][3] in Et[b[j]]) )
       if k!=i && !(n1[k][3] in Et[b[j]])
        n2=n2+1
       end
      end
     end
     println(n2)
     if n2>3
      println("hello",n1[i][1],n1[i][2])
      deleteat!(a[n1[i][1]],[n1[i][2]])
      println("bye")
     end
    end
  end
  println(a)
  for i=1:length(a)
   if length(a[i])>0
      Edge_contract(Bh[i],H,Et[a[i][1]],n)
      N[i] = Et[a[i][1]][2]
   end
  end
  for i=1:length(a)
     if length(a[i])>1
      for k=2:length(a[i])
      for j=1:length(N)
        if i!=j && N[j] == Et[a[i][k]][2]
          Bh[i]=Bh[i]*Bh[j]
          Bh[j]=Bh[i]
        end
      end
      end
     end
  end
  l=0
  for i=1:length(a)
   if length(a[i])>0
    for j=1:length(a[i])
     deleteat!(Et,(a[i][j]-l))
     l=l+1
    end
   end
  end

  end

  
 for i=1:n
  index=inds(Bh[i])
  Bh[i]=Par_Trac(Bh[i], index[i])
    #println("o= ",order(A[j]))
 end
 return Bh
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
