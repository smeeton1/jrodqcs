#
#
# Contains functions for manipulating networks of tensors.
#
#


module tensor_fun

using ITensors
include("component_def.jl")
include("parcing.jl")
include("measure.jl")


export Ten_Add, Ten_split, line_mps, Contract_Lines, Contract_Node, Q_Meas,  Search_edge, Search_edge_1, Edge_contract
 




function Tproduct(T)
  h=T[1]
  for i=2:length(T)
    h=h*T[i]
  end
  return h
end

function measurebit(T,M,H,C,verbose)

  indexs=inds(H)
  
  if dim(indexs[1])==2
    form= "Wave"
  else
    form= "Density"
  end
  if indexs[1] in inds(T)
    qn = indexs[1]
  else
    qn = indexs[2]
  end
  if verbose
   println(form)
  end
  
  T,C[M[3]]=measure.measure_out(T,qn,100,form,verbose)
  return T,C

end

##########
#
# Functions for working with graph edges.
#
##########

# Finds an edge that contains node n.
function Search_edge(T,n)
  a=[]
  for i=1:length(T)
   if n in T[i]
    push!(a,i)
   end
  end
  return a
end

# Finds an edge that starts with node n.
function Search_edge_1(T,n)
  a=[]
  for i=1:length(T)
   if n == T[i][1]
    push!(a,i)
   end
  end
  return a
end

# Contracts a node into the state B according to the edge E, with H being the array of gates and n the number of qubits.
function Edge_contract(B,H,E,n,M,C,verbose)
 i=E[2]-n
 if M[i][1]=="m"
   B,C = measurebit(B,M[i],H[i],C,verbose)
 
 elseif M[i][1] != -1
   if C[M[i][1]] == 1
    B=B*H[i]
   end
 else
    B=B*H[i]
 end
 return B,C
end
 
##########
#
# Functions for manipulating tensor networks.
#
##########
 
# Adds tensor which share indexes. 
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

# Splits tensors that act on more then one qubit.
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

# Performs a line contraction. This function contracts along a single index.
function line_mps(Q,T,M,C,verbose)
 N=size(T,1)
 index=inds(Q)
 for i=1:N
   if hasinds(T[i],index)
     if M[i][1] == "m"
       Q,C=measurebit(Q,M[i],T[i],C,verbose)
     elseif M[i][1] !=-1
       if C[M[i][1]] == 1
         Q=Q*T[i]
       end
     
     else
        Q=Q*T[i]
     end
  
   end
  
 end

 return Q

end

# Performs line contractions for each qubit combines the lines and then traces out the final state for each qubit. 
function Contract_Lines(Q,H,M,Cb,verbose)
 N=size(Q,1)
 A=[]
 C=copy(Q)
 for i=1:N
   push!(A,line_mps(C[i],H,M,Cb,verbose))
 end
 B=A[1]*A[2]
 for i=3:N
   B=B*A[i]
 end
#  index=inds(B)
#  #println(index)
#  #println(length(index))
#  for j=1:N
#   A[j]=measure.Par_Trac(B, index[j])
#   #println("o= ",order(A[j]))
#  end
 return B

end

#Performs a nodal contraction by depth and traces out the final state for each qubit.
function Contract_Node(Q,H,E,M,C,verbose,depth)
#need to contract by column/ depth
#need to skip rows when no gate in (check multi row gates to see if skip)
 N=length(Q)+length(H)
 n=length(Q)
 Et=copy(E)
 Bh=copy(Q)

 N=[]
 for i=1:n
  push!(N,i)
 end
 D=1
 while length(Et)>0 && D < depth
  if verbose == true
   println("Depth=",D)
   println("N=",N)
   println("Edge List Length=",length(Et))
   println("Gate List Length=", length(Bh))
   println("Edge List=",Et)
  end
  D=D+1
  a=[]
  n1=[]
  rm=0
  for i=1:length(N)
   if length(N[i])>1
    for k=1:length(N[i])
     push!(a,Search_edge_1(Et,N[i][k]))
    end
   else
    push!(a,Search_edge_1(Et,N[i]))
   end
  end
  if verbose == true
   println("a=",a)
  end
  
  for i=1:length(a)
   for j=1:length(a[i])
     push!(n1,[i,j,Et[a[i][j]][2]])
   end
  end
  if verbose == true
    println(n1)
  end
  for i=1:length(n1)
   if order(H[n1[i][3]-n])>2
    b= Search_edge(Et,n1[i][3])
    n2=true

     for j=2:length(b)
      for k=1:length(n1)
       #println(k!=i && !(n1[k][3] in Et[b[j]]) )
       if k!=i && (n1[k][3] in Et[b[j]])
        n2=false
       end
      end
     end
     #println(n2)
     if n2
      deleteat!(a[n1[i][1]],[n1[i][2]])
     else
     for j=1:length(b)
      for k=1:length(n1)
       #println(k!=i && !(n1[k][3] in Et[b[j]]) )
       if k!=i && (n1[k][3] in Et[b[j]]) && !(parcing.isin(a,b[j]))
         push!(a,[b[j]])
       end
      end
     end
     end
    end

  end
  #println(a)
  for i=1:length(a)
   if length(a[i])>0
     for j=1:length(N)
     if Et[a[i][1]][1] in N[j] && !(parcing.isin(N,Et[a[i][1]][2]) ) #(isassigned(N,Et[a[i][1]][2] ))
      Bh[j],C=Edge_contract(Bh[j],H,Et[a[i][1]],n,M,C,verbose)
      if length(N[j]) > 1
       for k=1:length(N[j])
        if Et[a[i][1]][1] == N[j][k]
         N[j][k] = Et[a[i][1]][2]
        end
       end
      else
        N[j] = Et[a[i][1]][2]
      end
     else
      for k=1:length(N)
       if j!=k && Et[a[i][1]][1] in N[j] && Et[a[i][1]][2] in N[k]
         Bh[j]=Bh[j]*Bh[k]
         N[j]=[N[j],N[k]]
         deleteat!(Bh,k)
         rm=rm+1
         for l=k:length(N)-1
          N[l]=N[l+1]
         end
         #deleteat!(N,k)
       end
     end
     end
     end
   end
  end
  for i=1:rm
    pop!(N);
  end
  l=0
  sort!(a)
  for i=1:length(a)
   if length(a[i])>0
    for j=1:length(a[i])
     if length(Et)>1
      deleteat!(Et,(a[i][j]-l))
      l=l+1
     else
      pop!(Et);
     end
    end
   end
  end
  for i=1:length(N)
    if length(N[i])>1
     N[i]=parcing.flattenA(N[i])
    end
  end
  #println("N_e=",N)
  end

 A=Bh[1]
 if length(Bh) > 1
  for i=2:length(Bh)
    A=A*Bh[i]
  end
 end
 return A
end

end
