#
#
# Contains functions for manipulating networks of tensors.
#
#


module tensor_fun

using ITensors
include("component_def.jl")
include("parcing.jl")


export Ten_Add, Ten_split, line_mps, Contract_Lines, Contract_Node, Q_Meas,  Search_edge, Search_edge_1, Edge_contract
 
# Performs a partial trace for index I of tensor T.  
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

# Performs a measurement on qubit Q returning 1 or 0 for state of the qubit
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
function Edge_contract(B,H,E,n)
 i=E[2]-n
 B=B*H[i]
 return B
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

# Performs line contractions for each qubit combines the lines and then traces out the final state for each qubit. 
function Contract_Lines(Q,H)
 N=size(Q,1)
 A=[]
 C=copy(Q)
 for i=1:N
   push!(A,line_mps(C[i],H))
 end
 B=A[1]*A[2]
 for i=3:N
   B=B*A[i]
 end
#  index=inds(B)
#  #println(index)
#  #println(length(index))
#  for j=1:N
#   A[j]=Par_Trac(B, index[j])
#   #println("o= ",order(A[j]))
#  end
 return B

end

#Performs a nodal contraction by depth and traces out the final state for each qubit.
function Contract_Node(Q,H,E,verbose,depth)
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
   println("D=",D)
   println("N=",N)
   println("E_L=",length(Et))
   println("Q_L=", length(Bh))
   println("E=",Et)
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
      Bh[j]=Edge_contract(Bh[j],H,Et[a[i][1]],n)
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
