module openQASM

using ITensors
using Random
using DelimitedFiles

include("component_def.jl")
include("parcing.jl")
include("tensor_fun.jl")
include("wave.jl")
include("wave_comp.jl")
include("noteb.jl")


##################################################################
#                                                                #
# Structures to hold the data for gate function and registrars   #
#                                                                #
##################################################################

mutable struct registar
   name
   l
   qorc
   init
end

mutable struct Name_gate
    name 
    paramaters
    inputs
    gate
    Name_gate(v) = (name = v, paramaters = [], inputs = [], gate = [])
end


##################################################################
#                                                                #
# Basic function for sorting through input files                 #
#                                                                #
##################################################################

function findchar(a,s)
   i=1
    while (a[i]!=s) && i<length(a)
        i=i+1
    end
    if a[i] == s
        return i-1
    else
        return i
    end
end

function findChArray(a,s)
    go=true
    j=1
    while go && j<length(a)
           if findchar(a[j],s)<length(a[j]) 
                go=false
           else
                j=j+1
           end
    end
    return j
end


function splitequal(a,st = ';')
    s="0"
    n=length(a)
    #println(isa(a,Array))
    if n>1 && isa(a,Array)
        m=findChArray(a,'=')
        if findchar(a[m],'=')+2<length(a[m])
            s=a[m][findchar(a[m],'=')+2:findchar(a[m],st)]
        else
            if m < length(a)
                s=a[m+1][1:findchar(a[m+1],st)]
            end
        end
    else
      if isa(a,Array)
        m=findchar(a,'=') 
        if m+2<length(a[1])
            s=a[1][m+2:findchar(a[1],st)] 
        end
      else
        m=findchar(a,'=') 
        if m+2<length(a)
            s=a[m+2:findchar(a,st)] 
        end
      end
    end
    return s
end

function get_contained_string(a,n,target=',')
    n1=findchar(a[n+2:end],target)+n+1
    s=a[n+2:n1]
    
   return s,n1 
end

function split_onchar(a,s)
             n=findchar(a,s) 
             s=a[1:n]
             ch=a[n+2:end]
    return s, ch
end

function compress_string(a)
   b=""
   for i=1:length(a) 
        b=b*a[i]
   end
   return b
end

function ispattern(a,s)
   n=findfirst(s, a)
   if isnothing(n)
     return false     
   else
     return true
   end
end

function findStArray2D(a,s)
    go=true
    j=1
    while go && j<length(a)
           #println(j)
           for i=1:length(a[j])
               if a[j][i]==s || ispattern(a[j][i],s)
                go=false
               else
                j=j+1
               end
           end
    end
    return j
end

function findStArray1D(a,s)
    go=true
    j=1
    while go && j<length(a)
           #println(j)
               if a[j]==s || ispattern(a[j],s)
                go=false
               else
                j=j+1
               end
           
    end
    return j
end

##################################################################
#                                                                #
# Functions to prepare data                                      #
#                                                                #
##################################################################

function splitlines(a)
    #println(a)
    i=1
    #println(i)
    m=length(a)
    #println(m)
    while i < m+1
        n=findChArray(a[i],';')
        #println(n)
        if  n < length(a[i])
            b=a[i][1:n]
            c=a[i][n+1:end]
            #println(b)
            #println(c)
            deleteat!(a[i], n+1:length(a[i]))
            #println(a[i])
            insert!(a,i+1,c)
            #println(a[i+1])
            m=length(a)
        end
        i=i+1
    end
    
    return a
    
end

function remove_comment(a)
    i=1
    n=length(a)
    while i<n+1
        #println(i)
        if length(a[i][1])>1
            if a[i][1]=="//" || (a[i][1][1:2]=="//")
               deleteat!(a, i) 
               n=length(a)
               i=i-1
            else
                if a[i][1]=="/*" || (a[i][1][1:2]=="/*")
                   m=findStArray2D(a[i:end],"*/") 
                   #println(m)
                   for j=1:m
                       deleteat!(a, i)  
                   end
                   i=i-m
                   n=length(a)
                end
            end
        end
        if i>0
            if length(a[i])>1
            m = findStArray1D(a[i],"//")
            if m < length(a[i]) || ispattern(a[i][m],"//")
                if a[i][m][1:2]=="//"
                    deleteat!(a[i],m:length(a[i]))
                else
                    if m+1<length(a[i]) || m+1==length(a[i])
                        deleteat!(a[i],m+1:length(a[i]))
                    end
                    l=findfirst("//", a[i][m])
                    a[i][m]=chop(a[i][m],head=0,tail=length(a[i][m])-l[1]+1)
                    
                end
            end
            end
        end
        if i>0
            if a[i][1] == "" && length(a[i]) ==1
                    deleteat!(a, i)
                    i=i-1
                    n=length(a)
            end
        end
        if i>0
            j=1
            k=length(a[i])
            while j<k+1 && k > 0
                if a[i][j] == ""
                    deleteat!(a[i], j)
                    k=length(a[i])
                else
                    j=j+1
                end
            end            
        end   
        if i>0
            if isempty(a[i])
                    deleteat!(a, i)
                    i=i-1
                    n=length(a)                
            end
        end
        i=i+1
        if i<1
            i=1
        end
        #println(n)
        #println(i)
    end
    return a
end


##################################################################
#                                                                #
# Functions to set up gate and registrars data                   #
#                                                                #
##################################################################

function qregis(a)
    n = findchar(a[2],'[')
    if n < length(a[2]) 
        n1 = findchar(a[2],']')
        m=parse(Int64,a[2][n+2:n1])
        A = registar(a[2][1:n],m,true,'0'^m)
    else
        n = findchar(a[1],'[')
        n1 = findchar(a[1],']')
        m=parse(Int64,a[1][n+2:n1])
        A = registar(a[2][1:end-1],m,true,'0'^m)
    end
    s=splitequal(a)
    if s != "0"
        A.init=s
    end
   return A 
end

function cregis(a)
    n = findchar(a[2],'[')
    if n < length(a[2]) 
        n1 = findchar(a[2],']')
        m=parse(Int64,a[2][n+2:n1])
        A = registar(a[2][1:n],m,false,'0'^m)
    else
        n = findchar(a[1],'[')
        n1 = findchar(a[1],']')
        m=parse(Int64,a[1][n+2:n1])
        A = registar(a[2][1:end-1],m,false,'0'^m)
    end
    s=splitequal(a)
    if s != "0"
        A.init=s
    end
   return A 
end

function add_gate(a)
    m=findchar(a[1][2],'(')
    A=Name_gate(a[1][2][1:m])
    k=findChArray(a[1],')')
    b=compress_string(a[1][2:k])
    n=findchar(b,'(')
    l=length(b)
    go =true
    while go && n<l
      s,n1= get_contained_string(b,n,':')
      if n1<l
            n=n1
      end
      s,n= get_contained_string(b,n,',')
      push!(A.paramaters,s[1:findchar(s,')')]) 
      if findchar(s,')')<length(s)
        go =false
      end
      
    end

    
    if a[1][end]=="{" || findchar(a[1][end],'{') < length(a[1][end])
        n=length(a[1])-1
    else
        n=length(a[1])
    end
    
    for i=k+1:n
       push!(A.inputs,a[1][i]) 
    end
    
    
    i=2
    j=true
    while j==true
        if a[i][end]=="}"
            j=false
            if length(a[i])>1
             push!(A.gate,a[i][1:end-1])
            end
        else
            push!(A.gate,a[i])
            i=i+1
        end
    end
    
    return A
end



##################################################################
#                                                                #
# Functions for working with registar                            #
#                                                                #
##################################################################

function set_init(a,r)
    s=splitequal(a)
    if s != "0"
        r.init=s
    end
end

function give_qbit(a,r)
    k=findchar(a,'[')
    l=findchar(a,']')
    b = a[1:k]
    n=0
    m=1
    while b != r[m].name && m<length(r)
            m = m+1
    end
    if m > 1
        for i=1:m-1
            n=n+r[i].l
        end
    end
    n=n+parse(Int64,a[k+2:l])
    return n
end


##################################################################
#                                                                #
# Functions to set up gate for network                           #
#                                                                #
##################################################################

function get_n(a,r, mod = "none" )
    if a[1][1]=='U'
           i=findChArray(a,')') 
           b=compress_string(a[1:i])
            #println(b)
           m=findchar(b,'(')
           m2=findchar(b,',')
           s1=b[m+2:m2]   
           m3=findchar(b[m2+2:end],',')+m2+1 
           s2 = b[m2+2:m3]
           m4 = findchar(b,')')
           s3 = b[m3+2:m4]
            #println("s1= ",s1)
            #println("s2= ",s2)
            #println("s3= ",s3)
            #println(a[end][1:findchar(a[end],';')])
        if mod == "inv"
            n=[give_qbit(a[end][1:findchar(a[end],';')],r),eval(Meta.parse(s1)),eval(Meta.parse(s2)),eval(Meta.parse(s3))]
        else
            if mod == "ctrl"
                n=[give_qbit(a[end-1][1:findchar(a[end-1],';')],r),give_qbit(a[end][1:findchar(a[end],';')],r),-eval(Meta.parse(s1)),-eval(Meta.parse(s2)),-eval(Meta.parse(s3))]
            else
                n=[give_qbit(a[end][1:findchar(a[end],';')],r),-eval(Meta.parse(s1)),-eval(Meta.parse(s2)),-eval(Meta.parse(s3))]
            end
        end
            
        else
          n=[] 
          for i=2:length(a)
               push!(n, give_qbit(a[i][1:findchar(a[i],';')],r))#give_qbit(a[i][1:findchar(a[i],';')],r))
          end 
        end
    return n
end

function apply_g(a,r,g,gate,qbit)
    if a[1] == "ctrl" || a[1] == "inv" || a[1][1:findchar(a[1],'[')] == "pow"
        for j=1:length(g)
            if a[2] == g[j].name
               if a[1][1:findchar(a[1],'[')] == "pow"
                   i=findChArray(a,']')
                   b=compress_string(a[1:i]) 
                   gate,qbit= apply_ng(a[i:end],r,g[j],gate,qbit,mod= b)
                    
                else
                   gate,qbit= apply_ng(a[2:end],r,g[j],gate,qbit,mod= a[1])
                end
            end
        end
        if a[1][1:findchar(a[1],'[')] == "pow"
           i=findChArray(a,']') 
           b=compress_string(a[1:i])
           m=eval(Meta.parse(b[findchar(b,'[')+2:findchar(b,']')]))
            #println(a[i+1:end])
           for k=1:m
                apply_g(a[i+1:end],r,g,gate,qbit) 
           end
            
        end
        
        if a[1] == "inv"
           push!(gate,a[2][1:findchar(a[2],'(')]) 

           push!(qbit,get_n(a[2:end],r,"inv")) 
        end
        
        if a[1] == "ctrl"
           push!(gate,"c"*a[2][1:findchar(a[2],'(')]) 

           push!(qbit,get_n(a[2:end],r,"ctrl")) 
        end
        
    else
       push!(gate,a[1][1:findchar(a[1],'(')]) 

       push!(qbit,get_n(a,r))
    end
    return gate,qbit
end

function apply_ng(a,r,g,gate,qbit,mod= "none")
  #println(a)
  if mod[1:findchar(mod,'[')] == "pow"
        m=eval(Meta.parse(mod[findchar(mod,'[')+2:findchar(mod,']')]))
        for j=1:m-1
            apply_ng(a,r,g,gate,qbit)
        end
        
  end
  i=findChArray(a,')') 
  #println(i)

  if i < length(a)
    c=compress_string(a[1:i])
    #println(c)
    n=findchar(c,'(')
    #println(n)
    #println(length(c))
    #println(length(g.gate))
    for k=1:length(g.gate)
      #println(g.gate[k])
      if mod =="inv"
          b = copy( g.gate[length(g.gate)-k+1])     
      else
          b = copy( g.gate[k])
      end
      #println(b)
      for j=1:length(g.paramaters)
         #println(g.paramaters[j])
         s,ch=split_onchar(g.paramaters[j],'=')
         if n < length(c)
             s1,n1=get_contained_string(c,n) 
             s1,h=split_onchar(s1,'=')
             if h == ""
                  ch=s1[1:findchar(s1,')')]   
                  n=n1
             else
               if s==s1
                  ch = h[1:findchar(h,')')]
                  n = n1        
               end
                      
             end
         end
         #println(ch)
         if ch != ""
             for l =1:length(b)
                b[l]=replace(b[l],s=>ch)
             end
         else
             println(s," is not defined.")        
         end
      end
      for j=1:length(g.inputs)
         for l=1:length(b)
                 b[l]=replace(b[l],g.inputs[j]=>a[j+i][1:findchar(a[j+i],';')])
         end
      end
      gate,qbit = apply_g(b,r,g,gate,qbit)
      if mod == "inv"
                gate[end]="i"*gate[end]      
      end
    end
    
  else
    for k=1:length(g.gate)
      if mod =="inv"
          b = copy( g.gate[length(g.gate)-k+1])     
      else
          b = copy( g.gate[k])
      end
      if length(g.paramaters)>0
          for j=1:length(g.paramaters)
             #println(g.paramaters[j])
             s,ch=split_onchar(g.paramaters[j],'=')
             #println(s)
             #println(ch)
             for l=1:length(b)
                 b[l]=replace(b[l],s=>ch)
             end
          end        
      end
      for j=1:length(g.inputs)
         for l=1:length(b)
             b[l]=replace(b[l],g.inputs[j]=>a[1+j][1:findchar(a[1+j],';')])
         end
      end      
      gate,qbit = apply_g(b,r,g,gate,qbit)
      if mod == "inv"
                gate[end]="i"*gate[end]      
      end
    end
  end

  return gate,qbit
end

##################################################################
#                                                                #
# Functions to set up network                                    #
#                                                                #
##################################################################


function create_networkinputs(a,c,r,g,v=false)
    gate=[]
    qbit=[]
    #println(length(a))
    i=1
    while i<length(a)+1
        #println(a[i])
        for j=1:length(r)
            if a[i][1][1:findchar(a[i][1],'=')]==r[j].name
                set_init(a[i],r[j])
            end
        end
        for j=1:length(c)
            if a[i][1][1:findchar(a[i][1],'=')]==c[j].name
                set_init(a[i],c[j])
            end
        end
        for j=1:length(g)
            if a[i][1][1:findchar(a[i][1],'(')] == g[j].name
               gate,qbit= apply_ng(a[i],r,g[j],gate,qbit)
            end
        end
        if a[i][1] == "gate"
            go = true
            while go
                for j=1:length(a[i])
                   if findchar(a[i][j],'}')<length(a[i][j])
                        go=false
                   end
                end
                if go
                    i=i+1
                end
            end
        end
        if a[i][1] == "CX" || a[i][1] == "CCX" || a[i][1][1] == 'U' || a[i][1] == "ctrl" || a[i][1] == "inv" || a[i][1][1:findchar(a[i][1],'[')] == "pow"
           gate,qbit= apply_g(a[i],r,g,gate,qbit)
        end
        if v
            println(gate)
            println(qbit)
        end
        i=i+1
    end
    
    
    
    n=r[1].l
    init=r[1].init
    for i=2:length(r)
       n=n+r[i].l 
       init=init*r[i].init
    end
    cn=c[1].l
    cinit=c[1].init
    for i=2:length(c)
       cn=cn+c[i].l 
       cinit=cinit*c[i].init
    end
    if v
        println(init)
        println(cinit)
    end
    
    A=noteb.qc_network(n);
    #println(A)
    noteb.Set_init(A,init);
    #println(A)
    noteb.Set_cinit(A,cinit);
    #println(A)
    noteb.Add_gate(A,gate,qbit);
    return A
end

function InterpitopenQasm(a,v=false)
  qreg=[]
  reg=[]
  g_def=[]
  for i=1:length(a)
        b=a[i][1][1:findchar(a[i][1],'[')]
        #println(b)
        if b == "gate"
                l=i
                j=true
                while j==true
                    if a[l][end]=="}"
                        j=false
                    else
                        l=l+1
                    end
                end
            push!(g_def,add_gate(a[i:l]))
        end
        if b == "qreg"||b == "qubit"
            push!(qreg,qregis(a[i]))
        end
        if b == "creg"||b == "bit"
            push!(reg,cregis(a[i]))
        end
  end
  
  if v  
   println(qreg)
   println(reg)
   println(g_def)
  end
    
  return create_networkinputs(a,reg,qreg,g_def,v)
    
end

##################################################################
#                                                                #
# Call Function                                                  #
#                                                                #
##################################################################

function ReadopenQasm(InPutFile,v=false)
   a = parcing.Read_InPutFile(InPutFile)
   a = remove_comment(a)
   a = splitlines(a)
   if v
       for i=1:length(a)
           println(a[i])
       end
   end
   return InterpitopenQasm(a,v)
end


end
