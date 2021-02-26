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

function splitlines(a)
    i=1
    m=length(a)
    while i < m
        n=findChArray(a[i],';')
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


function splitequal(a)
    n=length(a)
    m=findchar(a[2],'=')
    s="0"
    if m<length(a[2])
        if m<length(a[2])+2
            s=a[2][m+2:end-1]
        else
            s=a[3][1:end-1]
        end
    end
    if n>2
        for i=3:n
           m=findchar(a[i],'=')
           if m<length(a[i])+2
               s=a[i][m+2:end-1] 
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
    b = a(1:k)
    n=0
    m=0
    for i=1:length(r)
        if b == r[i].name
            m = i
        end
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

function get_n(a,r, inv = false )
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
        if !inv
            n=[give_qbit(a[end][1:findchar(a[end],';')],r),eval(Meta.parse(s1)),eval(Meta.parse(s2)),eval(Meta.parse(s3))]
        else
            n=[give_qbit(a[end][1:findchar(a[end],';')],r),-eval(Meta.parse(s1)),-eval(Meta.parse(s2)),eval(Meta.parse(s3))]
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

           push!(qbit,get_n(a[2:end],r,true)) 
        end
        
        if a[1] == "ctrl"
           push!(gate,"c"*a[2][1:findchar(a[2],'(')]) 

           push!(qbit,get_n(a[2:end],r)) 
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
      println(g.gate[k])
      if mod =="inv"
          b = copy( g.gate[length(g.gate)-k+1])     
      else
          b = copy( g.gate[k])
      end
      for j=1:length(g.paramaters)
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
         println(ch)
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
      if length(g.paramaters)>1
          for j=1:length(g.paramaters)
             s,ch=split_onchar(g.paramaters[j],'=')
             println(s)
             println(ch)
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


end
