

module wave_comp

using ITensors

include("component_def.jl")
include("wave_comp.jl")

function Index_setupW(a)

    d=Index[]
    for i=1:a
        append!(d,[Index(2),Index(2)])
    end
    d=reshape(d,(a,2))
    return d

end

function Inital_StateW(N,d,a)
 Q=ITensor[]
 for i=1:N
  push!(Q,wave_comp.Init_stW(d[i,2],a[1+i]))
 end
 return Q
end

function gate_setW(s,n,d)


   if s == "H"

       Ham=wave_comp.HGateW(d[n,1],d[n,2]) 
     
   end
   
   if s == "X"

       Ham=wave_comp.XGateW(d[n,1],d[n,2]) 
     
   end
   
   if s == "CN" || s == "CX"

       Ham=wave_comp.CNotGateW(d[n[1],1],d[n[1],2],d[n[2],1],d[n[2],2]) 
     
   end

 
   if s == "Y"

       Ham=wave_comp.YGateW(d[n,1],d[n,2]) 
     
   end
   
   if s == "Z"

       Ham=wave_comp.ZGateW(d[n,1],d[n,2]) 
     
   end
   
   if s == "I"

       Ham=wave_comp.ZGateW(d[n,1],d[n,2]) 
     
   end
   
   if s == "Rx"

       Ham=wave_comp.RxGateW(d[Int64(n[1]),1],d[Int64(n[1]),2],n[2]) 
     
   end
   
   if s == "Ry"

       Ham=wave_comp.RyGateW(d[Int64(n[1]),1],d[Int64(n[1]),2],n[2]) 
     
   end
   
   if s == "Rz"

       Ham=wave_comp.RzGateW(d[Int64(n[1]),1],d[Int64(n[1]),2],n[2]) 
     
   end

   if s == "SW"

       Ham=wave_comp.SWGateW(d[n[1],1],d[n[1],2],d[n[2],1],d[n[2],2]) 
     
   end 
   
   if s == "CZ"

       Ham=wave_comp.CZGateW(d[n[1],1],d[n[1],2],d[n[2],1],d[n[2],2]) 
     
   end

   if s == "CP"

       Ham=wave_comp.CPhaseGateW(d[Int64(n[1]),1],d[Int64(n[1]),2],d[Int64(n[2]),1],d[Int64(n[2]),2],n[3]) 
     
   end
   
   if s == "CR"

       Ham=wave_comp.CRGateW(d[n[1],1],d[n[1],2],d[n[2],1],d[n[2],2],n[3]) 
     
   end
   
 return Ham
end


end
