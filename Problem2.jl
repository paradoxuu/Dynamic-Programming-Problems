using Plots
using BenchmarkTools

function vfsolvex(vnew, kgrid, tolerance, imax, σ=1.5)
    β=0.9
  
    v=vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i=1
  
    c = kgrid .- kgrid'
    c[c .< 0] .= 0
    if σ==1
      u = log.(c)
    else
      u = (c.^(1-σ).-1)./(1-σ)
    end
    u[kgrid .- kgrid' .< 0] .= -Inf
  
    while maximum(abs.(v-vnew)) > tolerance && i <= imax
      v = vnew;
  
      (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
      i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v=vnew, kprime=kgrid[kprimeindex], kprimeindex=kprimeindex)
  end
  
  
  kupper = 5
  klower = 0.01
  n = 500
  kgrid = collect(range(klower, stop=kupper, length = n))
  (v, kprime, kprimeindex) = vfsolvex(zeros(n), kgrid, 0.001, 1000)
  #(b)
    kgrid = collect(range(klower, stop=kupper, length=500))
    (v, kprime, kprimeindex) = vfsolvex(zeros(500), kgrid, 0.001, 1000);
    plot(kgrid, v, label = "v, σ=1.5")
  
  #(c)
  plot()
  for x in [1.25, 1.1, 1.001, 1.00001, 1]
    kgrid = collect(range(klower, stop=kupper, length=500))
    (v, kprime, kprimeindex) = vfsolvex(zeros(500), kgrid, 0.001, n, x);
    display(plot!(kgrid, v, label = "v, σ=$x"))
  end
f(x, β)= (1/(1-β))log(1-β) + (β/(1-β)^2)log(β) + (1/(1-β))log(x)
g(x) = f(x, 0.9)
plot!(a, g, kgrid, label = "vlog")

(d)
function policy(x)
    e=getindex(findall(kgrid .== x),1)
    return kprime[e]  
end 

function findindex(x)
    a=findmin(abs.(x .- kgrid), dims=1)
    return a[2]
end 

function findkpath(T, k_0)
    w = findindex(k_0)
    kpath = zeros(T)
    k = kgrid[w]
        for i in 1:T
            k = policy(k)
            kpath[i] = k
        end 
    time = collect(range(1, T, length = T))
    plot(time, kpath, label = "kpath")
end 
    
findkpath(100, 5)

