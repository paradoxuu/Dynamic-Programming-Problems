using Plots
using BenchmarkTools


#Problem 1(a)
function vfsolve(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex, M=i)
end


kupper = 5
klower = 0.05
n = 500
kgrid = collect(range(klower, stop = kupper, length = n))
(v, kprime, kprimeindex, M) = vfsolve(zeros(n), kgrid, 0.001, 1000);
M
#M=95 the output value function is close enough to the true value function
#1(b)
#define a new function that take imax as the only control for the number of iterations
function vfsolve1(vnew, kgrid, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 1.5

    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while i <= imax

        (vnew, cartesianindex) = findmax(u .+ β*vnew', dims = 2);
        i += 1;

    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

(v1, kprime1, kprimeindex1) = vfsolve1(zeros(n), kgrid, 1);
(v5, kprime5, kprimeindex5) = vfsolve1(zeros(n), kgrid, 5);
(v10, kprime10, kprimeindex10) = vfsolve1(zeros(n), kgrid, 10);
plot(kgrid, v1, label = "v1")
plot!(kgrid, v5, label = "v5")
plot!(kgrid, v10, label = "v10")
plot!(kgrid, v, label = "v")

#1(c)
function findindex(x)
    a=findmin(abs.(x .- kgrid), dims=1)
    return a[2]
end 


function policy(x)
    e=getindex(findall(kgrid .== x),1)
    return kprime[e]  
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

findkpath(100, 0.5)

#1(d)
findkpath(100, 4.5)

#(e)

#(e)(a)
function vfsolve_E(vnew, kgrid, tolerance, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 0.5

    v = vnew .+ 2*tolerance
    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while maximum(abs.(v - vnew)) > tolerance && i <= imax
        v = vnew;
        
        (vnew, cartesianindex) = findmax(u .+ β*v', dims = 2);
        i += 1;
    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex, n=i)
end
kupper = 5
klower = 0.05
n = 500
kgrid = collect(range(klower, stop = kupper, length = n))
(v_E, kprime_E, kprimeindex_E, n_E) = vfsolve_E(zeros(n), kgrid, 0.001, 1000)

#(e)(b)
function vfsolve1_E(vnew, kgrid, imax)
    α = 0.37
    δ = 0.1
    β = 0.95
    σ = 0.5

    cartesianindex = Array{CartesianIndex{2}, length(v)}
    i = 1
    c = kgrid.^α + (1-δ)*kgrid .- kgrid'
    c[c .< 0] .= 0
    u = (c .^ (1-σ) .- 1) ./ (1-σ)
    u[kgrid.^α + (1-δ)*kgrid .- kgrid' .< 0] .= -Inf

    while i <= imax

        (vnew, cartesianindex) = findmax(u .+ β*vnew', dims = 2);
        i += 1;

    end
    kprimeindex = getindex.(cartesianindex, 2)
    return (v = vnew, kprime = kgrid[kprimeindex], kprimeindex = kprimeindex)
end

(v1_E, kprime1_E, kprimeindex1_E) = vfsolve1_E(zeros(n), kgrid, 1);
(v5_E, kprime5_E, kprimeindex5_E) = vfsolve1_E(zeros(n), kgrid, 5);
(v10_E, kprime10_E, kprimeindex10_E) = vfsolve1_E(zeros(n), kgrid, 10);
plot(kgrid, v1_E, label = "v1_E")
plot!(kgrid, v5_E, label = "v5_E")
plot!(kgrid, v10_E, label = "v10_E")
plot!(kgrid, v_E, label = "v_E")

#(e)(c)
function policy_E(x)
    e=getindex(findall(kgrid .== x),1)
    return kprime_E[e]  
end 

function findkpath_E(T, k_0)
    w = findindex(k_0)
    kpath = zeros(T)
    k = kgrid[w]
        for i in 1:T
            k = policy_E(k)
            kpath[i] = k
        end 
    time = collect(range(1, T, length = T))
    plot(time, kpath, label = "kpath")
end 

findkpath_E(100, 0.5)

#1(d)
findkpath_E(100, 4.5)