using NLsolve

β = 0.995
y_H = 1
y_L = 0.75
z_H = 0.3
z_L = 0
ϕ = 0.55
s = 0.02
c = 1
p = 0.35 

function f!(F, v)
    ω_H = v[1]
    ω_L = v[2]
    θ = v[3]
    J_H = v[4]
    J_L = v[5]
    F[1] = β * (0.5*(p/θ)) * (J_H + J_L) - c
    F[2] = (y_H -ω_H) / (1-β*(1-s)) - J_H
    F[3] = (y_L -ω_L) / (1-β*(1-s)) - J_L
    F[4] = ϕ*y_H + (1-ϕ)*z_H + β*p*ϕ*J_H -ω_H
    F[5] = ϕ*y_L + (1-ϕ)*z_L + β*p*ϕ*J_L -ω_L
end

k = [1 0.5 0.5 1.5 1.5]'
soln = nlsolve(f!, k, autodiff = :forward)

@show soln.zero