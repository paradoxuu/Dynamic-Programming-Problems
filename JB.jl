θ = (0.38)^(1/0.6)
α = 0.4
s = 0.0225
u_T = zeros(20)
u = 0.1
for i in 1:20
    
    u = (1-θ^(1-α)) * u + s * (1-u)
    u_T[i] = u
end 
u_T
u_v = vcat([0.1], u_T)
time = collect(range(1, 21, length = 21))
using Plots
plot(time, u_v, label = "unemployment rate time series")

