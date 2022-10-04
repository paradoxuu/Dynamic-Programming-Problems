using SymPy
function solveram(α, g_1, g_2)
c_1, c_2, λ, μ, l = symbols("c_1 c_2 λ μ l")
vars = (c_1, c_2, λ, μ, l)
system = [c_1 - μ^-1, c_2 - μ^-1, λ * (α/(1-l)^2) - μ + (α/(1-l)), l - 2/(2+α), c_1 + c_2 - l + g_1 + g_2]
a = nonlinsolve(system, vars)
return a
end 

solveram(1, 0.05, 0.15)