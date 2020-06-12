using Optim, NLsolve

function Himmelblau(in::Array{Float64,1})
    x = in[1]
    y = in[2]

    out = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
end


function g(G, in::Array{Float64,1})
    x = in[1]
    y = in[2]

    G[1] = (4*x)*(x^2 + y - 11) + 2*(x + y^2 - 7)
    G[2] = 2*(x^2 + y - 11) + (4*y)*(x + y^2 - 7)
end

guess = [4.0, 2.0]
opt_newton = @elapsed optimize(Himmelblau, guess, Newton())
opt_nlsolve = @elapsed nlsolve(g, guess)
