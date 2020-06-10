# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/09/2020
# lec3.jl
# ------------------------------------------------------------------------------

using Optim, Plots, Distributions, Random

## one dimension
# Brends Method
f(x) = (x-1)^2
opt = optimize(f, -5.0, 5.0) #min and max
opt.minimum #min
opt.minimizer #argim

# only optimize over one var
f(x,y) = (x-y)^2
opt = optimize(x-> f(x,1.0), -5.0, 5.0) #fixing y = 1.0 (arrow)
opt.minimum #min
opt.minimizer #argim


## two-dimension (Rosebrock)
function Rosenbrock(vec:: Vector{Float64})
    x1, x2 = vec[1], vec[2]

    val = (1 - x1)^2 + 100*(x2 - x1^2)^2
    return val
end

#eval function
x_grid = collect(-3:0.01:3)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:nx
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Rosenbrock(guess)
end

Plots.surface(x_grid, x_grid, z_grid, camera = (50, 50))
Plots.contourf(x_grid, x_grid, log.(1 .+ z_grid)) #filled contor

 #lets optimze
 guess = [0.0, 0.0]
 opt = optimize(Rosenbrock, guess)
 guess_toofar = [1000000000.0, 10000000.0]
 opt = optimize(Rosenbrock, guess_toofar)

opt = optimize(Rosenbrock, guess, LBFGS()) #telling it to run LBFGS

# Newtons method (somethings wrong here ... fix later)
function g(G, vec:: Vector{Float64})
    x1, x2 = vec[1], vec[2]

    G[1] = -2*(1-x1) - 400*(x2 - x1^2)*x1
    G[2] = 200*(x2-x1^2)

    return G
end

function h(H, vec:: Vector{Float64})
    x1, x2 = vec[1], vec[2]

    H[1,1] = 2 - 400*x2 + 1200*x1^2
    H[1,2] = -400*x1
    H[2,1] = -400*x1
    H[2,2] = 200

    return H
end

guess = [0.0, 0.0]
opt = optimize(Rosenbrock, g, h, guess)

## Many local min function (Greiwank)
function Greiwank(x::Vector{Float64})
    val = (1/400).*sum(x.^2) + 1 - prod((cos.(x./sqrt(length(x)))))

    return val
end

#eval function
x_grid = collect(-5:0.01:5)
nx = length(x_grid)
z_grid = zeros(nx, nx)

for i = 1:nx, j = 1:nx
    guess = [x_grid[i], x_grid[j]]
    z_grid[i,j] = Greiwank(guess)
end

Plots.surface(x_grid, x_grid, z_grid, camera = (50, 50))
Plots.contourf(x_grid, x_grid, log.(1 .+ z_grid)) #filled contor


guess_init = [3.0, 3.0]
opt = optimize(Greiwank, guess_init)

guess = [2.0, 2.0]
opt = optimize(Greiwank, guess)

# checkig many diff starting points
function multi_start()
    x_grid = collect(-5:0.5:5)
    nx = length(x_grid)
    min, argmin = 100, [100, 100]


    for i = 1:nx, j = 1:nx
        guess = [x_grid[i], x_grid[j]]
        opt = optimize(Greiwank, guess)

        if opt.minimum < min
            min = opt.minimum
            argmin = opt.minimizer
        end
    end

    return min, argmin
end

min, argmin = multi_start()

## OLS with Nelder Mead
dist = Normal(0,1)
cBET0 = 1.0
cBET1 = 2.0
cBET2 = 3.0

n = 100000
x = rand(n).*10
x2 = x.^2

Random.seed!(1234)
cEPS = rand(dist,n)

Y = cBET0  .+ cBET1.*x .+ cBET2.*x2 + cEPS
X = hcat(ones(n), x, x2)
cBET_ols = inv(X'*X)*(X'*Y)

function OLS_NM(cBET)
    Random.seed!(1234)

    cBET0, cBET1, cBET2 = cBET[1], cBET[2], cBET[3]
    cEPS = rand(dist,n)
    Y_true = 1.0 .+ 2.0.*x .+ 3.0.*x2 + cEPS
    Y_predict = cBET0 .+ cBET1.*x .+ cBET2.*x2

    error = sum((Y_true .- Y_predict).^2)
    return error
end

guess_init = [1.0, 2.0, 3.0]
opt = optimize(OLS_NM, guess_init)
