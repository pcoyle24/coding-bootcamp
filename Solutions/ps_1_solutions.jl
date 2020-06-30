########Problem 1
function factorial2(n::Int64)
    val = 1 #initialize value
    for i = 1:n
        val*=i #multiply by all integers between 1 and n
    end
    val #return
end

########Problem 2
function p(x::Float64, coeff::Vector{Float64})
    val = 0 #initializes
    for (i, coef) in enumerate(coeff)
        val += x^(i-1) * coef
    end
    val
end
p(1.0, [2.0, 3.0, 4.0])


########Problem 3
function Monte_Pi(iter::Int64)
    avg = 0
    for i = 1:iter
        x, y = rand(2) #bivariate uniform draw
        dist = sqrt((x - 0.5)^2 + (y - 0.5)^2) #distance from middle of square; pythagorean theorem
        avg+=(dist<0.5)
    end
    avg/=iter #approximation of area
    avg*=4 #π = A / r^2, where r = 1/2. So here π = 4A
    avg
end
approx = Monte_Pi(1000000)

########Problem 4
using Distributions, Plots

function Simulate_Q4(M::Int64)
    a, b, c, d, σ = 0.1, 0.2, 0.5, 1.0, 0.1 #true parameter values
    dist = Normal(0,1)
    ols_vec = Any[]
    x_1, x_2 = rand(dist, 50), rand(dist, 50) #random draws

    #loop over simulations
    for i = 1:M
        w = rand(dist, 50) #randomly draw values of w

        #construct X and Y
        X = hcat(x_1, x_1.^2, x_2, ones(50))
        Y = a .* x_1 .+ b.*x_1.^2 .+ c.*x_2 .+ 1.0 .+ σ .* w
        ols = inv(X' * X) * X' * Y
        push!(ols_vec, ols)
    end
    ols_vec = hcat(ols_vec...)
end

#Simulate and plot
ols_vec = Simulate_Q4(200)
p1 = histogram(ols_vec[1,:], title="Estimates of a")
p2 = histogram(ols_vec[2,:], title="Estimates of b")
p3 = histogram(ols_vec[3,:], title="Estimates of c")
p4 = histogram(ols_vec[4,:], title="Estimates of d")
Plots.plot(p1, p2, p3, p4, layout = (2,2), legend=:none)
Plots.savefig("C:/Users/Garrett/Desktop/ps1_4.png")


########Problem 5
using Distributions, Plots

function Simulate_Q5(n::Int64, α::Float64)
    first_times = zeros(n)
    dist = Normal(0,1)
    σ = 0.2

    for i = 1:n #loop over simulations
        cross, x_t = false, 1
        for t = 1:200
            x_tp = α * x_t + σ * rand(dist)
            if x_tp<=0 #first crossing
                first_times[i] = t
                cross = true
                break
            end
            x_t = x_tp
        end
        if !cross #never crosses
            first_times[i] = 200
        end
    end
    first_times #return
end

#1.0
temp = Simulate_Q5(10000, 1.0)
histogram(temp)
Plots.savefig("C:/Users/Garrett/Desktop/ps1_5a.png")

#0.8
temp = Simulate_Q5(10000, 0.8)
histogram(temp)
Plots.savefig("C:/Users/Garrett/Desktop/ps1_5b.png")

#1.2
temp = Simulate_Q5(10000, 1.2)
histogram(temp)
Plots.savefig("C:/Users/Garrett/Desktop/ps1_5c.png")

########Problem 6
function Newton(f, fp, x_0::Float64, tol::Float64, maxiter::Int64)
    error, root = 100, 0
    while error>tol
        root = x_0 - f(x_0)/fp(x_0)
        error = abs(root - x_0)
        x_0 = root
    end
    root
end

f(x) = (x-1)^3
fp(x) = 3 * (x-1)^2
root = Newton(f, fp, 4.0, 1e-3, 100)

#different function: f(x) = ln(x) + 3x - 7
f(x) = 2*log(x) + 3*x - 7
fp(x) = (2/x) + 3
root = Newton(f, fp, 4.0, 1e-3, 100)

########Problem 7
using Parameters, Plots #read in necessary packages

#a struct to hold model primitives.
@with_kw struct Primitives
    β::Float64 = 0.99 #discount rate.
    θ::Float64 = 0.36 #capital share
    δ::Float64 = 0.025 #capital depreciation
    k_grid::Array{Float64,1} = collect(range(1.0, length = 1000, stop = 45.0)) #capital grid
    nk::Int64 = length(k_grid) #number of capital elements
    markov::Array{Float64,2} = [0.977 0.023; 0.074 0.926] #markov transition process
    z_grid::Array{Float64,1} = [1.25, 0.2] #productivity state grid
    nz::Int64 = length(z_grid) #number of productivity states
end

#mutable struct to hold model results
mutable struct Results
    val_func::Array{Float64,2}
    pol_func::Array{Float64,2}
end

#function that executes the model and returns results
function Solve_model()
    prim = Primitives() #initialize primitives
    val_func = zeros(prim.nk, prim.nz) #preallocate value function as a vector of zeros
    pol_func = zeros(prim.nk, prim.nz) #preallocate value function as a vector of zeros
    res = Results(val_func, pol_func) #initialize results
    V_iterate(prim, res) #value function iteration
    prim, res #return deliverables
end

#value function iteration program. Note that we do EVERYTHING in functions, handling as few global
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-3)
    error = 100 #starting error
    n = 0 #counter
    while error>tol #main convergence loop
        n+=1
        v_next = Bellman(prim, res) #execute Bellman operator
        error = maximum(abs.(v_next - res.val_func)) #reset error term
        res.val_func = v_next #update value function held in results vector
    end
    println("Value functions converged in ", n, " iterations.")
end

#Bellman operator
function Bellman(prim::Primitives, res::Results)
    @unpack β, δ, θ, nz, nk, z_grid, k_grid, markov = prim #unpack parameters from prim struct. Improves readability.
    v_next = zeros(nk, nz)

    for i_k = 1:nk, i_z = 1:nz #loop over state space
        candidate_max = -1e10 #something crappy
        k, z = k_grid[i_k], z_grid[i_z] #convert state indices to state values
        budget = z*k^θ + (1-δ)*k #budget given current state. Doesn't this look nice?

        for i_kp = 1:nk #loop over choice of k_prime
            kp = k_grid[i_kp]
            c = budget - kp #consumption
            if c>0 #check to make sure that consumption is positive
                val = log(c) + β * sum(res.val_func[i_kp,:].*markov[i_z, :])
                if val>candidate_max #check for new max value
                    candidate_max = val
                    res.pol_func[i_k, i_z] = kp #update policy function
                end
            end
        end
        v_next[i_k, i_z] = candidate_max #update next guess of value function
    end
    v_next
end

prim, res = Solve_model() #solve for policy and value functions

#unpack our results and make some plots
@unpack val_func, pol_func = res
@unpack k_grid = prim

#plot value function
Plots.plot(k_grid,val_func, title="Value Functions", label = ["Good State" "Bad State"])
Plots.savefig("C:/Users/Garrett/Desktop/Value_Functions.png")

#plot policy functions
Plots.plot(k_grid, pol_func, title="Policy Functions", label = ["Good State" "Bad State"])
Plots.savefig("C:/Users/Garrett/Desktop/Policy_Functions.png")
