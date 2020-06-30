###show off compilation
function f(a, b)
    y = (a + 8b)^2
    return 7y
end
f(1, 2)
@code_native f(1, 2)  #view machine code generate by Julia

###globals vs not
a = 3.0 #bad global variable
@elapsed for i = 1:1000000
    global a
    a+=i
end

##
function timetest()
    a = 3.0 #good function-contained variable
    for i = 1:1000000
        a+=i
    end
    return a
end

@elapsed timetest() #much faster!

#####showing off profiler
using Profile

function test(x::Float64)
    for i = 1:10000000
        x+=i #do a bunch of things to x
        x-=i
        x*=i
        x/=i
        x = abs(x)
        x = log(x)
    end
    x
end

Profile.clear() #empty the profile
@profile test(2.0) #profile a call of test function
Profile.print() #print profile

####Optimal Savings
using Parameters, Plots

#struct to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.99 #discount factor
    θ::Float64 = 0.36 #production
    δ::Float64 = 0.025 #depreciation
    k_grid::Array{Float64,1} = collect(range(0.1, length = 1800, stop= 45.0)) #capital grid
    nk::Int64 = length(k_grid) #number of capital grid states
end

mutable struct Results
    val_func::Array{Float64,1} #value function
    pol_func::Array{Float64,1} #policy function
end

#function to solve the model
function Solve_model()
    #initialize primitives and results
    prim = Primitives()
    val_func, pol_func = zeros(prim.nk), zeros(prim.nk)
    res = Results(val_func, pol_func)


    error, n = 100, 0
    while error>0.0001 #loop until convergence
        n+=1
        v_next = Bellman(prim, res) #next guess of value function
        error = maximum(abs.(v_next .- res.val_func)) #check for convergence
        res.val_func = v_next #update
        println("Current error: ", error)
    end
    println("Value function converged in ", n, " iterations")
    prim, res
end

#Bellman operator
function Bellman(prim::Primitives, res::Results)
    @unpack β, δ, θ, nk, k_grid = prim #unpack primitive structure
    v_next = zeros(nk) #preallocate next guess of value function

    for i_k = 1:nk #loop over state space
        max_util = -1e10
        k = k_grid[i_k] #value of capital
        budget = k^θ + (1-δ)*k #budget

        for i_kp = 1:nk
            kp = k_grid[i_kp] #value of k'
            c = budget - kp #consumption
            if c>0 #check if postiive
                val = log(c) + β * res.val_func[i_kp] #compute utility
                if val > max_util #wow new max!
                    max_util = val #reset max value
                    res.pol_func[i_k] = kp #update policy function
                end
            end
        end
        v_next[i_k] = max_util #update value function
    end
    v_next
end

@elapsed prim, res = Solve_model() #solve the model.
Plots.plot(prim.k_grid, res.val_func) #plot value function
###
