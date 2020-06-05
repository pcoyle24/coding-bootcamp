# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/04/2020
# lec2.jl
# ------------------------------------------------------------------------------

## Compilation
function f(a,b)
    y = (a + 8b)^2
    return 7y
end

f(1,2)
@code_native f(1,2) # show off assembly code

## Global vars v not (and some macros)
a = 3.0
# run time macro @time or @elapsed
@elapsed for i = 1:1000000
    global a
    a += i
end

function timetest()
    a = 3.0
    for i = 1:1000000
        a += i
    end
    return a
end
@elapsed timetest() # runs like HELLLLLLA faster

# profile macro (where are bottlenecks)
using Profile
function test(x::Float64) #declaring type of var
    for i = 1:100000000
        x += i
        x -= i
        x *= i
        x /= i
        x = abs(x)
        x = log(x)
    end
    return x
end

Profile.clear() #clear profile (initialize)
@profile test(2.0) #macro will store the important info that we want to see
Profile.print() #output will tell us where bottlefunctions are (investigate later)

## A simple economic model
# Optimal Savings

using Parameters, Plots
# Parameters Pkg (to avoid listing EVERY param in input to function)
#with_kw macro allows to give values to primitives
@with_kw struct Primitives
    cBET::Float64 = 0.99
    cTHETA::Float64 = 0.36
    cDEL::Float64 = 0.025
    nk::Int64 = 1800
    k_grid::Array{Float64,1} = collect(range(0.01,length = nk, 45.0))
end

mutable struct Results #allows us to update struct
    val_func::Array{Float64,1}
    pol_func::Array{Float64,1}
end

function solve_model()
    prim = Primitives()
    val_func, pol_func = zeros(prim.nk), zeros(prim.nk)
    res = Results(val_func, pol_func)

    error, n = 100, 0

    while error > 0.0001
        n += 1
        v_next = Bellman(prim, res)

        error = maximum(abs.(v_next .- res.val_func))
        res.val_func = v_next
        #println("Current Error: ", error)
    end
    println("At iteration: ", n)

    return prim, res
end

function Bellman(prim::Primitives, res::Results) #define types of structs
    @unpack cBET, cTHETA, cDEL, nk, k_grid = prim #so we dont need prim.__
    v_next = zeros(nk)

    for i_k = 1:nk
        k = k_grid[i_k]
        max_util = -1e10
        budget = k^cTHETA + (1-cDEL)*k

        for i_kp = 1:nk
            kp = k_grid[i_kp]
            c = budget - kp
            if c > 0
                val = log(c) + cBET*res.val_func[i_kp]
                if val > max_util
                    max_util = val
                    res.pol_func[i_k] = kp
                end
            end
        end
        v_next[i_k] = max_util
    end
    return v_next
end

@elapsed prim,res = solve_model()

Plots.plot(prim.k_grid, res.val_func)
