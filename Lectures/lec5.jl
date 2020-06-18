# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/16/2020
# lec5.jl
# ------------------------------------------------------------------------------

## Fibonacci numbers, 2 ways

# Counting Upwards
function fib1(n::Int64)
    fib_vec = zeros(n+1)
    fib_vec[2] = 1

    for i = 3:n+1
        fib_vec[i] = fib_vec[i-1] + fib_vec[i-2]
    end

    return fib_vec[n+1]
end

fib1(8)

# Recursive Method
function fib2(n::Int64)
    if n <= 1
        return n
    else
        return fib2(n-1) + fib2(n-2)
    end
end

fib2(8)

# Which one is faster?
@elapsed fib1(45)
@elapsed fib2(45)

# Why is fib2 slower?
# fib2(4) = fib2(3) + fib2(2)
# fib2(3) = fib2(2) + fib2(1) + fib2(2)
# fib2(2) = fib2(1) + fib2(0) + fib2(2)
# .... (we have to run fib2(2) TWICE .... grows if we have even larger numbers)
# It requires  A LOT of function calls
@elapsed fib1(1000000)
@elapsed fib2(1000000) #impossible....the universie will end first (lol)

# fib1() requires no additional function calls. so it can run quickly. It builds up.

## Egg dropping problem
# Recursive formulation
function egg_drop_rec(n::Int64, k::Int64)
    if k <= 1 || n == 1
        return k
    else
        candidate_min = 1e100 # a bad number
        for j = 1:k # Looping over possible floors
            val = max(egg_drop_rec(n-1, j-1), egg_drop_rec(n, k-j))
            if val < candidate_min
                candidate_min = val
            end
        end
    end
    return candidate_min + 1 #To account for the initial drop from floor j
end

# This is super slow because it requires A TON of function calls.

# Dynamic Programming
function egg_drop_dp(n::Int64, k::Int64)
    trials = zeros(n,k) #stores number of trials for any combos of n eggs and k floors
    # Fill in easy elements first and then build off that

    # The trivial solitions
    trials[:,1] .= 1 #Only 1 floor
    # for i = 1:k
    #     trials[1,i] = i
    # end
    trials[1,:] = 1:k #Only 1 egg

    # Looping over higher eggs/floors
    for i_n = 2:n
        for i_k = 2:k
            candidate_min = 1e100

            for i_j = 1:i_k # loop up to number of floors you are considering...what happens if eggg does/doesn not break
                # Eggs Breaks
                val_1 = 0 # Number of additioal trials if egg breaks
                if i_j > 1
                    val_1 = trials[i_n-1, i_j-1]
                end

                # Egg Survives
                val_2 = 0 # Number of additioal trials if egg survives
                if i_j < i_k
                    val_2 = trials[i_n, i_k - i_j]
                end

                res = max(val_1,val_2) + 1
                if res < candidate_min
                    candidate_min = res
                end
            end

            trials[i_n, i_k] =  candidate_min
        end
    end

    return trials[n,k]
end



egg_drop_dp(2,36)
# How to solve in 8 trials?
# If drop at 18, then floor 1-17 to go.
# First drop should be at floor 8.
# If egg breaks, you have 7 floors to look at.  # Adding up by how many floors you maybe have left.
# If egg does not break, you have 36-8. Next drop is floor 15
# If egg breaks(9-14)
# If egg does not break, 21
# If egg does not break, 26
# If egg does not break, 30
# If egg does not break, 33
# If egg does not break, 35
# If egg does not break, 36



## Matrix chain multiplication
# Express chain as a vector.
# Example: vec = [10, 20, 30, 40] corresponds to 3 matricies [(10x20) (20X30) (30x40)]
# In general, ith matrix is of dim vec[1]xvec[1+i]

#What is optimal placement of parethesis?
function  Matrix_multiply(mat_vec::Vector{Int64,1})
    num_mat = length(mat_vec) - 1
    mult_mat = zeros(n,n) # Fill in elements that correspond to easier sols and build off to find harder sols
    # mult_mat[i,j] = mimimum number of oeprations to compute matrix i through matrix j
    # dimension of mat i is mat_vec[i]xmat_vec[i+1]

    # Trivial Cases
    for i = 1:num_mat
        mult_mat[i,i] = 0
    end

    # Loop over length of chains
    # Consider mat A, B, C, D. How to compute whole thing the fastest?
    # Check pairs: (AB), (CD), (BC)...
    # From pairs, we can then figure out triplets of matricies....and build out from there.

    # First: Compute how costly it is to multily all pairs
    # Second: figure out how to multiply any trio
    # .... continue building out until the end.

    for L = 2:num_mat
        # CONSIDER ALL PAIRS. loop over rows
        for i = 1:(n - L + 1)
            j = i + L - 1 #How far we are going...

            candidate_min = 1e10
            for k = i:(j-1)#Where to place a () in current chain of matricies.
                cost = mult_mat[i,k] + mult_mat[k+1,j] + mat_vec[i]*mat_vec[k+1]*mat_vec[j+1]
                if cost < candidate_min
                    candidate_min = cost
                end
            end
            mult_mat[i,j] = candidate_min
        end

    end
    mult_mat[i,n]
end
