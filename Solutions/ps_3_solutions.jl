############Problem 1: Making Change
function count(coins::Vector{Int64}, n::Int64)
    num = length(coins)
    table = zeros(n+1, num) #base case for n is n = 0

    #only one way to make change for zero
    table[1,:] .= 1

    #fill remainder
    for i = 2:n+1 #loop over target values
        for j = 1:num

            #number of solutions with
            count_with, count_without = 0, 0
            if i - coins[j]>0
                count_with = table[i - coins[j], j] #after including one of coin i, compute numgber of solutions including coin i for REMAINING value
            end

            if j!=1 #if not on last coin
                count_without = table[i, j-1] #number of solutions without coin j
            end
            table[i,j] = count_with + count_without #add together
        end
    end
    return table[n+1, num]
end

coins = [2,5,3,6]
n = 100
count(coins, n)

#########Problem 2: Rod Cutting
function cut_rod(prices)
    n = length(prices)
    val = zeros(n+1) #include length zero, the value of which is assumed to be zero.

    for i = 2:n+1 #+1 #loop over remaining lengths
        max_val = -1e10

        for j = 1:(i-1)
            max_val = max(max_val, prices[j] + val[i - j]) #choose the best option via DP
        end
        val[i] = max_val
    end
    val[n+1] #return
end

prices = [1,5,8,9,10,17,17,20]
cut_rod(prices)

#########Problem 3: 0-1 Knapsack
function Knapsack(weights::Vector{Int64}, vals::Vector{Int64}, C::Int64) #accepts knapsack capacity, weights, and values
    n = length(weights)
    opts = zeros(C+1, n+1)
    opts[1,:].=0 #can't obtain any value with a sack of zero capacity!
    opts[:,1].=0 #can't obtain any value with zero items

    for w = 2:C+1, i = 2:n+1  #begin looping
        opts[w, i] = opts[w, i-1] #default: item i is not used
        if weights[i-1]<=(w-1) #possible to include item i. Have to adjust because slot 1 stands for 0 weight.
            opts[w,i] = max(opts[w, i-1], vals[i-1] + opts[w - weights[i-1], i-1]) #choose the best option via DP
        end
    end

    return opts[C+1, n+1]
end

weights = [5,4,6,3]
vals = [10,40,30,50]
C = 10
Knapsack(weights, vals, C)

weights = [10, 20, 30]
vals = [60, 100, 120]
C = 50
Knapsack(weights, vals, C)









##############
