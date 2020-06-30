#########Dynamic programming examples
###Fibonacci numbers
#counting up
function fib1(n::Int64)
    fib_vec = zeros(n+1)
    fib_vec[2] = 1

    for i = 3:n+1
        fib_vec[i] = fib_vec[i-1] + fib_vec[i-2]
    end
    return fib_vec[n+1]
end

#simple implementation
#fib1(1), fib1(0)
#fib1(2) = fib1(1) + fib1(0)
#fib1(3) = fib1(2) + fib1(1)
#...

##recursion
function fib2(n::Int64)
    if n <= 1
        return n
    else
        return fib2(n-1) + fib2(n-2)
    end
end

#fib2(4)
#fib2(3) + fib2(2)
#fib2(2) + fib2(1) + fib2(1) + fib2(0)
#fib2(1) + fib2(0) + fib2(1) + fib2(1) + fib2(0) = 3

#fib2(6)
#fib2(5) + fib2(4)
#fib2(4) + fib2(3) + fib2(4) uh oh!

@elapsed fib1(10000000)
#@elapsed fib2(10000000) don't even try to run this.


#####Egg Dropping

#recursive
function egg_drop_rec(n::Int64, k::Int64)

    if k == 1 || k == 0 || n == 1
        return k
    else
        candidate_min = 1e10 #something bad
        for x = 1:k #loop over floors to drop from
            val = max(egg_drop_rec(n-1, x-1), egg_drop_rec(n, k-x))
            if val<candidate_min
                candidate_min = val
            end
        end
        return candidate_min + 1
    end
end


####Dynamic Programming
function egg_drop_dp(n::Int64, k::Int64)
    trials = zeros(n, k) #matrix that stores number of trials required for any combination of eggs/floors up through and including [n,k]

    trials[:,1] .= 1
    for i = 1:k
        trials[1,i] = i
    end

    for i_n = 2:n, i_k = 2:k
        candidate_min = 1e10

        for i_x = 1:i_k

            val_1 = 0 # Number of additional trials if egg breaks. If floor = 1, then done
            if i_x>1
                val_1 = trials[i_n-1, i_x-1]
            end

            val_2 = 0 #number of additional trials if egg does not breka. If floor = i_k, then done
            if i_x<i_k
                val_2 = trials[i_n, i_k - i_x]
            end

            res = 1 + max(val_1, val_2)
            if res<candidate_min
                candidate_min = res
            end
        end
        trials[i_n, i_k] = candidate_min
    end
    trials[n,k]
end

egg_drop_dp(2, 36)

####How???
#First trial: drop from floor 8. If egg breaks, then 1 egg and 7 flooors
#If not break: drop from floor 15. If breaks, check floors 9-14 with one egg.
#If not break: floor 21
#If not break: floor 26
#If not break: floor 30
#If not break: floor 33
#If not break: floor 35
#If not break: floor 36

#########Matrix Chain
#vector of lengths. ith matrix is of dimension p[i]*p[i+1]

function Matrix_multiply(mat_vec::Vector{Int64})
    n = length(mat_vec)-1 #eg five elements correspnds to 4 matrices and so forth.
    mult_mat = zeros(n, n)

    #mult_mat[i,j]: minimum number to compute A[i] * . . .  A[j]
    #dimension of A[i] is mat_vec[i] * mat_vec[i+1]
    for i = 1:n
        mult_mat[i,i] = 0
    end

    #What we essentially do here is loop over placement of parenthesis.
    #EX: if we have matrices ABCD, then our options for first placement outer side are
    #A(BCD), (AB)CD, or (ABC)D

    #we first loop over lengths of possible chains. Each time, we consider all the ways
    #to multiply from matrix i to matrix j. If j-i =1, then there is only one way.

    ####BASIC idea
    #First we compute how costly it is to multiply all pairs of matrices
    #From this, we can figure out how to optimally multiply any trio of matrices
    #and so on and so forth. Build up until we're done.

    for L = 2:n #loop over chain lengths
        for i = 1:(n-L+1) #going down rows
            j = i + L -1  #how far we're going to
            mult_mat[i,j] = Inf #bad starting values
            for k = i:(j-1) #possible choices of splitting up the vector
                cost = mult_mat[i,k] + mult_mat[k+1,j] + mat_vec[i]*mat_vec[k+1]*mat_vec[j+1]
                if cost<mult_mat[i,j]
                    mult_mat[i,j] = cost
                end
            end
        end
    end
    mult_mat[1,n] #return minimum cost to go ALL THE WAY
end

mat_vec = [1,2,3,4]
Matrix_multiply(mat_vec)

mat_vec = [40, 20, 30, 10, 30]
Matrix_multiply(mat_vec)

#####
