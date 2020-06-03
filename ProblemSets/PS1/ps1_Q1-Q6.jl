# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 05/26/2020
# ps1_Q1-Q6.jl
# ------------------------------------------------------------------------------

# clearconsole()
dir = "/Users/philipcoyle/Documents/School/University_of_Wisconsin/SecondYear/Summer_2020/CodingBootcamp/ProblemSets/PS1/"

## Question 1
function factorial2(n)
    nfac = 1;
    if n > 0
        for i in 1:n
            nfac = nfac*i
        end
    end
    return nfac
end
fac_7 = factorial2(7)

## Question 2
function p(x,coeff)
    sum = 0;
    for (index,value) in enumerate(coeff)
        sum = sum + value*(x^(index-1))
    end
    return sum
end

Σ = p(5,[-3., 1., 4.])

## Question 3
function π_mc(sim)
    x = rand(sim);
    y = rand(sim);

    z2 = x.^2 + y.^2;
    count = length(z2[z2 .<= 1.])

    π_apx = 4*(count/sim);
    return π_apx
end

π_apx = π_mc(10000000)


## Question 4
using LinearAlgebra
using Plots

# Housekeeping
a = 0.1;
b = 0.2;
c = 0.5;
d = 1.0;
σ = 0.1;
N = 50;
sim = 20;

# Main Code
β_ols = zeros(4,sim);
for i in 1:20
    # Random Draws for X and X
    R = randn(N,3);
    x_1 = R[:,1];
    x_2 = R[:,2];
    w = R[:,3];

    y = a*x_1 + b*(x_1.^2) + c*x_2 + ones(N,1)*d + σ*w;

    X = hcat(x_1, x_1.^2, x_2, ones(N,1));
    X_tr = transpose(X);
    β_ols[:,i] = inv(X_tr*X)*(X_tr*y);
end

# Plotting Histograms
sp1 = histogram(β_ols[1,:],title="a", bins = 5);
sp2 = histogram(β_ols[2,:],title="b", bins = 5);
sp3 = histogram(β_ols[3,:],title="c", bins = 5);
sp4 = histogram(β_ols[4,:],title="d", bins = 4);
H₁ = plot(sp1,sp2,sp3,sp4,layout=(2,2),legend=false)
savefig(dir*"Q4_Hist.pdf")


## Question 5
# Functions
function first_time_passage(α,a,T)
    ϵ = randn(200,1);
    X = zeros(T+1,1);
    for t in 1:T-1
        X[t+1] = α*X[t] + σ*ϵ[t]
    end
    X_neg = (LinearIndices(X))[findall(x->x < a, X)];
    if isempty(X_neg)
        out = 0;
    else
        out = X_neg[1];
    end

    return out
end

# Housekeeping
α_grid = [0.8, 1.0, 1.2];
a = 0;
σ = 0.2;
T = 200;
sim = 100;
T_0 = zeros(length(α_grid),sim);

# Main Code
for (i,α) in enumerate(α_grid)
    #println(i)
    #println(α)
    for s in 1:sim
        T_0[i,s] = first_time_passage(α,a,T)
    end
end

# Plotting Histograms
α_08 = histogram(T_0[1,:],title="α = 0.8", bins = 5);
α_1 = histogram(T_0[2,:],title="α = 1", bins = 5);
α_12 = histogram(T_0[3,:],title="α = 1.2", bins = 5);
H₂ = plot(α_08,α_1,α_12,layout=(1,3),legend=false)
xlabel!("T_0")
savefig(dir*"Q5_Hist.pdf")

## Question 6
function newtons_method(func, fprime, x_in, tol, maxit)
    it = 1;
    converged = 0;

    x = x_in;
    while converged == 0 && it < maxit
        x_up = x - func(x)/fprime(x);
        global x_up

        diff = broadcast(abs,x_up - x)
        if diff < tol
            converged = 1;
        end
        x = x_up;
    end

    x_out = x_up;
    return x_out
end

# Housekeeping
tol = 1e-10;
maxit = 10000;
x_0 = 2;

# Main Code
f_1(x) = (x-1)^3;
fpr_1(x) = 3(x-1)^2;
xroot_1 = newtons_method(f_1, fpr_1, x_0, tol, maxit)

f_2(x) = x^3 - 4;
fpr_2(x) = 3x^2;
xroot_2 = newtons_method(f_2, fpr_2, x_0, tol, maxit)
