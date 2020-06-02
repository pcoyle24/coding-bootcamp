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
fac_7 = factorial2(7);
println(fac_7)

## Question 2
function p(x,coeff)
    sum = 0;
    for (index,value) in enumerate(coeff)
        sum = sum + value*(x^(index-1))
    end
    return sum
end

Σ = p(5,[-3., 1., 4.]);
println(Σ)

## Question 3
function π_mc(sim)
    x = rand(sim);
    y = rand(sim);

    z² = x.^2 + y.^2;
    count = length(z²[z² .<= 1.])

    π_apx = 4*(count/sim);
    return π_apx
end

π_apx = π_mc(10000000);
println(π_apx)


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
β̂ = zeros(4,sim);
for i in 1:20
    # Random Draws for X and X
    R = randn(N,3);
    x₁ = R[:,1];
    x₂ = R[:,2];
    w = R[:,3];

    y = a*x₁ + b*(x₁.^2) + c*x₂ + ones(N,1)*d + σ*w;

    X = hcat(x₁, x₁.^2, x₂, ones(N,1));
    Xᵀ = transpose(X);
    β̂[:,i] = inv(Xᵀ*X)*(Xᵀ*y);
end

# Plotting Histograms
sp1 = histogram(β̂[1,:],title="a", bins = 5);
sp2 = histogram(β̂[2,:],title="b", bins = 5);
sp3 = histogram(β̂[3,:],title="c", bins = 5);
sp4 = histogram(β̂[4,:],title="d", bins = 4);
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
T₀ = zeros(length(α_grid),sim);

# Main Code
for (i,α) in enumerate(α_grid)
    #println(i)
    #println(α)
    for s in 1:sim
        T₀[i,s] = first_time_passage(α,a,T)
    end
end

# Plotting Histograms
α_08 = histogram(T₀[1,:],title="α = 0.8", bins = 5);
α_1 = histogram(T₀[2,:],title="α = 1", bins = 5);
α_12 = histogram(T₀[3,:],title="α = 1.2", bins = 5);
H₂ = plot(α_08,α_1,α_12,layout=(1,3),legend=false)
xlabel!("T₀")
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
x₀ = 2;

# Main Code
f₁(x) = (x-1)^3;
fpr₁(x) = 3(x-1)^2;
xroot₁ = newtons_method(f₁, fpr₁, x₀, tol, maxit);
println(xroot₁)

f₂(x) = x^3 - 4;
fpr₂(x) = 3x^2;
xroot₂ = newtons_method(f₂, fpr₂, x₀, tol, maxit);
println(xroot₂)
