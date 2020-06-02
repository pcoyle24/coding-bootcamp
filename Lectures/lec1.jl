# ------------------------------------------------------------------------------
# Author: Philip Coyle
# Date Created: 06/02/2020
# lec1.jl
# ------------------------------------------------------------------------------

## Hello World!
println("Hello World");

## Variables (Types)
x = 2; # An Integer
y = 4.0 # A float
typeof(x)
typeof(y)

z = "foo"; # A String
typeof(z)

## Basic Operations
x*y #multiplies two different types, but outputs a float
x^2

z2 = "bar";
z*z2 #Concatinating strings

## Arrays and Vectors
# Horizontal and vertical vectors
# Vert: with commas
vec_vert = [1., 2, "Foo"] #Types differ
vec_vert = [1., 2., 3.] #Types the same
vec_horz = vec_vert = [1. 2. 3.] #No commas

#identifying elements (must be int)
vec_vert[3]

# Array
array = [1 1.; "foo" "bar" ]
array = [1 1.;
        "foo" "bar" ]

#identifying elements
array[1,2] #requires two ints

## Range and fancier index
X = collect(1:1:100)
Y = collect(range(1,length = 100, stop = 100))

# @show (macro) like println()
@show X[2:50]
@show X[75:end]

## Loops
# For Loops (loops over an array)
languages = ["Julia", "Python", "Stata"]
for l in languages
        println(l)
end

for ll in 1:3
        println(languages[ll])
end

# While Loops
val = 100
tol = 1
while val > tol
        global val # So julia knows its same val above
        val = val/2
end

## Boolean
x = 1
x == 1 # == to ask boolean
x == 2

x == 1 && x == 2 # and
x == 1 || x == 2 # or

if x == 2
        println("foo")
elseif x != 2
        println("bar")
end

## Functions (required function; customizable options)
function foo(x; a = 2)
        y = x^2*a
        return x,y
end

x,y = foo(2)
x,y = foo(2;a = 4)


## Applying to all elements (Broadcast with a period)
x = collect(1:1:100)
y = log(x)
y = log.(x)

## Random number generation and plots
using Plots, Distributions, Random

dist = Normal(0,1) #(μ,σ)
n = 1000
ϵ = rand(dist,n) #randomize given dist above

plot(1:n,ϵ)
histogram(ϵ)

## OLS
β_0  = 1.0
β_1  = 2.0
β_2  = 3.0
n = 10000

X = rand(n)*10 #uniform [0,10]
X2 = X.^2
ϵ = rand(dist,n) #randomize given dist above

Y = β_0*ones(n) + β_1*X + β_2*X2 + ϵ
XX = hcat(ones(n),X,X2) # Design Matrix

β_ols = inv(XX'*XX)*(XX'*Y)



##
