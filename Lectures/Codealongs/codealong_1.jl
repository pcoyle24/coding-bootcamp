println("hello world!")

######Variables
#different types
x = 2
y = 4.0
z = "foo"
typeof(x)
typeof(y)
typeof(z)

#basic operations
x*y
x^2
x+y

####Arrays, vectors
#horizontal and vertical vectors
vec_vert = [1, 1.0, "foo"]
vec_hor = [1 1.0 "foo"]
vec_vert[1]
vec_hor[1]

#matrix
array = [1 1.0; "foo" "bar"]
array[1,1]
array[2,2]

#tuples
x = ("foo", "bar", 1)
x[3]
x[3] = 2

#ranges, fancier indexing
x = collect(1:1:100)
y = collect(range(1, length = 100, stop = 100))
@show x[2:50]
x[75:end]

#####For, while loops
#loop over an array
languages = ["julia", "python", "stata"]
for l in languages
    println(l)
end

#iterator
for i = 1:100
    println(i)
end

#while loop
val = 1.0
tol = 0.002
while val > tol
    global val
    val = val / 2
end
val

#####Conditionals
x = 1
x == 1
x == 2
x == 1 && x == 2
x == 1 || x == 2

if x == 2
    println("true")
elseif x!=2
    println("false")
end


#####Making functions
function foo(x; a = 2)
    y = x^2 * a
    x, y
end

x, y = foo(2)
x, y = foo(2, a = 3)

f(x) = x^2
f(4)

#broadcasting
x = collect(1:1:100)
y = log.(x)


######Example 1: White Noise
using Plots, Distributions

dist = Normal(0,1)
n = 1000
ϵ = rand(d, n)
plot(1:n, ϵ)
histogram(ϵ)

#######Example 2: OLS
β_0 = 1.0
β_1 = 2.0
β_2 = 3.0
n = 10000
x = rand(n).*10
x2 = x.^2
ϵ = rand(dist, n)
Y = β_0 .+ β_1.*x + β_2.*x2 .+ ϵ
X = hcat(ones(n), x, x2)
β_ols = inv(X' * X) * X' * Y



#####################
