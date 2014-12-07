# Teclo.jl: Julia code to calculate eigenvalue spectrum (density of states) for disordered molecular crystal 

println("Teclo: DoS of Molecular Crystals")
println("Long goes the night
Longer the day
Teclo your death
Will send me to my grave")

# Jiahao Chen's Random Matrices package: https://github.com/jiahao/RandomMatrices.jl
#using(RandomMatrices)

# Distributions package: http://distributionsjl.readthedocs.org/en/latest/univariate.html
using(Distributions)

# Normally distributed distro...
gaussian_x=Normal(0.0, 0.05)
println(rand(gaussian_x,20))

# Physical Constants
#units eV
kB=8.6173324E-5

# Simulation Constants
N=10 # Number of sites
E0=0.126  # Potential energy curve
T=300     # Temperature, Kelvin
disorder=0.05 #Gaussian disorder... units eV?

# Logfiles
outfile=open("DoS.dat","w+")
potfile=open("potential.dat","w+")

B=1/(T*kB) #300K * k_B in eV

# Useful functions
function randexp(N) # random exponential with the ln(X) 0<X<1 method
    return (log(rand(N)))
end

# Harmonic energy...
U(theta)=( E0 * sin(theta)^2 ) #P3HT like
Z=1 # false partition function?

function randH(disorder,B,Z,U) # Old long-snake-moan version!
# Random Trace / diagonal elements
    D=5.0 + disorder*randn(N)
# Random Off-diag elements
#E=0.1+0.05*randn(N-1)

#Generate thetas...
#    thetas=randexp(N-1)
    thetas=Float64[]
    for i=1:N-1 #number of thetas we need for off-diagonal element
        theta=0.0
        while true # this is a do-while loop, Julia styleee
            theta=2*pi*rand()     #random theta
            p=exp(-U(theta)*B)/Z  #probability by stat mech
            p>rand() && break     #rejection sampling of distribution
        end
        push!(thetas,theta)       #tack theta onto end of list of thetas
    end

#Transfer integral from
    J0=0.800 #Max Transfer integral
    E=(J0*cos(thetas)).^2
# Squared (element wise) for the Sturm algorithm...
    E= E.^2
    return (D,E)
end

D,E=randH(disorder,B,Z,U)

println("Full square matrix H");
H=diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1) #build full matrix from diagonal elements; for comparison
println(H)

println("Eigenvalues")
println(eigvals(H))
println("Min Eigenvalue")
println(eigmin(H))

println("I time eigvals(H)")
@time eigvals(H) 
