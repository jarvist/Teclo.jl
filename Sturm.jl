# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

#using Roots, Polynomial

using Calculus

println("Sturm Library...")

N=1000000

# Calculates number of eigenvalues less than 'sigma' in tridiagonal matrix 
# described by: diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1)
# Nb: Off digonal elements E are supplied squared for computational speed
function sturm(D,E,sigma)
    t=0.0
    countnegatives=0
    
    t=D[1]-sigma #first term of sequence, to avoid needing E[0], t[0]
    if t<0.0
        countnegatives=countnegatives+1
    end

    for i=2:length(D)
        t=D[i]-sigma-E[i-1]/t   # Sturm sequence, overwriting temporary files...
        if t<0.0                     # if t<0, we've found another eigenvalue
            countnegatives=countnegatives+1 
        end
    end

    return countnegatives
end

#units eV
kB=8.6173324E-5

function randexp(N) # random exponential with the ln(X) 0<X<1 method
    return (log(rand(N)))
end


function randH(disorder,B,Z,U)
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

