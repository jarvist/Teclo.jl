# Sturm_Drang.jl
# Sturm and Urge - Calculation of electronic densities of states for conjugated polymer chains, with off-diagonal disorder parameterised by statistical mechanics.
# https://github.com/jarvist/Teclo

push!(LOAD_PATH,"./") # Temporary version of modules in PWD
using Sturm 
using ApproxFunFromFile

#Bolztmann's constant in units of eV - thereby all the potentials (of functional form or tabulated data) in units eV
kB=8.6173324E-5

N=10^6 # Length of tridiagonal Hamiltonian constructed; changes the statistics

println("Sturm und Drang: DoS by Sturm sequences")

# Outputfiles, streamed to with C-like @printf
DoSfile=open("DoS.dat","w+")
@printf(DoSfile,"#T sigma Cumulative-DoS DoS-this-bin\n#Trying plotting like: gnuplot> p \"DoS.dat\" u 2:4 \n")
onsetfile=open("onset.dat","w+")
potfile=open("potential.dat","w+")
@printf(potfile,"# Try plotting like: gnuplot> p \"potential.dat\" u 1:2, \"\" u 1:(0.1*\$3) w l \n")

P=0.0 # No idea what this does... only appears in line 65 in a printf ?
disorder=0.0 # Energetic disorder, Gaussian form, for the site energies of the polymer

# Size of potential from Raos P3HT ForceField paper (Moreno et al. J.Phys.Chem.B 2010)
# 'full' potential energy Fig 4.a. puts a barrier at 90 degress of ~3.0 kCal / mol = 126 meV
E0=0.126

# The @sync @parallel magic makes this forloop execute in parallel, at least on Linux.
@sync @parallel for T=200.0:10:400.0 #T=100.0:100:400 #:0.1:1
    B=1/(T*kB) #300K * k_B in eV
    #U(theta)=( E0 * sin(theta*pi/180.0)^2 ) #P3HT like
    #U(theta)=( E0 * (-sin(theta*pi/180.0)^2 - sin(2*theta*pi/180.0)^2 ) ) # PFO like
    U,raw=ApproxFunVandermonde("INDT-modred-eV.dat",25) # Use Vandermonde interpolation to load an ApproxFun function
    # See Figure 5.7, Page 213: https://dx.doi.org/10.6084/m9.figshare.91370.v1
   
    Z,epsilon=quadgk(theta -> exp(-U(theta)*B),0.0,360.0) # Now using Julia language (>0.4) built in quadgk numeric integration
    Z=Z/30 # to improve statistics on rejection sampling of thetas 
    # Calculation partition function Z; particular to this potential energy surface + temperature 
    println("Integration of Z via quadgk method, estimated upper bound on absolute error: ",epsilon)
    println("Partition function for Z(E0=",E0,",T=",T,")=",Z)
    

# Following checks the partition function code, outputting p(robability) as a fn(theta) for varying P
    @printf(potfile,"# Potential and probability density at T=%f",T)
    @printf(potfile,"# theta U(theta) probability(theta)")
    for theta = 0.0:1.0:360.0
#    println("Partition function Z=",Z)
        p=exp(-U(theta)*B)/Z
#        println(t," ",p)
        @printf(potfile,"%f %f %f\n",theta,U(theta),p)
    end
    
    # generates separate (D)iagonal and (E)-offdiagonal terms of Tight Binding Hamiltonian
    D,E=randH(5.0,disorder,B,Z,U,N)
    @printf("Calculated with E0= %f Z= %f\n",E0,Z)
#println("STURM sequence method...")
    pDoS=0
    pold=0
    onset=false
    for sigma=3.5:0.01:6.5 # Energy window (eV) over which to bin eigenvalues
        # Sturm sequence returns number of eigenvalues of Hamiltonian {D,E} below sigma
        pDoS=sturm(D,E,sigma)
        @printf(".") # Progress indicator - marching dots ('......')
        @printf(DoSfile,"%f %f %f %f\n",T,sigma,pDoS, pDoS-pold)
        if (pDoS-pold>10.0 && onset==false)
            @printf(onsetfile,"%f %f %f\n",P,sigma,pDoS-pold)
            onset=true
        end
        pold=pDoS
    end
    @printf("\n") # New line after all DoS calculated
end

close(DoSfile)
close(onsetfile)
close(potfile)

end

# This code was used for debugging + profiling

#println("Elements...(offdiag^2, diag))");
#println([E;D])

#println("Full square matrix H");
H=diagm(E.^0.5,-1)+diagm(D)+diagm(E.^0.5,1) #build full matrix from diagonal elements; for comparison
println(H)

println("Eigenvalues")
println(eigvals(H))
println("Min Eigenvalue")
println(eigmin(H))

println("I time a single histogram (Es<sigma) point:")
@time sturm(D,E,5.0)
println("I time eigvals(H)")
@time eigvals(H) 
