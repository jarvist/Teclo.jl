# Julia code to calculate spectrum of tridiagonal 'DoS' TB matrices with Sturm sequences
# See section 5, p38: http://arxiv.org/pdf/math-ph/0510043v2.pdf

#using Roots, Polynomial

using Calculus

require("Sturm")

println("Sturm und Drang: DoS by Sturm sequences")

N=1000000

outfile=open("DoS.dat","w+")
onsetfile=open("onset.dat","w+")
potfile=open("potential.dat","w+")

P=0.0
disorder=0.0

# OK; Raos FF paper (Moreno et al. J.Phys.Chem.B 2010)
# 'full' potential energy Fig 4.a. puts a barrier at 90 degress of ~3.0 kCal / mol = 126 meV
E0=0.126

@sync @parallel for T=200.0:10:400.0 #T=100.0:100:400 #:0.1:1
    B=1/(T*kB) #300K * k_B in eV
    U(theta)=( E0 * sin(theta)^2 ) #P3HT like
    Z=integrate(theta -> exp(-U(theta)*B),-pi,pi, :monte_carlo ) # recalculate Z now that P is changing
    println("Partition function for Z(E0=",E0,")=",Z)

# Following checks the partition function code, outputting p(robability) as a fn(theta) for varying P
    for t = -pi:(pi/180.0):pi
#    println("Partition function Z=",Z)
        p=exp(-U(t)*B)/Z
#        println(t," ",p)
#        @printf(potfile,"%f %f %f\n",t,U(t),p)
    end
e
    
    D,E=randH(disorder,B,Z,U)
    @printf("Calculated with P= %f Z= %f\n",P,Z)
#println("STURM sequence method...")
    sigma=4.0
    pDoS=0
    pold=0
    onset=false
    for sigma=6:0.001:6.4 #sigma=3:0.010:7 #sigma=6:0.001:6.4 #6 to 6.4 covers right lobe
        pDoS=sturm(D,E,sigma)
        @printf("%f %f %f %f\n", T, sigma , pDoS, pDoS-pold)
#        @printf(outfile,"%f %f %f %f\n",T,sigma,pDoS, pDoS-pold)
        if (pDoS-pold>10.0 && onset==false)
#            @printf(onsetfile,"%f %f %f\n",P,sigma,pDoS-pold)
            onset=true
        end
        pold=pDoS
    end
end

close(outfile)
close(onsetfile)

end

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
