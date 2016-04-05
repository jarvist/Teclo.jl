Sturm und Drang / Storm and Urge
=====

Calculation of electronic densities of states for conjugated polymer chains, with off-diagonal disorder parameterised by statistical mechanics.

* Given an arbitrary effective potential energy landscape U=f(\theta). [The Urge]
  * Specified as a functional form ( `U(theta)=( E0 * sin(theta*pi/180.0)^2 ) #P3HT like` )
  * OR read in from tabulated data (i.e. a Quantum-Chemical potential energy scan), fitting a Chebyshev Polynomial with a Vandermonde matrix via the excellent ApproxFun package.
* Integrate (monte carlo direct sampling) to get a (statistical mechanical) partition function, Z=\sum e^{U/kBT}
* Use this partition function to generate random samples of \theta
* Use a model for the transfer integral between two neighbouring units (i.e. monomers in a polymer chain, J~cos(theta)) to build a tridiagonal tight-binding Hamiltonian
* Solve this tridiagonal Hamiltonian with `Sturm sequence` methods, which are linear in time and require no memory. [The Storm]

Method developed in these codes are discussed in these talk slides, but I'm afraid it's pretty incoherent without the talk (and not much better with...).

https://speakerdeck.com/jarvist/2016-03-pvcdt-jarvistmoorefrost-from-atoms-to-solar-cells?slide=82

The only published application of this method we applied it to the P3HT system, treating the P3HT as non-interacting chains. Unfortunately, we found that it didn't agree particularly well with the detailed Molecular Dynamics of the rest of the paper. I suspect this is due to the poor model for the inter-monomer potential (we need an effective potential that includes steric hindrance + entropic effects of the sidechains).

Parameter free calculation of the subgap density of states in poly(3-hexylthiophene)  
Jarvist M. Frost,   James Kirkpatrick,   Thomas Kirchartz and   Jenny Nelson  
Faraday Discuss., 2014,174, 255-266  
http://dx.doi.org/10.1039/C4FD00153B

# Teclo

Extension of these methods to molecular crystals. 

* If tri-diagonal, use linear scaling Sturm sequencies to generate a historgrammed DoS
* If not - currently just standard linear algebra methods. Though perhaps Arrowhead methods for 2D/3D systems in the future?
