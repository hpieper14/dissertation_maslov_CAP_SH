This directory contains the source code that generates the results presented [here](results). The
first 3 modules are as follow and contain the results discussed in Part 1 of the dissertation:  
* [Fourier Series](source/@FourierSeries): contains code for a `FourierSeries` object 
and operations.
* [Pulse Solution](source/@PulseSolution): contains code for a `PulseSolution`
object. This approximates the pulse solution to the Swift-Hohenberg equation 
using a Fourier series and Newton's method. 
* [Conjugate Points](source/@ConjugatePoints): contains code for a `ConjugatePoints` 
object. This computes the conjugate points associated to the pulse solution. 


The following modules contain code that accompanies Part 2 of the dissertation. 
* [Bistable Bundles](source/BistableBundles): contains code to compute the resonant and non-resonant vector bundles of the bistable equation discussed in Chapter 11.2. 
* [Invariant Manifolds](source/InvariantManifolds): contains code to compute the invariant stable and unstable manifolds for the Swift-Hohenberg equation. This directory also contains a computer assisted proof of their existence. This is discussed in Chapter 9.1. 
* [L minus](source/Lminus): contains code to compute the bound on $L_-$ discussed in Chapters 8.1 and 10.2. 
* [Pseudo Arclength Continuation](source/PseudoArclengthContinuation): contains code to compute branches of pulse solutions for the 1 D Swift-Hohenberg equation. This was exploratory work and is not included in the dissertation. 
* [Pulse Validation](source/PulseValidation): contains code to compute a stationary pulse solution of the Swift-Hohenberg equation and a computer assisted proof of its existence. This is discussed in Chapter 9. 
* [Sequence Space](source/SequenceSpace): contains helper functions for operations in sequence space. 
* [Swift-Hohenberg Bundles](source/SwiftHohenbergBundles): contains code to compute the resonant and non-resonant solutions of the Swift-Hohenberg equation as discussed in Chapter 11.3. This directory is self-contained. See [README.md](source/SwiftHohenbergBundles/README.md).
* [Vector Field](source/VectorField): contains helper functions for expressing the Swift-Hohenberg equation, its vector field, and linearization about the origin. 