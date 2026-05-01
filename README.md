# FiniteFieldFrames.jl

[![Build Status](https://github.com/jorqueraian/FiniteFieldFrames.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jorqueraian/FiniteFieldFrames.jl/actions/workflows/CI.yml?query=branch%3Amain)


Julia implementation (using Oscar.jl) to work with frames over finite fields.
Everything is done on the level of gram matrices.
The Oscar.jl documentation can be found [here](https://docs.oscar-system.org/stable/) with installation instructions [here](https://www.oscar-system.org/install).

Activate and load this package: In the FiniteFieldFrames directory
```julia
pkg> activate .
pkg> instantiate  # Press ']' to enter the Pkg REPL mode.
```
Or by running `using Pkg; Pkg.activate("."); Pkg.instantiate();`.

Then 
```julia
julia> using FiniteFieldFrames; using Oscar
```

## Instructions for Case U
Pick base field using Oscar(Oscar is loaded by FiniteFieldFrames.jl but not imported). Examples below:
```julia
using Oscar
base_f = GF(5,3,"b");
base_f = GF(5);
```
Build quadradic extension
```julia
Kx, x = base_f["x"];
# specify degree 2 irreducible polynomial.
ff, a = finite_field(x^2+x+1, "a");
```
Specify `case=:U` when calling functions.
An example is the following:
```julia
base_f = GF(5);
Kx, x = base_f["x"];
ff, a = finite_field(x^2+x+1, "a");
hessa_sic = matrix(ff, [
        1    1       1    -1 -1  -1    0   0      0;
        0    0       0     1  a  a^2  -1 -1*a   -1*a^2;
       -1  -1*a^2  -1*a    0  0   0    1  a^2     a
]);
hessa_gram = conjugate_transpose(hessa_sic)*hessa_sic;

binder_finder(hessa_gram, case=:U)
```

## Instructions for Case O
Specify finite field, of odd characteristic, using Oscar.
```julia
using Oscar

ff = GF(5,3,"a");
```
Then specify `case=:O` when calling functions.

```julia
ff = GF(3);
Phi = matrix(ff, [
   0 0 0  0  1  1  1  1  1  1;
   0 0 2  1  0  0  1  1  2  2;
   1 2 0  0  0  0  1  2  1  2;
   1 1 2  2  1  2  0  0  0  0     
])
gram = transpose(Phi)*diagonal_matrix([ff(1),ff(1),ff(1),ff(2)])*Phi

binder_finder(3, gram, case=:O)
```

## Basic Functionality
To check is a provided Gram matrix `gram` is the gram matrix of a frame in either case :O or :U use

```julia
(frame_bool, d) = is_frame(gram; case=case)
```
Returns a boolean `frame_bool` if the provided matrix is the gram matrix of a frame for some `d`-dimensional space. if `frame_bool=false` then `d=nothing`.

To check if a given Gram matrix `gram` is the gram matrix of a collection of equiangular vectors is case :O or :U use
```julia
(equiangular_bool, a, b) = is_equiangular(gram; case=case)
```
Returns a boolean `equiangular_bool` and parameters `a`, the common magnatudes of the vectors, and `b` the angle. If `equiangular_bool=false` then both `a` and `b` will be `nothing`.

To check if the Gram matrix of a frame `gram`, is the Gram matrix of a tight frame use
```julia
(tight_bool, c) = is_frame_tight(gram)
```
Returns a boolean `equiangular_bool` and parameters `c` where `G^2=cG`.

To check if the matrix `gram` is the Gram matrix of an ETF use
```julia
(etf_bool, a, b, c, d) = is_ETF(gram; case=case)
```
which calls all of the above functions.

If `gram` is the Gram matrix of a frame, you can use 
```julia
Phi = reconstruct_frame_from_gram(gram; case=case)
```
which will attempt to construct the corresponding frame. I am not sure I have got this one completly working yet, but I think it is atleast basically correct. In case O, the output `Phi` may live in a field extension of the field in which `gram` was defined in. The output will always be in the real model, so the the corresponding scalar product is just the dot product.


## Constructions 
The following would construct a 12x78 ETF in Case O for a finite field of characteristic 5.
```julia
gram = etf_from_triangular_graph(12, 5);
```

The following would construct a 41x82 ETF in Case O for a finite field of characteristic 7.
```julia
gram = real_etf_to_case_O(real_dx2d_etf_from_prime_power(3^4), Int((3^4+1)/2), 7);
```

The following would construct a 7x14 ETF in Case U over the finite field of 27^2 elements.
```julia
Phi = etf_from_pmod_diff_set([0,4,6,7,8,11,13], 14, 27);
gram = conjugate_transpose(Phi)*Phi;
```
or the following produces a 7x18 ETF in case U over the finite field of 125^2 elements
```julia
k=2
p=3*k-1 # 5
q=p^3
# sorting not needed
D = sort([[x for x in 0:3:(9k-1)]; 1])
Phi = etf_from_pmod_diff_set(D, 9k, q);
gram = conjugate_transpose(Phi)*Phi;
```
This method does not verify that D is a p-modular difference set.



## Binder Finder
Binder Finder finds the binder of a frame in Case O or U.


## Referneces
[OSCAR] OSCAR -- Open Source Computer Algebra Research system, Version 1.7.2, The OSCAR Team, 2026. (https://www.oscar-system.org)

[OSCAR-book] Wolfram Decker, Christian Eder, Claus Fieker, Max Horn, Michael Joswig, eds. The Computer Algebra System OSCAR: Algorithms and Examples, Algorithms and Computation in Mathematics, Springer, 2025. (https://link.springer.com/book/9783031621260)

[1] Greaves, G., Iverson, J., Jasper, J., & Mixon, D. (2022). Frames over finite fields: basic theory and equiangular lines in unitary geometry. Finite Fields Appl., 77, Paper No. 101954, 41.

[2] Greaves, G., Iverson, J., Jasper, J., & Mixon, D. (2022). Frames over finite fields: equiangular lines in orthogonal geometry. Linear Algebra Appl., 639, 50–80.

[3] Iverson, J., King, E., & Mixon, D. (2021). A note on tight projective 2-designs. J. Combin. Des., 29(12), 809–832.

[4] Joseph W. Iverson, & Dustin G. Mixon. (2025). Asymmetric SICs over finite fields. 

[5] Ian Jorquera, & Emily J. King. (2025). On the Structure of Frames and Equiangular Lines over Finite Fields and their Connections to Design Theory. 

[6] Fickus, M., Jasper, J., King, E., & Mixon, D. (2018). Equiangular tight frames that contain regular simplices. Linear Algebra Appl., 555, 98–138.