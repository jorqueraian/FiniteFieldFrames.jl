# FiniteFieldFrames.jl

[![Build Status](https://github.com/jorqueraian/FiniteFieldFrames.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jorqueraian/FiniteFieldFrames.jl/actions/workflows/CI.yml?query=branch%3Amain)


Julia implementation (using Oscar.jl) to work with frames of finite fields
Everything is done on the level of gram matrices.
The Oscar.jl documentation can be found [here](https://docs.oscar-system.org/stable/) with installation instructions [here](https://www.oscar-system.org/install).

Activate and load this package. In the FiniteFieldFrames directory
```julia
pkg> activate .  # Press ']' to enter the Pkg REPL mode.
```
then 
```julia
julia> using FiniteFieldFrames
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
Specify `case="U"` when calling functions.
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

binder_finder(hessa_gram, "U")
```

## Instructions for Case O
Specify finite field, of odd characteristic, using Oscar.
```julia
using Oscar

ff = GF(5,3,"a");
```
Then specify `case="O"` when calling functions.

```julia
ff = GF(3);
Phi = matrix(ff, [
   0 0 0  0  1  1  1  1  1  1;
   0 0 2  1  0  0  1  1  2  2;
   1 2 0  0  0  0  1  2  1  2;
   1 1 2  2  1  2  0  0  0  0     
])
gram = transpose(Phi)*diagonal_matrix([ff(1),ff(1),ff(1),ff(2)])*Phi

binder_finder(3, gram, "O")
```

## Basic Functionality
Better documentation writeup coming soon

To check is a provided Gram matrix `gram` is the gram matrix of a frame in either case "O" or "U" use

```julia
(frame_bool, d) = is_frame(gram, case)
```
Returns a boolean `frame_bool` if the provided matrix is the gram matrix of a frame for some `d`-dimensional space. if `frame_bool=false` then `d=nothing`.

To check if a given Gram matrix `gram` is the gram matrix of a collection of equiangular vectors is case "O" or "U" use
```julia
(equiangular_bool, a, b) = is_equiangular(gram, case)
```
Returns a boolean `equiangular_bool` and parameters `a`, the common magnatudes of the vectors, and `b` the angle. If `equiangular_bool=false` then both `a` and `b` will be `nothing`.

To check if the Gram matrix of a frame `gram`, is the Gram matrix of a tight frame use
```julia
(tight_bool, c) = is_frame_tight(gram)
```
Returns a boolean `equiangular_bool` and parameters `c` where `G^2=cG`.

To check if the matrix `gram` is the Gram matrix of an ETF use
```julia
(etf_bool, a, b, c, d) = is_ETF(gram, case)
```
which calls all of the above functions.

If `gram` is the Gram matrix of a frame, you can use 
```julia
Phi = reconstruct_frame_from_gram(gram, case)
```
which will attempt to construct the corresponding frame. I am not sure I have got this one completly working yet, but I think it is atleast basically correct. In case O, the output `Phi` may live in a field extension of the field in which `gram` was defined in. The output will always be in the real model, so the the corresponding scalar product is just the dot product.


## Constructions 

## Binder Finder
This is not guaranteed to return a complete binder, or maybe it can be proven that it does in some cases, idk. wait I think i did that. TODO: did i do this?


Referneces
https://arxiv.org/pdf/2012.12977
https://arxiv.org/pdf/2012.13642
https://arxiv.org/abs/2505.12175



