# FiniteFieldFrames
Julia implementation (using Oscar.jl) to work with frames of finite fields
Everything is done on the level of gram matrices.

## Instructions for Case U
Pick base field using Oscar. Examples below:
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
Specify finite field, of odd characteristic.
```julia
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



