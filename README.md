# KroneckerTools
Functions for working with Kronecker sum systems and as well as polynomials expressed using Kronecker products.

Used by various functions in both the QQR and NLbalancing software repositories.

# License
This software is covered by the MIT license.

## Installation Notes
Download and untar, or clone this repository:
```
  git clone https://www.github.com/jborggaard/KroneckerTools
```

Optional: 
Get the Matlab functions for efficiently solving linear systems with a special Kronecker sum structure (laplace-like structure) at https://anchp.epfl.ch/index-html/software/misc and detailed in the preprint: _Recursive blocked algorithms for linear systems with Kronecker product structure_, by Minhong Chen and Daniel Kressner.  Place the directory "tensor_recursive" within the main directory.

## Summary of Functions (in the src directory)
```
 [w] = calTTtimesv(T,m,k,v)
```
Given a cell array T containing m matrices where the dimension of T{i} is n-by-n^i, and the dimension of v is n-by-n^m, computes the product
calT_{m,k}.'*v.  The term calT_{m,k} is the sum of Kronecker products of matrices kron(T{i_1},kron(T{i_2},...,T(i_m)) where sum_{j=1}^m i_j = k.

```
 [x] = kronPolyEval(c,z,degree)
```
Evaluates a Kronecker product polynomial out to a given degree.  For example, 

> x = c{1}z + c{2}kron(z,z) + ... + c{degree}kron(z,kron(z,...kron(z,z)).

```
 [Y] = kroneckerLeft(M,B)
```
Given a cell array of matrices M where M(k) has dimension n1(k)-by-n2(k), and matrix B of dimension prod(n2)-by-m, computes the product:

> Y = kron(M{1},kron(M{2},kron(M{3},...))) B.

```
 [Y] = kroneckerRight(B,M)
```
Given a matrix B of dimension m-by-prod(n1), and a cell array of matrices M where M(k) has dimension n1(k)-by-n2(k), computes the product:

> Y = B kron(M{1},kron(M{2},kron(M{3},...))).

```
 [x] = KroneckerSumSolver(A,b,degree,M)
```
Solves the linear system

> (kron(A{d},eye(n^(d-1))) + ... + kron(eye(n^(d-2)),kron(A{2}),eye(n)) + kron(eye(n^(d-1)),A{1}) + diag(M))x = b.

```
 [Mv] = LyapProduct(M,v,d)
```

Computes the product of an N-way Lyapunov matrix with a vector v.  Namely, 

> Mv = ( M⊗I⊗...⊗I + I⊗M⊗I⊗...⊗I + ... + I⊗...⊗I⊗M ) v

where there are d terms in the Kronecker product expressions above, M is an n-by-n matrix, and I=eye(n).

```
 [S] = perfectShuffle(p,q)
```
Creates the perfect shuffle matrix coffesponding to (p,q).  This is a matrix that is useful for Kronecker product permutation operations.  For example, S 
takes a deck of pq cards, splits it into p piles, then produces a new deck by including the top card from each pile in cyclic fashion.  

If A is m-by-n and B is p-by-q, then

> kron(B,A) = perfectShuffle(m,p)*kron(A,B)*perfectShuffle(n,q).'
