function [Y] = calTTtimesv(T,m,k,v)
%calTTtimesv Calculates the term \cal{T}_{m,k}.'v
%         w = \|i\|=k kron(T_i1.',kron(T_i2.',kron(T_i3,...,T_im))*v
%
%  Usage:  w = calTTtimesv(T,m,k,B)
%
%  Variables:   T    a cell array with matrices of dimension n times n^p
%               v    a vector of dimension (n^m)
%
%  Contains a subfunction (similar to LyapProduct.m) that computes one term
%  in the above sum.
%
%  This is multiplication performed recursively using Kronecker product rules.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Reference:  Nonlinear Balanced Truncation: Part 2-nonlinear manifold model
%              reduction, Kramer, Gugercin, and Borggaard (arXiv).
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%

  [~ ,n] = size(M);
  [nB,m] = size(B);

  if ( n^2==nB )
    Y = zeros(n^2,m);
    for j=1:m
      Y(:,j) = reshape(M*reshape(B(:,j),n,n)*M.',n^2,1);
    end
    
  else
    Y = zeros(nB,m);
    for j=1:m
      T = reshape(B(:,j),nB/n,n)*M.';
      
      Z = zeros(nB/n,n);
      for k=1:n
        Z(:,k) = kroneckerLeft(M,T(:,k));
      end
      Y(:,j) = Z(:);
    end
  end
  
end

