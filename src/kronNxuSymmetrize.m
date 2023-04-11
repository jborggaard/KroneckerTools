function [Mout] = kronNxuSymmetrize(M,n,d)
%kronNxuSymmetrize  Symmetrizes the columns of a matrix multiplying a
%                   Kronecker monomial kron(x^(d),u), where x^(d) is
%                   the d-times Kronecker product of x.  This is useful for
%                   systems with bilinear and higher-order terms.
%
%  Usage: 
%           Nxu = kronNxuSymmetrize(Nxu,n,d);
%
%  Input Variables:
%     Nxu - a matrix multiplying a Kronecker monomial kron(x^(d),u)
%     n   - length of the vector x
%     d   - degree of the x terms in the multinomial
%
%  Returns: A symmetrized version of Nxu.  Thus, 
%      Nxu*kron(x^(d),u) on input = Nxu*kron(x^(d),u) on output, but the
%      values of Nxu multiplying terms x(i)*x(j)*...*x(l) are all the same
%      for any permutation of (i,j,...,l).
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository.
%%

  [m1,m2] = size(M);
  Mout = zeros(m1,m2);

  %  We achieve the symmetrization by working with all of the columns associated
  %  with a particular entry of the vector "u".  First we find the size of u:
  m = m2/(n^d);

  for j=1:m
    idx = j:m:m2;

    % loop over the rows
    for i=1:m1
      Mout(i,idx) = kronMonomialSymmetrize(M(i,idx),n,d);
    end
  end

end