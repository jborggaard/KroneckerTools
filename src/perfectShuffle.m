function [S] = perfectShuffle(p,q)
%perfectShuffle  Creates the perfect shuffle matrix corresponding to (p,q)
%
%  A matrix that is useful for Kronecker product permutation operations.  
%  For example, S takes a deck of p*q cards, splits it into p piles then
%  produces a new deck by including the top card from each pile in cyclic
%  fashion.  
%
%  See, e.g.,
%    van Loan, The ubiquitous Kronecker product, JCAM (123), pp. 85-100,
%    2000.
%
%  Part of the QQR library.

  r = p*q;
  I = speye(r);
  
  S = sparse(r,r);
  
  for qq=1:q
    S( ((qq-1)*p+1):qq*p, : ) = I(qq:q:r,:);
  end
end

