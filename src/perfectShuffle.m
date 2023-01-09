function [S] = perfectShuffle(p,q)
%perfectShuffle  Creates the perfect shuffle matrix corresponding to (p,q)
%
%  A matrix that is useful for Kronecker product permutation operations.  
%  For example, S takes a deck of p*q cards, splits it into p piles (with q
%  cards in each pile) then produces a new deck (from the top down) by 
%  including the top card from each pile in cyclic fashion.  A perfect
%  shuffle matrix is a permutation matrix, so shares those usual identities
%  as well as some additional properties shown in the example below. 
%
%  See, e.g.,
%    van Loan, The ubiquitous Kronecker product, JCAM (123), pp. 85-100,
%    2000.
%
%  Examples using perfect shuffle matrices (verifying identities):
%
%    m1 = 3; n1 = 4; m2 = 2; n2 = 5;
%    A = rand(m1,n1); B = rand(m2,n2);
%    Smm = perfectShuffle(m1,m2);  Snn = perfectShuffle(n1,n2);
%    norm(kron(A,B) - Smm.'*kron(B,A)*Snn) % should be zero
%
%    vec = @(x) x(:);    % create a convenience function
%    S = perfectShuffle(m1,n1);
%    norm(S.'*vec(A) - vec(A.')) % should be zero
%
%  Part of the KroneckerTools library.
%%

  r = p*q;
  I = speye(r);
  
  S = zeros(r,r);
  
  for qq=1:q
    S( ((qq-1)*p+1):qq*p, : ) = I(qq:q:r,:);
  end

  S = sparse(S);
end

