function [Ms] = kronMatrixSymmetrize(M,n,d)
%kronMatrixSymmetrize  Symmetrizes the columns of a matrix multiplying a
%                      Kronecker monomial (kron(x,kron(x,...,x)).

  [m1,m2] = size(M);
  Ms = zeros(m1,m2);

  % get the degree of the monomial term if it isn't provided
  if (nargin<3)
    d = round(log(m2)/log(n));
  end

  % loop over the rows
  for m=1:m1
    Ms(m,:) = kronMonomialSymmetrize(M(m,:),n,d);
  end

end