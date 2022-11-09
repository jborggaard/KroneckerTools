function [c] = kronPolySymmetrize(c,n,degree)
%kronPolySymmetrize Produces symmetric Kronecker polynomial coefficients
%
%  Usage:
%     c = kronPolySymmetrize(c,n,degree);
%
%  Input Variables:
%     c - coefficients of the Kronecker polynomial (vector of coefficients)
%     n - number of variables in the multinomial
%     degree - degree of the multinomial terms (c should have n^degree terms)
%
%  Returns: symmetrized coefficients of the multinomial
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository.
%%
  
  dc2 = size(c,2);

  if (dc2>1)
    isRow = true;
    c = c(:);
  else
    isRow = false;
  end
  
  if (nargin<3)
    degree = round( log(length(c))/log(n) );
  end
  
  % an inefficient way to approach this and requires functions in the util
  % subdirectory
  S = Kron2CT(n,degree);
  C = CT2Kron(n,degree);

  c = C*S*c;
  
  if (isRow)
    c = c';
  end

end % function kronPolySymmetrize
