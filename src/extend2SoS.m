function [v] = extend2SoS(v)
%extend2SoS extends a given polynomial to make it a sum-of-squares (thus PD).
%
%  Given a polynomial (with no constant or linear term and SPD quadratic term)
%  this function returns a higher degree polynomial where the additional
%  terms guarantee positive (semi-)definiteness of the polynomial by ensuring
%  it is a sum-of-squares.
%
%  Usage:  [v] = extend2SoS(v)
%  
%  The input describes a polynomial of the form:
%     p(x) = v{2}.'*kron(x,x) + v{3}.'*kron(kron(x,x),x) + ... v{degree}.'*...
%
%  The output is a higher degree polynomial that is written as
%     pbar(x) = (vbar{1}.'*x + vbar{2}.'*kron(x,x) + ... )^2
% 
%  where vbar{k}.' is n-by-n^k and pbar(x) and p(x) have the same coefficients
%  up to degree=degree.  The polynomial pbar is expanded and the coefficients
%  of the 2*(degree-1) polynomial are returned.
%
%  Authors:  Hamza Adjerid and Jeff Borggaard
%
%  License:  MIT
%
%  Part of the KroneckerTools repository.
%%

  degree = length(v);           % the degree of the original polynomial
  n      = sqrt(length(v{2}));

  vbar = cell(1,degree-1);      % number of intermediate variables

  vec = @(X) X(:);              % define the vec function for convenience

  vbar{1} = chol( reshape(v{2},n,n), 'lower' );

  for k=2:degree-1
    rhs = reshape(v{k+1},n,n^k);
    for j=2:k-1
      rhs = rhs - reshape(vec(v{j}*v{k+1-j}.'),n,n^(j+1));
    end
    vbar{k} = 0.5*vbar{1}\rhs;
  end

  % form the new monomial terms to make the polynomial v a sum-of-squares.
  if (degree==3)
    v{4} = vec(vbar{2}.'*vbar{2});

  elseif (degree==4)
    v{5} = vec(vbar{3}*vbar{2}.') + vec(vbar{2}*vbar{3}.');
    v{6} = vec(vbar{3}*vbar{3}.');
  
  else
    warning('yet to be implemented in general')
 
  end

  %for d=degree+1:2*(degree-1)
  %  v{d} = vec(v{})

end
