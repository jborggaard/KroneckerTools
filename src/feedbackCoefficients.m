function [k] = feedbackCoefficients(v,R,B,g)
%feedbackCoefficients Calculates the coefficients of a feedback law.
%  Given coefficients of an energy function 
%
%     v(x) = 0.5*( v{2}*kron(x,x) + ... + v{d+1}*kron(kron(...((x,x),...x) )
%
%  this function provides the coefficients of a feedback law
%
%     k(x) = -R^{-1}(B+g(x))^T(\nabla v(x))^T
%          = k{1}*x + k{2}*kron(x,x) + ... k{d}*kron(kron...),x)
%
%  Usage:
%    [k] = feedbackCoefficients(v,R,B,g);
%
%  Variables:
%    v  - a cell array containing the coefficients of an energy function given
%         as a Kronecker polynomial
%    R  - the control weighting matrix
%    B  - the input matrix (constant portion of the control inputs)
%    g  - coefficients of input terms that are degree one or higher
%         (not implemented currently)
%
%  Author:
%    Jeff Borggaard, Virginia Tech
%
%  License:
%    MIT
%
%  Part of the KroneckerTools repository
%%
  if (nargin==4)
    warning('degree 1 and higher input terms are not yet implemented')
  end

  d = length(v)-1;
  n = size(B,1);

  k = cell(1,d);
  for i=1:d
    k{i} = -((i+1)/2)*(R\B.'*reshape(v{i+1},n,n^i));
  end

end % function feedbackCoefficients