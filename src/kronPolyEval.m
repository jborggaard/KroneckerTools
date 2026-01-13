function [x] = kronPolyEval(c,z,degree)
%kronPolyEval evaluates a Kronecker polynomial out to degree=degree
%
%  Usage:
%     x = kronPolyEval(c,z,degree);
%
%  Input Variables:
%     c - coefficients of the polynomial (cell array)
%     z - value to calculate the polynomial at (vector)
%     degree - polynomial will be evaluated out to degree=degree
%              (default is the length of c)
%
%  Returns:  x = c{1}*z + c{2}*kron(z,z) + ... c{degree}*kron(...)
%
%  or optionally,
%
%            x = c{1}.'*z + c{2}.'*kron(z,z) + ...
%
%  We optionally handle the special case c{1} = [] by starting with the
%  quadratic term for applications such as energy functions or value functions.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository.
%%
  
  d = length(c);
  if (nargin<3)
    degree = d;
  else
    if (d<degree)
      error('kronPolyEval: not enough entries in the coefficient cell array')
    end
  end
  
  %% Transpose the coefficients if required
  % Check dimensions of the last entry in c to see which dimension is the
  % appropriate size. To multiply, the second dimension of c{k} has to be
  % length(z)^k; instead of checking all of them, just check the last
  % entry, with k=d=length(c) already having been defined.
  [nRows, nCols] = size(c{d});
  
  if (length(z)^d == nCols)
      % No need to transpose
  elseif (length(z)^d == nRows)
      % Need to transpose, otherwise dimensions will not work
      c = cellfun(@transpose,c,'UniformOutput',false); 
  else
      % Probably won't ever happen but good contingency
      warning('Dimensions of polynomial coefficients are not consistent')
  end
  
  %% Perform polynomial evalauation

  %  Special case if the linear term is an empty cell
  if isempty(c{1})
    n = size(c{2},1);
    x = zeros(n,1);
  else
    x = c{1}*z;
  end
  
  zk = z;
  for k=2:degree
    zk = kron(zk,z);
    x  = x + c{k}*zk;
  end
  
end % function kronPolyEval
