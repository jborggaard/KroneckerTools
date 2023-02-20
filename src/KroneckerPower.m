function [vd] = KroneckerPower(v,d)
%KroneckerPower Recursively computes the d^th power Kronecker product of v
%
%     Given a vector (or matrix) v, compute it's d^th power Kronecker product 
%
%   Usage:  [vd] = KroneckerPower(v,d);
%
%     vd = (v \otimes v \otimes ... \otimes v)
%           |--           d terms         --| 
%                          
%  Author: Ali Bouland, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%

  % Check inputs
  validateattributes(d,{'numeric'},{'positive','scalar','integer'});

  % Begin recursion
  % Base cases
  if d==0
    vd = zeros(size(v));

  elseif d==1
    vd = v;

  elseif d==2
    vd = kron(v,v);

  % Recursive step
  else
    if (mod(d,2)==1) % the odd case
      tmp = KroneckerPower(v,(d-1)/2);
      vd  = kron(kron(v,tmp),tmp);

    else % d is even
      tmp = KroneckerPower(v,d/2);
      vd  = kron(tmp,tmp);
    end
  end

end % function KroneckerPower

