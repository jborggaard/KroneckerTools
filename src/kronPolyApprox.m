function [A,B,Nxx,Nxu] = kronPolyApprox(f,g,n,m,degree,x0)
%kronPolyApprox Uses symbolic tools to find approximate polynomial system
%   Given a non-linear system of the form
%
%        \dot{x} = f(x) + g(x)u,    f:R^n->R^n   g:R^n-> R^(n x m)
%
%  this function uses symbolic tools to find a polynomial approximation to this
%  system in Kronecker form, namely
%
%        \dot{x} = Ax + Bu + Nxx{2}kron(x,x) + ... + Nxx{d}kron(x,kron(x,...,x)
%                  + Nxu{1}kron(x,u) + ... Nxu{d}kron(x,kron(x,...kron(x,u))
%
%  Note:  This is an easy implementation, but not efficient, esp. for large d.
%
%  Usage:
%         [A,B,Nxx,Nxu] = kronPolyApprox(f,g,n,m,d)
%
%  Variables:
%         f - a function handle
%         g - a function handle
%         n - number of states
%         m - number of control inputs
%         d - a small integer
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository.
%%

  if (nargin<6)      % set optional argument if not expanding about nonzero x0.
    x0 = zeros(n,1);
  end

  if (m>2 && degree>2)
    fprintf('multiple control inputs haven''t been implemented yet\n')
    fprintf('Nxu will be returned as empty\n')
  end

  syms('x',[n,1])

  J = jacobian(f(x),x);
  A = double(subs(J,x,x0));

  B = double(subs(g(x),x,x0));

  Nxx = cell(1,degree);
  for d=2:degree
    %Jnext = J;
    for j=1:n
      for i=1:n^(d-1)
        Jnext(:,i+n^(d-1)*(j-1)) = diff(J(:,i),x(j));
      end
    end
    Jnext = Jnext/d;
    Nxx{d} = double(subs(Jnext,x,x0));
    J = Jnext;
  end


  if (m==1)
    Nxu = cell(1,degree-1);
    J = jacobian(g(x),x);  Jnext = J;
    Nxu{1} = double(subs(J,x,x0));
    for d=2:degree-1
      for j=1:n
        for i=1:n^(d-1)
          Jnext(:,i+n^(d-1)*(j-1)) = diff(J(:,i),x(j));
        end
      end
      Jnext = Jnext/d;
      Nxu{d} = double(subs(Jnext,x,x0));
      J = Jnext;
    end

  elseif (degree==2)
    Nxu = cell(1,degree-1);
    G = g(x);
    for k=1:m
      J = jacobian(G(:,k));
      Nxu{1}(:,k:m:n*m) = double(subs(J,x,x0));
    end

  else
    Nxu = [];
    
  end

end