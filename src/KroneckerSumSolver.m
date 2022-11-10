function [x] = KroneckerSumSolver(A,b,degree,M)
%KroneckerSumSolver Provides an efficient solution for a Kronecker sum system
%
%  Implements an N-Way version of the Bartels-Stewart algorithm for a 
%  special Kronecker sum system (the special case all matrices are the same)
%  with an additional
%
%  (kron(A{d},eye(n^(d-1))) + ... + kron(eye(n^(d-1),A{1}))+diag(M))x = b
%
%  where A is a cell array containing d matrices each of size (n,n) 
%  and b and x are size (n^d,1).  M is an nxn matrix that is added to
%  block diagonal:  diag(M) = kron(eye(n^(d-1)),M)
%
%  Usage:  [x] = KroneckerSumSolver(A,b,d,M)
%
%  For now, we assume each A{i} is the same (an N-Way Lyapunov equation) but
%  this can be relaxed by performing a Schur factorization to each A{i}
%  and performing the proper change of variables.
%
%  This is has been developed for the NLbalancing repository and its 
%  functionality is consistent with the function KroneckerSumSolver that 
%  was developed for the QQR problem (without the matrix M).
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%

  n = size(A{1},1);

  
  %  As a simplification, we assume all A{i} are the same size.  We furthermore
  %  assume they are all the same (so we can use the kroneckerLeft function as is)
  %  Both of these can be easily relaxed if an application arises.
  [U,T] = schur(A{1},'complex');

  b = kroneckerLeft(U',b);

  L = length(A);
  for l=1:L
    A{l} = T;
  end
  
  if ( nargin<4 )
    UMU = sparse(n,n);
  else
    UMU = U'*M*U;
  end

  X = zeros(n,n^(degree-1));
  B = reshape(b,n,n^(degree-1));

  jIndex = n*ones(1,degree);  jIndex(degree) = n+1;
  jRange = cell(1,degree);

  last = false;
  while( ~last )
    [jIndex,last] = decreaseByOne(jIndex,n);

    diagA = 0;
    for i=2:degree
      diagA = diagA + A{i}(jIndex(i),jIndex(i));
    end
    At = A{1} + diagA*eye(n) + UMU;

    colIdx = jIndex(2);
    for i=3:degree
      colIdx = colIdx + (jIndex(i)-1)*n^(i-2);
    end

    rhs = B(:,colIdx);

    %  Backsubstitution steps
    for i=2:degree    % this could be done by decreaseByOne as well
      jRange{i} = (jIndex(i)+1):n;
    end

    for i=2:degree
      if (~isempty(jRange{i}))
        shift = 1;
        for j=2:(i-1)
          shift = shift + (jIndex(j)-1)*n^(j-2);
        end
        for j=(i+1):degree
          shift = shift + (jIndex(j)-1)*n^(j-2);
        end
        jIdx = shift + (jRange{i}-1)*n^(i-2);

        rhs = rhs - X(:,jIdx)*A{i}(jIndex(i),jRange{i}).';
      end
    end

    X(:,colIdx) = At\rhs;

  end

  x = X(:);
  x = real(kroneckerLeft(U,x));
end


function [jIndex,last] = decreaseByOne(jIndex,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  dp1 = length(jIndex);
  d   = dp1-1;

  if (jIndex(dp1)==1)
    jIndex(dp1) = n;
    jIndex(1:d) = decreaseByOne(jIndex(1:d),n);

  else
    jIndex(dp1) = jIndex(dp1)-1;

  end

  last = norm(jIndex(2:end)-ones(1,d))==0;

end
