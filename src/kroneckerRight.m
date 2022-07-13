function [Y] = kroneckerRight(B,M)
%kroneckerRight multiplies a matrix by a special Kronecker product matrix.
%         Y = B*kron(M{1},kron(M{2},kron(M{3},...))
%
%  where B has m rows and each matrix M{k} has size(n1(k),n2(k)).
%
%  Usage:  Y = kroneckerRight(B,M)
%
%  Variables:   B    a matrix with size (m,prod(n1))
%               M    a cell array of matrices where M{k} has size (n1(k),n2(k))
%
%  This is multiplication performed recursively using Kronecker product rules.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%
  [m, nB] = size(B);

  if (iscell(M))
    d = length(M);
    n1 = zeros(1,d); n2 = zeros(1,d);
    for i=1:d
      [n1(i),n2(i)] = size(M{i});
    end

    if ( d==2 )
      Y = zeros(m,n2(1)*n2(2));
      for i=1:m
        Y(i,:) = reshape( M{2}.'*reshape(B(i,:),n1(2),n1(1))*M{1}, 1,n2(1)*n2(2) );
      end

    else
      n2Y = prod(n2);
      Y = zeros(m,n2Y);
      for i=1:m
        T = M{d}.'*reshape(B(i,:),n1(d),nB/n1(d));

        Z = zeros(n2(d),prod(n2(1:end-1)));
        for j=1:n2(d)
          Z(j,:) = kroneckerRight(T(j,:),M(1:end-1));
        end
        Y(i,:) = Z(:).';
      end
    end

  else % we consider the case where M{i} = M for each i:
    Y = kroneckerLeft(M.',B.').';
  
  end
end

