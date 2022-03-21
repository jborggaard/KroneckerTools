function [Y] = kroneckerLeft(M,B)
%kroneckerLeft multiplies a Kronecker product of matrices with a matrix.
%         Y = kron(M{1},kron(M{2},kron(M{3},...kron(M{d-1},M{d})))*B
%
%  Usage:  Y = kroneckerLeft(M,B)
%
%  Variables:   M    a cell array of matrices where M{k} has size (n1(k),n2(k))
%               B    a matrix with size (prod(n2),m)
%
%  This is multiplication performed recursively using Kronecker product rules.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Part of the KroneckerTools repository: github.com/jborggaard/KroneckerTools
%%
  d = length(M);
  n1 = zeros(1,d); n2 = zeros(1,d);
  for i=1:d
    [n1(i),n2(i)] = size(M{i});
  end

  [nB ,m] = size(B);

  if ( d==2 )
    Y = zeros(n1(1)*n1(2),m);
    for j=1:m
      Y(:,j) = reshape(M{2}*reshape(B(:,j),n2(2),n2(1))*M{1}.',n1(1)*n1(2),1);
    end
    
  else
    n1Y = prod(n1);
    Y = zeros(n1Y,m);
    for j=1:m
      T = reshape(B(:,j),nB/n2(1),n2(1))*M{1}.';
      
      Z = zeros(prod(n1(2:end)),n1(1));
      for k=1:n1(1)
        Z(:,k) = kroneckerLeft(M(2:end),T(:,k));
      end
      Y(:,j) = Z(:);
    end
  end
  
end

