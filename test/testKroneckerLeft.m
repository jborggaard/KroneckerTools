%  This script tests the kroneckerLeft function using random matrices
%
%  The first tests are for matrices of the form
%    ( rand(m1,n) \otimes rand(m2,n) \otimes rand(m3,n) ) * w
%  where w is a vector, and
%    ( rand(m1,n) \otimes rand(m2,n) \otimes rand(m3,n) ) * B
%  where B is a matrix.
%
%  The second test is for the constant matrix case (M is m1 \times n)
%    ( M \otimes M \otimes M ) * w
%  where w is a vector, and
%    ( M \otimes M \otimes M ) * B
%
%  For compatibility, w = rand(n^3,1) and B = rand(n^3,nB)

addpath('../src/')

% Generate test matrices
m1 = 4; m2 = 3; m3 = 5; n = 5; nB = 7;
clear M
M{1} = rand(m1,n);
M{2} = rand(m2,n);
M{3} = rand(m3,n);
w    = rand(n^3,1);
B    = rand(n^3,nB);

% Perform test 1
[Mw] = kroneckerLeft(M(1:2),w(1:n^2));
absErr = norm(Mw - kron(M{1},M{2})*w(1:n^2),inf);
fprintf('testKroneckerLeft: in test 1a, ||absErr|| = %g\n',absErr);

[Mw] = kroneckerLeft(M,w);
absErr = norm(Mw - kron(kron(M{1},M{2}),M{3})*w,inf);
fprintf('testKroneckerLeft: in test 1b, ||absErr|| = %g\n',absErr);

[Mw] = kroneckerLeft(M,B);
absErr = norm(Mw - kron(kron(M{1},M{2}),M{3})*B,inf);
fprintf('testKroneckerLeft: in test 1c, ||absErr|| = %g\n\n',absErr);

% Perform test 2
M2 = rand(m1,n);
[Mw] = kroneckerLeft(M2,w);
absErr = norm(Mw - kron(kron(M2,M2),M2)*w,inf);
fprintf('testKroneckerLeft: in test 2a, ||absErr|| = %g\n',absErr);

M2 = rand(n,n);
[Mw] = kroneckerLeft(M2,w);
absErr = norm(Mw - kron(kron(M2,M2),M2)*w,inf);
fprintf('testKroneckerLeft: in test 2b, ||absErr|| = %g\n',absErr);

M2 = rand(m1,n);
B  = rand(n^3,nB);
[Mw] = kroneckerLeft(M2,B);
absErr = norm(Mw - kron(kron(M2,M2),M2)*B,inf);
fprintf('testKroneckerLeft: in test 2c, ||absErr|| = %g\n',absErr);
