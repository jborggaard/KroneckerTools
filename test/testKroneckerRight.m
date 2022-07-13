%  This script tests the kroneckerRight function using random matrices
%
%  The first tests are for multiplications of the form
%     w * ( rand(n,m1) \otimes rand(n,m2) \otimes rand(n,m3) )
%  and
%     B * ( rand(n,m1) \otimes rand(n,m2) \otimes rand(n,m3) )
%
%  The second tests are for the constant matrix case (M is n \times m1)
%     w * ( M \otimes M \otimes M )
%  and
%     B * ( M \otimes M \otimes M )
%
%  where w = rand(1,n^3) and B = rand(nB,n^3)

addpath('./src')

% Generate test matrices
m1 = 4; m2 = 3; m3 = 5; n = 5; nB = 7;
clear M
M{1} = rand(n,m1);
M{2} = rand(n,m2);
M{3} = rand(n,m3);
w    = rand(1,n^3);
B    = rand(nB,n^3);

% Perform test 1
[wM] = kroneckerRight(w(1:n^2),M(1:2));  % test the base case
absErr = norm(wM - w(1:n^2)*kron(M{1},M{2}),inf);
fprintf('testKroneckerRight: in test 1a, ||absErr|| = %g\n',absErr);

[wM] = kroneckerRight(w,M);
absErr = norm(wM - w*kron(kron(M{1},M{2}),M{3}),inf);
fprintf('testKroneckerRight: in test 1b, ||absErr|| = %g\n',absErr);

[BM] = kroneckerRight(B,M);
absErr = norm(BM - B*kron(kron(M{1},M{2}),M{3}),inf);
fprintf('testKroneckerRight: in test 1c, ||absErr|| = %g\n\n',absErr);

% Perform test 2
M2 = rand(n,m1);
[wM] = kroneckerRight(w,M2);
absErr = norm(wM - w*kron(kron(M2,M2),M2),inf);
fprintf('testKroneckerRight: in test 2a, ||absErr|| = %g\n',absErr);

M2 = rand(n,n);  % test the square matrix case
[wM] = kroneckerRight(w,M2);
absErr = norm(wM - w*kron(kron(M2,M2),M2),inf);
fprintf('testKroneckerRight: in test 2b, ||absErr|| = %g\n',absErr);

M2 = rand(n,m1);
B  = rand(nB,n^3);
[BM] = kroneckerRight(B,M2);
absErr = norm(BM - B*kron(kron(M2,M2),M2),inf);
fprintf('testKroneckerRight: in test 2c, ||absErr|| = %g\n',absErr);
