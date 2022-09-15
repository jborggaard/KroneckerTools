%  Script to test the LyapProduct function

addpath('../src')

n = 4;
M = rand(n,n);

% two term test
d = 2;
v = rand(n^d,1);
Mv = LyapProduct(M,v,d);
fullMv = (kron(M,eye(n))+kron(eye(n),M))*v;
e2 = norm(Mv - fullMv);
assert(e2<1e-14) % should scale tolerance with n^d

% three term test
d = 3;
v = rand(n^d,1);
Mv = LyapProduct(M,v,d);
fullMv = (kron(M,eye(n^2))+kron(eye(n),kron(M,eye(n))) + kron(eye(n^2),M))*v;
e3 = norm(Mv - fullMv);
assert(e3<1e-14) % should scale tolerance with n^d

% four term test
d = 4;
v = rand(n^d,1);
Mv = LyapProduct(M,v,d);
fullMv = (kron(M,eye(n^3))+kron(eye(n),kron(M,eye(n^2))) + ...
          kron(eye(n^2),kron(M,eye(n))) + kron(eye(n^3),M))*v;
e4 = norm(Mv - fullMv);
assert(e4<1e-13) % should scale tolerance with n^d

% five term test
d = 5;
v = rand(n^d,1);
Mv = LyapProduct(M,v,d);
fullMv = (kron(M,eye(n^4)) + ...
          kron(eye(n),kron(M,eye(n^3))) + ...
          kron(eye(n^2),kron(M,eye(n^2))) + ...
          kron(eye(n^3),kron(M,eye(n))) + ...
          kron(eye(n^4),M))*v;
e5 = norm(Mv - fullMv);
assert(e5<1e-13) % should scale tolerance with n^d
