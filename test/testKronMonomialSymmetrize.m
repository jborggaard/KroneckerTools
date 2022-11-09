function [] = testKronMonomialSymmetrize()
%  This function runs a number of tests of kronMonomialSymmetrize
%
%  The comparison is made with the composition of kron2CT and CT2kron
%  which is now found in kronPolySymmetrize_old.
%
  addpath('../src')
  addpath('../util')
  
  
  n = 4;
  d = 3;
  c = rand(1,n^d);
  
  symError = tester(c,n,d);
  
  assert( symError < 1e-15*sqrt(n^d) )
  
  n = 5;
  d = 7;
  c = rand(1,n^d);
  
  symError = tester(c,n,d);
  
  assert( symError < 1e-15*sqrt(n^d) )
end

function val = tester(c,n,d) 
  fprintf('Evaluating a testcase with n=%d, d=%d:\n',n,d)
  tic, cSym = kronMonomialSymmetrize(c,n,d); CPUtime = toc;
  fprintf('  CPU time for symmetrizing using gather-scatter:  %g\n',CPUtime)

  tic, cOld = kronPolySymmetrize_old(c,n,d); CPUtime = toc;
  fprintf('  CPU time for symmetrizing using kron2CT CT2kron: %g\n',CPUtime)

  val = norm(cSym-cOld);
  fprintf('  The error is %g\n\n',val)
end