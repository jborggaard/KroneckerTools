function [] = testKronMonomialSymmetrize()
%  This function runs a number of tests of kronMonomialSymmetrize
%
%  The comparison is made with the file kronMonomialSymmetrize which
%  existed pre-2026.
%
%  The comparison made with the composition of kron2CT and CT2kron
%  found in kronPolySymmetrize_old is now commented out.
%
  addpath('../src')
  addpath('../util')
  
  
  n = 4;
  d = 3;
  c = rand(1,n^d);
  
  symError = tester(c,n,d);
  
  assert( symError < 1e-15*sqrt(n^d) )
  
  n = 9;
  d = 8;
  c = rand(1,n^d);
  
  symError = tester(c,n,d);
  
  assert( symError < 1e-15*sqrt(n^d) )
end

function val = tester(c,n,d) 
  fprintf('Evaluating a testcase with n=%d, d=%d:\n',n,d)
  tic, cSym = kronMonomialSymmetrize(c,n,d); CPUtime = toc;
  fprintf('  CPU time for symmetrizing using optimized version:  %g\n',CPUtime)

  % tic, cOld = kronPolySymmetrize_old(c,n,d); CPUtime = toc;
  % fprintf('  CPU time for symmetrizing using kron2CT CT2kron: %g\n',CPUtime)
  tic, cOld = kronMonomialSymmetrize_old(c,n,d); CPUtime = toc;
  fprintf('  CPU time for symmetrizing using gather-scatter: %g\n',CPUtime)

  val = norm(cSym-cOld);
  fprintf('  The error is %g\n\n',val)
end