function [c] = kronMonomialSymmetrize(c, b, d, verbose)
%kronMonomialSymmetrize Symmetrization of Kronecker coefficients.
%
%  Input Variables:
%     c - coefficients of the Kronecker polynomial (row vector of coefficients)
%         c has size (1,b^degree) or (b^degree,1), where
%     b - number of variables in the multinomial
%     degree - degree of the multinomial terms (c should have n^degree terms)
%
%  Returns: symmetrized coefficients of the multinomial (same shape as c)
%
%  Authors: John Burkardt, University of Pittsburgh
%           Jeff Borggaard, Virginia Tech
%           Gemini Code Assist, used Dec 2025
%
%  Licence: MIT
%
%  Documented at: https://people.sc.fsu.edu/~jburkardt/m_src/vector/vector.html
%
%  Part of the KroneckerTools repository.
%%

  % preprocess the inputs
  if ( nargin < 4 ) 
    verbose = false; 
  end
  if ( nargin < 3 ) 
    error('kronMonomialSymmetrize expects 3 consistent inputs: coefficients, number of variables, and degree')
  end

  % test the inputs for consistency
  n = b^d;
  
  if numel(c) ~= n
      error('Input vector c has length %d, but expected b^d = %d', numel(c), n);
  end

  % Ensure c is a column vector for processing
  isRow = isrow(c);
  c = c(:);

  if verbose
      fprintf('Processing vector of size %d (b=%d, d=%d)...\n', n, b, d);
  end

  % Generate indices using the smallest data type required to store the
  % matrix of size (b^d x d) where each row is an index tuple.
  
  if b < 256
      intType = 'uint8';
  elseif b < 65536
      intType = 'uint16';
  else
      intType = 'uint32';
  end
  
  % Preallocate indices. 
  idx_mat = zeros(n, d, intType);
  
  % Create a temporary vector 0 to n-1
  tmp = 0:n-1;
  
  % Fill columns (Convert linear index to base-b subscripts)
  for k = d:-1:1
    idx_mat(:, k) = mod(tmp, b) + 1;
    tmp = floor(tmp / b);
  end

  % Identify similar monomials by sorting rows, e.g. [2, 1, 3] becomes [1, 2, 3]
  % Permutations of the same term will now have identical rows.
  idx_mat = sort(idx_mat, 2);

  % Map unique rows to group monomial terms
  % groupIDs is a vector size N where equivalent terms share the same ID.
  [~, ~, groupIDs] = unique(idx_mat, 'rows');

  % Average and Scatter (the averages are the balanced coefficients)
  % Sum values per group
  groupSums = accumarray(groupIDs, c);
  % Count elements per group
  groupCounts = accumarray(groupIDs, 1);
  
  % Calculate means
  groupMeans = groupSums ./ groupCounts;
  
  % Map means back to original positions
  c = groupMeans(groupIDs);

  % return the symmetrized coefficients in the same form (row or column)
  if isRow
    c = c.';
  end
  
  if verbose
    fprintf('Symmetrization complete.\n');
  end
end
