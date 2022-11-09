function [c] = kronMonomialSymmetrize(c,b,d,verbose)
%kronMonomialSymmetrize() symmetrizes Kronecker monomial coefficients.
%
%  A gather/scatter operation is used to redistribute the coefficients
%  of a Kronecker monomial term so that each coefficient is the same.
%  Usage:
%     c = kronPolySymmetrize(c,n,degree);
%
%  Input Variables:
%     c - coefficients of the Kronecker polynomial (vector of coefficients)
%     n - number of variables in the multinomial
%     degree - degree of the multinomial terms (c should have n^degree terms)
%
%  Returns: symmetrized coefficients of the multinomial
%
%  Author: John Burkardt, University of Pittsburgh
%          (01 November 2022)
%
%  Licence: MIT
%
%  Part of the KroneckerTools repository.
%%

%  Process the inputs
  if ( nargin < 4 )
    verbose = false;
  end

  if ( nargin < 3 )
    d = 3;
  end

  if ( nargin < 2 )
    b = 3;
  end

  n = b^d;

  if ( nargin < 1 ) %  Create data vector C.
    c = demo_data ( d, b );
  end

  if ( verbose )
    %
    %  Report C.
    %
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Original data C:\n' );
    a = [];
    more = false;
    i = 0;
    
    while ( true )
  
      [ a, more ] = vector_next ( d, b, a, more );
  
      if ( ~ more )
        break
      end
  
      i = i + 1;
      fprintf ( 1, '  #%2d  C=%g, [ %d, %d, %d ]\n', i, c(i), a );
  
    end
  end

  %
  %  Use sorted copy Z of vector A as a key.
  %
  c2 = zeros ( 1, n );
  z = [];
  more = false;
  j = 0;

  if ( verbose )
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Data gathered to representatives and averaged:\n' );
    fprintf ( 1, '\n' );
  end

  while ( true )
%
%  Compute the next Z value.
%
    [ z, more ] = vector_next_representative ( d, b, z, more );

    if ( ~ more )
      break
    end
%
%  Generate all vectors A whose sorted value is Z.
%  Gather those C values into one vector C2.
%
    j = j + 1;
    mult = 0;
    a = z;
    more2 = false;
    while ( true )
      mult = mult + 1;
      i = vector_rank ( d, b, a );
      c2(j) = c2(j) + c(i);
      [ a, more2 ] = vector_next_equivalent ( a );
      if ( ~ more2 )
        break;
      end
    end
%
%  Average C2 by number of contributors.
%
    c2(j) = c2(j) / mult;

    if ( verbose )
      fprintf ( 1, '  #%2d  mult=%d  C2=%g, [ %d, %d, %d ]\n', j, mult, c2(j), z );
    end
%
%  Scatter averaged value back to contributors.
%
    a = z;
    more2 = false;
    while ( true )
      i = vector_rank ( d, b, a );
      c(i) = c2(j);
      [ a, more2 ] = vector_next_equivalent ( a );
      if ( ~ more2 )
        break;
      end
    end

  end

  if ( verbose )
    %
    %  Report Symmetrized C.
    %
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Symmetric C:\n' );
    a = [];
    more = false;
    i = 0;
    
    while ( true )
  
      [ a, more ] = vector_next ( d, b, a, more );
  
      if ( ~ more )
        break
      end
  
      i = i + 1;
      fprintf ( 1, '  #%2d  C=%g, [ %d, %d, %d ]\n', i, c(i), a );
  
    end

  end

  return
end
function c = demo_data ( d, b )

%*****************************************************************************80
%
%% demo_data() defines a data vector for the gather/scatter demo.
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    01 November 2022
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer d: the vector dimension.
%
%    integer b: vector entries range from 1 to b.
%
%  Output:
%
%    integer c(b^d): a data value for every unique d-dimensional vector
%    with entries from 1 to b.  Lexicographic order is assumed.
%
  n = b^d;
  c = zeros ( 1, n );
  f = linspace ( 1, d, d );
%
%  Generate each vector A, and compute the corresponding C.
%
  a = [];
  more = false;
  i = 0;
  
  while ( true )

    [ a, more ] = vector_next ( d, b, a, more );

    if ( ~ more )
      break
    end

    i = i + 1;
    c(i) = f * a';

  end

  return
end

function [ a, more ] = vector_next ( d, b, a, more )
%*****************************************************************************80
%
%% vector_next() generates vectors in lexicographic order.
%
%  Discussion:
%
%    The vectors are produced in lexical order, starting with
%    (1,1,...,1),
%    (1,1,...,2),
%    ...
%    (B,B,...,B).
%
%  Example:
%
%    d = 2,
%    b = 3
%
%    1   1
%    1   2
%    1   3
%    2   1
%    2   2
%    2   3
%    3   1
%    3   2
%    3   3
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    31 October 2022
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer d: the vector dimension.
%
%    integer b: vector entries range from 1 to b.
%
%    integer A(D): except on the first call, this should
%    be the output value of A on the last call.
%
%    logical MORE: should be FALSE on the first call, and
%    thereafter should be the output value of MORE from the previous call.  
%
%  Output:
%
%    integer A(D): the next vector.
%
%    logical MORE: is TRUE if another vector was computed.
%    If MORE is FALSE on return, then ignore the output value A, and
%    stop calling the routine.
%
  if ( ~ more )

    a(1:d) = 1;
    more = true;

  else
      
    for i = d : -1 : 1

      a(i) = a(i) + 1;

      if ( a(i) <= b )
        return
      end

      a(i) = 1;

    end

    more = false;

  end

  return
end

function [ a, more ] = vector_next_equivalent ( a )
%*****************************************************************************80
%
%% vector_next_equivalent() generates equivalent vectors.
%
%  Discussion:
%
%    Two vectors are equivalent if one is simply a permuted copy of the 
%    other.  A vector may have no equivalents, or many.
%    This function produces the "next" equivalent vector, so it can
%    produce all the equivalent vectors one at a time.  To do so,
%    start with the vector whose entries are sorted in ascending order.
%
%    For example, here is for vectors of length 3, with digits 1, 2 or 3,
%    the possible sequences are:
%
%      #1  1 1 1 --> 1 1 1
%      #2  1 1 2 --> 1 1 2, 1 2 1, 2 1 1
%      #3  1 1 3 --> 1 1 3, 1 3 1, 3 1 1
%      #4  1 2 2 --> 1 2 2, 2 1 2, 2 2 1
%      #5  1 2 3 --> 1 2 3, 1 3 2, 2 1 3, 2 3 1, 3 1 2, 3 2 1
%      #6  1 3 3 --> 1 3 3, 3 1 3, 3 3 1
%      #7  2 2 2 --> 2 2 2
%      #8  2 2 3 --> 2 2 3, 2 3 2, 3 2 2
%      #9  2 3 3 --> 2 3 3, 3 2 3, 3 3 2
%     #10  3 3 3 --> 3 3 3
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    02 November 2022
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer a(d): a vector in the sequence.
%
%  Output:
%
%    integer a(d): the next equivalent vector.  But if the sequence is
%    exhausted, a is returned as [].
%
%    logical MORE: is TRUE if another vector was computed.
%    If MORE is FALSE on return, then ignore the output value A, and
%    stop calling the routine.
%
  more = true;

  n = length ( a );
%
%  Starting at the right, seek the highest index I for which A(I) < A(I+1).
%
  i = n - 1;

  while ( true )

    if ( i <= 0 )
      break
    end

    if ( a(i) < a(i+1) )
      break
    end

    i = i - 1;

  end
%
%  If no I could be found, then we have reach the final permutation,
%  N, N-1, ..., 2, 1.  Time to start over again.
%
  if ( i == 0 )
    a = [];
    more = false;
  else
%
%  Otherwise, let J be the greatest index after I such that A(I) < A(J).
%
    j = n;
    while ( a(j) <= a(i) ) 
      j = j - 1;
    end
%
%  Interchange elements I and J.
%
    t    = a(i);
    a(i) = a(j);
    a(j) = t;
%
%  Reverse the elements from I+1 to N.
%
    a(i+1:n) = a(n:-1:i+1);

  end

  return
end
function [ a, more ] = vector_next_representative ( d, b, a, more )

%*****************************************************************************80
%
%% vector_next_representative() generates vector representatives.
%
%  Example:
%
%    d = 3
%    b = 3
%
%    1: 1  1  1
%    2: 1  1  2
%    3: 1  1  3
%    4: 1  2  2
%    5: 1  2  3
%    6: 1  3  3
%    7: 2  2  2
%    8: 2  2  3
%    9: 2  3  3
%   10: 3  3  3
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    31 October 2022
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer d: the vector dimension.
%
%    integer B, the base to be used.  B = 2 will
%    give vectors of 1's and 2's, for instance.
%
%    integer A(D), except on the first call, this should
%    be the output value of A on the last call.
%
%    logical MORE, should be FALSE on the first call, and
%    thereafter should be the output value of MORE from the previous call.  
%
%  Output:
%
%    integer A(D), the next vector.
%
%    logical MORE, is TRUE if another vector was computed.
%    If MORE is FALSE on return, then ignore the output value A, and
%    stop calling the routine.
%
  if ( ~ more )

    a(1:d) = 1;
    more = true;

  else
      
    for i = d : -1 : 1

      a(i) = a(i) + 1;

      if ( a(i) <= b )
        a(i+1:d) = a(i);
        return
      end

    end

    more = false;

  end

  return
end
function rank = vector_rank ( d, b, a )

%*****************************************************************************80
%
%% vector_rank() ranks a vector.
%
%  Example:
%
%    d = 3,
%    b = 3
%    A = ( 2, 1, 3 )
%
%    RANK = 12
%
%  Licensing:
%
%    This code is distributed under the MIT license.
%
%  Modified:
%
%    31 October 2022
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer d: the vector dimension.
%
%    integer b, the upper limit for the array indices.
%
%    integer A(d), the vector to be ranked.
%
%  Output:
%
%    integer RANK, the rank of the index vector, or -1 if A
%    is not a legal index.
%
  rank = -1;

  for i = 1 : d
    if ( a(i) < 1 | b < a(i) )
      return;
    end
  end

  rank = 0;
  for i = 1 : d
    rank = b * rank + a(i);
  end

  rank = 1;
  range = 1;
  for i = 1 : d
    rank = rank + ( a(d+1-i) - 1 ) * range;
    range = range * b;
  end

  return
end

