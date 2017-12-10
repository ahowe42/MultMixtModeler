function [bin_comb, subs] = VarSubset(datamatrix)
%{
  [binary results, subsets] = VarSubset(data matrix)

  Compute all combinatorial nonempty subsets of columns of a data matrix.

  Where
  data matrix --- (nxp) matrix of data with p variables. You can also pass
    in just p, the number of dimensions.  If you do this, the second
    output will not be available.
  binary results --- (2^p-1)x(p+1) matrix; first column indicates number
    of variables in that subset, the rest of the row is a binary string
    that can be used to subset into the original data matrix.
  subsets --- cell array with all 2^p - 1 subsets (ignoring empty set)

  Example: dat = rand(10,4); [bin,subs] = VarSubset(dat); bin, subs{8}

  See Also ICSubTable

  Copyright (C) 2006 J. Andrew Howe; see below
%}

if (nargin ~= 1)
    % wrong number arguments
    fprintf('VarSubset: INVALID USAGE-Please read the following instructions!\n'), help VarSubset, return
end

if isscalar(datamatrix)
    p = datamatrix;
else
    p = size(datamatrix,2);
end

bin_comb = dec2bin([1:2^p]);    % make a binary string with all combinations of vars
bin_comb = bin_comb([1:(end-1)],[2:end]);   % drop the leftmost column and the last row
bin_comb = (bin_comb == '1');   % logical(char) doesn't work, but this does
bin_comb = sortrows([sum(bin_comb,2), bin_comb],-1);    % group by # of vars (this will unlogicalize the binaries

if isscalar(datamatrix); return; end;   % user passed in just p, so no data to subset

subs = cell(2^p-1,1);
for subcnt = 1:(2^p-1) % subtract 1, since we ignore the option of no variables
    subs{subcnt} = datamatrix(:,logical(bin_comb(subcnt,[2:end])));
end

%{
JAH 20060131, adapted for octave 3.4.3 20120305

Copyright (C) 2006 J. Andrew Howe

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}