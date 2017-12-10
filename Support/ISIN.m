function res = ISIN(srch_for, srch_in, ret)
%{
  search results = ISIN(data to find, where to find it, return type)

  This will search a numeric matrix for specific values.

  Where
  data to find --- scalar/vector/matrix of values to find
  where to find it --- vector/matrix to search
  return type --- 1:logical array indicating success or failure; 0: array
    indicating number of occurances
  search results --- array of same size as data to find with results

  Example: d = unidrnd(10,10,1), ISIN([2,5],d,0)

  Copyright (C) 2006 J. Andrew Howe
%}

if nargin ~= 3
    % what was wrong
    fprintf('ISIN: INVALID USAGE-Please read the following instructions!\n'), help ISIN, return
end

srchin = srch_in(:);
res = zeros(size(srch_for));
for srchcnt = 1:prod(size(srch_for))
    res(srchcnt) = sum(srch_for(srchcnt) == srchin);
end
if ret == 1; res = sign(res); end;

%{
JAH 20060222, adapted for octave 3.4.3 20120305

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
