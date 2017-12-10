function [Dp,mopeninv] = DupMatrix(p)
%{
  [dupmat, mopeninv] = DupMatrix(number dimensions)

  Returns Magnus and Neudecker's duplication matrix of size p, and the
  moore-penrose inverse of Dp - Dp+.

  Where
  number dimensions  --- scalar integer p
  dupmat --- (p^2 x 0.5*p*(p+1)) duplication matrix
  mopeninv --- moore-penrose inverse of duplication matrix

  [Dp,mopeninv] = DupMatrix(2)

  Copyright (C) 2006 J. Andrew Howe
%}

a = tril(ones(p));
i = find(a);
a(i) = 1:length(i);
a = a + tril(a,-1)';

a = a(:);
m = p*(p+1)/2;
Dp = zeros(p*p,m);
for r = 1:size(Dp,1)
    Dp(r, a(r)) = 1;
end

if nargout == 2; mopeninv = inv(Dp'*Dp)*Dp'; end;

%{
JAH 20061106, adapted for octave 3.4.3 20120305

Copyright (C) 2007 J. Andrew Howe

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