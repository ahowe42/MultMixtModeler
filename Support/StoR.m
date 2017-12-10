function R = StoR(S)
%{
  correlation matrix = StoR(covariance matrix)

  Convert a covariance matrix into correlation form.

  Where
  covariance matrix  --- (pxp) estimated covariance matrix
  correlation matrix --- (pxp) estimated correlation matrix

  Example: R = StoR([1,-2,0;-2,5,0;0,0,2])
  
  Copyright (C) 2007 J. Andrew Howe; see below
%}
  
[p,p2] = size(S);

if (nargin ~=1) || (p ~= p2)
    % wrong number arguments, input not square
    fprintf('StoR: INVALID USAGE-Please read the following instructions!\n'), help StoR, return
end

dg = diag(S);               % just get the diagonal elements
isdg = diag(1./sqrt(dg));   % reciprocal of standard deviations
R = isdg*S*isdg;

%{
JAH 20070812, adapted for octave 3.4.3 20120315

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