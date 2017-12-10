function C1 = EntComp(sqmat,IFIM)
%{
  complexity = EntComp(square matrix, type switch)

  Calculate the C1 or C1f Entropic Complexity Measure of a matrix 
  (presumably covariance matrix) from multivariate data.

  Where
  square matrix --- (pxp) covariance matrix for analysis
  type switch --- Scalar indicating how to do computation:
     1 = compute C1 complexity of matrix; 0 = compute complexity
     of IFIM depending on sigma and Dp by opening it up; -1 = compute C1F
     frobenius norm version complexity of matrix
  complexity --- Scalar complexity of the provided matrix

  Example: c1 = EntComp(randn(4),1)

  Copyright (C) 2006 J. Andrew Howe
%}

[p,p2] = size(sqmat);

if (p ~= p2) || nargin ~= 2
    % not square matrix, switch not passed
    fprintf('EntComp: INVALID USAGE-Please read the following instructions!\n'), help EntComp, return
end

cmdet = det(sqmat); cmtra = trace(sqmat); cmdim = rank(sqmat);

if (cmdet == 0) && (IFIM ~= -1) % this only matters if not C1F version
    disp('WARNING: Matrix is singular - determinant is 0.')
    C1 = Inf;
else
    if (IFIM == 1) || (isscalar(sqmat)) % force a 0 complexity for a scalar
        C1 = cmdim*0.5*log(cmtra/cmdim) - 0.5*log(cmdet);
    elseif IFIM == 0                    % covariance matrix passed, do the opened up version
        p1 = (2*p + p*(p+1))/4;
        p2 = cmtra + 0.5*trace(sqmat^2) + 0.5*cmtra^2 + sum(diag(sqmat).^2);
        p3 = (2*p + p*(p+1))/2;
        p4 = log(cmdet)*(p+2)/2 + p*log(2)/2;
        C1 = p1*log(p2/p3)-p4;
    elseif IFIM == -1                   % Frobenius Norm version
        ei = eig(sqmat); C1 = mean(ei); C1 = sum((ei-C1).^2)/(4*C1^2);
    end
end

%{
JAH 20060128, adapted for octave 3.4.3 20120315

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
