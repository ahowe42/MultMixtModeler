function probcode = MatrixProblem(sqmat,printmess,rclim)
%{
  problem code = MatrixProblem(square matrix, print message, rcond limit)

  This function will take in a square, presumably covariance, matrix and
  evaluate it for ill-conditionedness and positive definitedness, in that
  order.

  Where
  square matrix --- Square matrix to be evaluated.
  print message --- Optional logical, if 1, this will print a message.
  rcond limit --- Optional, smallest reciprocal condition number to still
   be considered ok; default is 1e-10.
  problem code --- 0 = no problem; 1 = ill-conditioned;  2 = not positive definite

  Example: MatrixProblem(rand(4),1)

  
  Copyright (C) 2007 J. Andrew Howe; see below
%}

[d1,d2] = size(sqmat);

if (d1 ~= d2) || (sum(nargin == [1,2,3]) ~= 1)
    % not a square matrix or wrong number arguments
    fprintf('MatrixProblem: INVALID USAGE-Please read the following instructions!\n'), help MatrixProblem, return
end

probcode = 0;
if nargin == 1; rclim = 1e-10; printmess = 0; end;
if nargin == 2; rclim = 1e-10; end;

if rcond(sqmat) < rclim
    if printmess
        disp('Matrix Is Ill-Conditioned')
    end
    probcode = 1; return;   % matrix is ill-conditioned
end

[R,p] = chol(sqmat);
if p > 0
    if printmess
        disp('Matrix Is Not Positive Definite')
    end
    probcode = 2; return;   % matrix is not positive definite
end

if printmess
    disp('Matrix Has No Detected Problems')
end

%{
JAH 20070305, adapted for octave 3.4.3 20120305

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