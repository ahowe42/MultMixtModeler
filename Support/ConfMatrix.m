function [cmat,errrats] = ConfMatrix(yactu,ypred,silent)
%{
  [confusion matrix, error rates] = ConfMatrix(actual labels, predicted labels,silent)

  Match actual values to predicted values, and form a matrix of correct 
  and false values and error rates.  In the confusion matrix, rows are for
  the actual labels, and columns are  predicted.  For example, if we have
  0     1     2     3     0
  1    76     6     3    85
  2    13    72     0    85
  3     4     0    81    85
  0    93    78    84   255
  the entry 13 says 13 actual 2's were classified as 1's.

  Where
  actual labels --- (1xn) vector of Ka actual labels
  predicted labels --- (1xn) vector of Kp predicted labels
  silent --- 1 (default) = just compute; 0 = print 
  confusion matrix --- square confusion matrix of size max(Ka,Kp)+2
  error rates --- (max(Ka,Kp) x max(Ka,Kp)) matrix with error rates (relative to n)

  Example: ConfMatrix(sort(unidrnd(5,100,1)),sort(unidrnd(5,100,1)),0);
    
  Copyright (C) 2006 J. Andrew Howe; see below
%}

act_n = length(yactu); act_p = length(ypred);
if (act_n ~= act_p) || (nargin < 2) || (nargin > 3)
    % mismatched lengths, wrong # of args
    fprintf('ConfMatrix: INVALID USAGE-Please read the following instructions!\n'), help ConfMatrix, return
end

if nargin == 2; silent = 1; end;

valus = unique([yactu(:);ypred(:)]);
vcnt = length(valus); cmat = zeros(vcnt);

for cntr1 = 1:vcnt
    ind = (yactu == valus(cntr1));
    for cntr2 = 1:vcnt
        cmat(cntr1,cntr2) = sum(ypred(ind) == valus(cntr2));        
    end
end
errrats = cmat/act_n; errrats(logical(eye(vcnt))) = 0;
cmat = [[cmat,sum(cmat,2)];[sum(cmat,1),act_n]];
cmat = [[0,valus',0];[[valus;0],cmat]];         % tack on labels

if silent == 0
	disp('Confusion Matrix')
	disp(cmat)
	disp('Rows are Actual, Columns are Predicted')
end

%{
JAH 20060413, 20120224

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