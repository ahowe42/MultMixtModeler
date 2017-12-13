function total_var = KMix_TotInClustEUDist(data, labels)
%{
  total distance = KMix_TotInClustEUDist(data matrix, mixture assignments)
  This will compute the total within cluster Euclidian variance for a mixture model.

  Where
  data matrix --- (nxp) data matrix
  mixture assignments --- n-valued vector of assignments 1:k
  total distance --- total within cluster variance

  Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}
  
[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 2)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMix_TotInClustEUDist: INVALID USAGE-Please read the following instructions!\n'), help KMix_TotInClustEUDist, return
end

ks = unique(labels); kcnt = length(ks); % must do this in case some k missing (shouldn't)
total_var = zeros(1,kcnt); means = zeros(1,p,kcnt);

% compute the cluster means
for mixcnt = 1:kcnt
    means(:,:,mixcnt) = mean(data(labels == ks(mixcnt),:),1);
end     % mixtures loop

% get the total variance for each cluster
for datcnt = 1:n
    thisclst = find(labels(datcnt) == ks);
    total_var(thisclst) = total_var(thisclst) + sum((data(datcnt,:) - means(:,:,thisclst)).^2);
end     % data loop
total_var = sum(total_var);

%{
JAH 20061215, adapted for octave 3.4.3 20120324

Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe

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