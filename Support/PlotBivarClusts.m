function PlotBivarClusts(data, labels, names)
%{
  PlotBivarClusts(data, labels, names)

  Plots bivariate clustered data with a different color and marker for
  each cluster.  It is assumed that the labels are in the {1,2,...,K} range,
  where K is the number of groups.  If names is passed in, the values will
  be used to mark the group centroids, otherwise 'Centroid k' will be used.

  Where
  data  --- bivariate data (nx2)
  labels --- cluster assignments (nx1)
  names --- optional cell array with one entry for each group

  See Also PlotTrivarClusts

  Copyright (C) 2006 J. Andrew Howe; see below
%}

[nd,pd] = size(data); nl = length(labels);

if (nd ~= nl) || ((nargin ~= 2) && (nargin ~= 3)) || (pd ~= 2) || ((nargin  == 3) && not(iscell(names)))
    % mismatching observations, wrong # of args, not bivariate, 3rd opt arg not a cell array
    fprintf('PlotBivarClusts: INVALID USAGE-Please read the following instructions!\n'), help PlotBivarClusts, return
end

uni = unique(labels); k = length(uni);
col = 'bgrmk'; mrk = '.s*o^+<x>vd';
cols = rem([1:k],5); cols(cols == 0) = 5; cols = col(cols);
mrks = rem([1:k],11); mrks(mrks == 0) = 11; mrks = mrk(mrks);

hold on;
for clstcnt = 1:k
    ind = (labels == uni(clstcnt));
    scatter(data(ind,1),data(ind,2),[],cols(clstcnt),mrks(clstcnt)) % JAH 20120214 change for octave
    mn = mean(data(ind,:),1);
    if nargin == 3
        if not(isempty(names{clstcnt}))
            TH = text(mn(1,1),mn(1,2),['*',names{clstcnt}]);
        end
    else
        TH = text(mn(1,1),mn(1,2),['*Centroid ',num2str(clstcnt)]); 
    end
    if exist('TH','var'); set(TH,'FontSize',11); end;
end
axis([min(data(:,1)),max(data(:,1)),min(data(:,2)),max(data(:,2))]);
hold off

%{
JAH 20061012, adapted for octave 3.4.3 20120214

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
