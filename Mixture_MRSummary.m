%{
  Usage: Mixture_MRSummary
  After performing many simulations/replications of M3, this script will
  provide summary results across all runs.  The workspace .mat file from the
  multirun must already be loaded in the workspace.

  Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe
%}

clc; 
if (exist([mydir,GAmat([1:(end-4)]),'_SMRY.out'],'file') == 2)
    delete([mydir,GAmat([1:(end-4)]),'_SMRY.out']);   % delete existing summary file
end
diary([mydir,GAmat([1:(end-4)]),'_SMRY.out'])
disp(repmat('#',1,50)), disp('Parameters'), Mixture_DispParms, disp(repmat('#',1,50))
% summarize by information criteria
sumsum = ones(length(MMM.InfCrit),4)*NaN; disp(' '), usd = [];
for infcnt = 1:length(MMM.InfCrit)    
    if isempty(MMM.InfCrit{infcnt}); continue; end;  % skip this one
    disp([repmat('+-',1,24),'+'])
    % compute summary
    thisic = (GAscores(:,2) == infcnt);
    tab1 = table2str({'Replic.';'k';MMM.InfCrit{infcnt};'%_Correct';'#_Clusters'},...
        [GAscores(thisic,[1,3,4,5,6])],{'%0.0f';'%0.0f';'%0.2f';'%0.3f';'%0.0f'},1);
    [val1,ind1] = min(GAscores(thisic,4)); val1 = GAscores(thisic,[4,5]);
    val1 = val1(ind1,:); tab1(ind1 + 3,1) = '*'; bstrun = GAfil(ind1,:);
    % compute tabulation
    tab2 = tabulate(GAscores(thisic,6)); [val2,ind2] = max(tab2(:,2)); val2 = tab2(ind2,[1,3]);
    tab2 = table2str({'#_Clusters';'%_Selected'},tab2(:,[1,3]),{'%0.0f','%0.2f'},1);
    tab2(ind2 + 3,1) = '*';
    % compute summary of summary
    sumsum(infcnt,:) = [val2,val1];
    % display
    disp(tab1), disp(['* best replication: ',bstrun])
    disp(repmat('-',1,50))
    disp(tab2), disp('* most frequently selected k')
    disp([repmat('+-',1,24),'+']), disp(' ');
    usd = [usd,infcnt];
end                 % information criteria loop
disp('SUMMARY: Multivariate Mixture Modeler')
disp(GAmat([1:(end-4)]))
disp([repmat('+-',1,24),'+'])
disp(table2str({'#_Clusters';'%_Selected';'Score';'%_Correct'},sumsum(usd,:),{'%0.0f';'%0.2f';'%0.2f';'%0.2f'},0,MMM.InfCrit(usd)))
tottim = etime(clock,MMM.totstt);
disp(sprintf('Total Elapsed Time:\n\t%1.0f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',tottim./[1,60,3600]));
disp([repmat('+-',1,24),'+']), diary off

%{
JAH 20071015, adapted for octave 3.4.3 20120319

Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe

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