function AIC = KMixMEC_AIC(data,labels,smth,type,ifprob)
%{
  AIC = KMixMEC_AIC(data matrix, mixture assignments, covariance smoother,...
  subtype, problem only)
   This will compute AIC for a mixture of elliptically contoured
   distriutions model.

   Where
   data matrix --- (nxp) data matrix
   mixture assignments --- n-valued vector of assignments 1:k
   covariance smoother --- alpha code to pass covsmooth
   subtype --- 'KT' (Kotz) or 'PVII' (Pearson Type VII)
   problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
   AIC --- value of AIC for MEC mixture model

  Copyright (C) 2008 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}   
   
[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 5)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMixMEC_AIC: INVALID USAGE-Please read the following instructions!\n'), help KMixMEC_AIC, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
mix_meanvecs = zeros(1,p,kcnt); mix_covrmats = zeros(p,p,kcnt);
mix_propors = zeros(1,kcnt);
novars = ones(1,kcnt);
if isequal(type,'KT')
    mix_ECparms = zeros(1,3,kcnt); % 3 parms for Kotz
elseif isequal(type,'PVII')
    mix_ECparms = zeros(1,2,kcnt); % 3 parms for Pearson VII
else
    fprintf('KMixMEC_AIC: TYPE INVALID-Please read the following instructions!\n'), help KMixMEC_AIC, return
end

% get the parameter estimates
for mixcnt = 1:kcnt
    ind = (labels == ks(mixcnt));
    if (sum(ind) <= 1); continue; end;  % empty or singleton cluster - no variance    
    % estimate the parameters
    mix_propors(mixcnt) = sum(ind)/n;                   % cluster proportion
    [mu,sigma,parms] = MECEstimate(data(ind,:),type,smth,ifprob);
    mix_meanvecs(:,:,mixcnt) = mu;                      % cluster means
    mix_covrmats(:,:,mixcnt) = sigma;                   % cluster covariances
    mix_ECparms(:,:,mixcnt) = parms;                    % cluster generator parameters
    if MatrixProblem(mix_covrmats(:,:,mixcnt)) == 0
        novars(mixcnt) = 0; % indicate that this cluster has estimable variance
    end
end     % mixtures loop

densities = zeros(n,kcnt);
for mixcnt = 1:kcnt    % estimate pi*pdf(data,muk,sigmak)
    if (novars(mixcnt) == 0)
        [loglike,densities(:,mixcnt)] = MECPDF(data,type,mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),mix_ECparms(:,:,mixcnt));
        densities(:,mixcnt) = mix_propors(mixcnt)*densities(:,mixcnt);
    end    
end     % mixtures loop
densities = sum(densities,2);
%densities(densities == 0) = eps;
loglike = sum(log(densities));
% number params per mixture: p means, p(p+1)/2 variances/covariances, 1 beta/nu, 1 mix proportion
penalty = 3*(kcnt*p + kcnt*p*(p+1)/2 + kcnt*1 + (kcnt-1));
AIC = -2*loglike + penalty;

%{
JAH 20080924, checked for octave 3.4.3

Copyright (C) 2008 Prof. Hamparsum Bozdogan & J. Andrew Howe

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