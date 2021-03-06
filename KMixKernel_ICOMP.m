function ICOMP = KMixKernel_ICOMP(data,labels,smth,htype,ifprob)
%{ 
  ICOMP = KMixKernel_ICOMP(data, labels, covariance smoother, ...
  bandwidth estimator, problem only)
   Compute ICOMP for a mixture of kernel density estimators model.

   Where
   data --- (nxp) matrix of data
   labels --- n-vector of mixtures assignments in range of 1:k
   covariance smoother --- alpha code to pass covsmooth
   bandwidth estimator --- code for MVKDE_Gauss for bandwidth matrix creation
   problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
   ICOMP --- value of ICOMP(IFIM) for a mixture of kernels model

  Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}

[n,p] = size(data);

if (nargin ~= 5) || (length(labels) ~= n)
    % wrong # of args, dimensional mismatch between data and labels
    fprintf('KMixKernel_ICOMP: INVALID USAGE-Please read the following instructions!\n'), help KMixKernel_ICOMP, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
posteriors = zeros(n,kcnt); hs = zeros(p,p,kcnt); slogknum = 0; ICOMP_mid = 0;
slogdetcov = 0; mix_hs = zeros(p,p,kcnt); mix_sqdetcov = zeros(1,kcnt);

for mixcnt = 1:kcnt                     % estimate pi*kernel(data,h)
    clust_ys = (labels == ks(mixcnt));  % datapoints in this mixture
    clust_no = (labels ~= ks(mixcnt));  % datapoints in this mixture - not!
    pik = sum(clust_ys)/n;
    if (sum(clust_ys) > 1) && not(isequal(cov(data(clust_ys,:)), zeros(p)))
        % get densities for datapoints in this cluster
        [f,hs(:,:,mixcnt),hest] = MVKDE_Gauss(data(clust_ys,:),data(clust_ys,:),htype,smth,ifprob);
        posteriors(clust_ys,mixcnt) = pik*f;
        % get densities for datapoints not in this cluster
        if (sum(clust_no) >= 1)
            posteriors(clust_no,mixcnt) = pik*MVKDE_Gauss(data(clust_ys,:),...
                data(clust_no,:),hs(:,:,mixcnt),smth,ifprob);
        end
    else
        continue;
    end
    
    if hs(:,:,mixcnt) == zeros(p)
        continue;
    end
    
    % computations with parameter estimates for ICOMP
    slogknum = slogknum + log(pik*n);   % sum(log(n_k))
    mix_trccov = trace(hs(:,:,mixcnt)); % trace of bandwidths
    mix_sqdetcov(mixcnt) = det(hs(:,:,mixcnt)); % determinate of bandwidths
%    mix_hs(:,:,mixcnt) = inv(hs(:,:,mixcnt));  % cluster bandwidths inverse
    % work on ICOMP    
    ICOMP_mid = ICOMP_mid + (1/pik)*mix_trccov + (trace(hs(:,:,mixcnt)^2) + ...
        mix_trccov^2)/2 + sum(diag(hs(:,:,mixcnt)).^2);
    slogdetcov = slogdetcov + log(mix_sqdetcov(mixcnt));
%    mix_sqdetcov(mixcnt) = 1/sqrt(mix_sqdetcov(mixcnt)); % inverse sqrt of det of bandwidths
end                 % mixtures loop
loglike  = sum(log(sum(posteriors,2))); % log-likelihood is sum(log(pi*kernel(data,h)))

% finalize ICOMP
ICOMP_mid = log(ICOMP_mid/(kcnt*hest));
penalty = kcnt*hest*ICOMP_mid - ((p + 2)*slogdetcov - p*slogknum) - (kcnt*p)*log(2*n);
ICOMP = -2*loglike + penalty;

%{
JAH 20070220, checked for octave 3.4.3

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