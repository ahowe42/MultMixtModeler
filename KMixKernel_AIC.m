function AIC = KMixKernel_AIC(data,labels,smth,htype,ifprob)
% AIC = KMixKernel_AIC(data, labels, covariance smoother, ...
% bandwidth estimator, problem only)
%  Compute AIC for a mixture of kernel density estimators model.
%
%  Where
%  data --- (nxp) matrix of data
%  labels --- n-vector of mixtures assignments in range of 1:k
%  covariance smoother --- alpha code to pass covsmooth
%  bandwidth estimator --- code for MVKDE_Gauss for bandwidth matrix creation
%  problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
%  AIC --- value of AIC for a mixture of kernels model
%
%  JAH 20070220
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT

[n,p] = size(data);

if (nargin ~= 5) || (length(labels) ~= n)
    % wrong # of args, dimensional mismatch between data and labels
    fprintf('KMixKernel_AIC: INVALID USAGE-Please read the following instructions!\n'), help KMixKernel_AIC, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
posteriors = zeros(n,kcnt);

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
    end
end                 % mixtures loop
loglike  = sum(log(sum(posteriors,2))); % log-likelihood is sum(log(pi*kernel(data,h)))

penalty = 3*(kcnt*hest + (kcnt-1));  % (k-1) mixing proportions, and some # of H elements per group
AIC = -2*loglike + penalty;
