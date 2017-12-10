function AIC = KMixPEPK_AIC(data,labels,smth,ifprob)
% AIC = KMixPEPK_AIC(data, labels, covariance smoother,problem only)
%  Compute AIC for a mixture model based on a multivariate product kernel
%  with the PE kernel.
%
%  Where
%  data --- (nxp) matrix of data
%  labels --- n-vector of mixtures assignments in range of 1:k
%  covariance smoother --- alpha code to pass covsmooth
%  problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
%  AIC --- value of AIC for a mixture of kernels model
%
%  JAH 20081006

[n,p] = size(data);

if (nargin ~= 4) || (length(labels) ~= n)
    % wrong # of args, dimensional mismatch between data and labels
    fprintf('KMixPEPK_AIC: INVALID USAGE-Please read the following instructions!\n'), help KMixPEPK_AIC, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
posteriors = zeros(n,kcnt);

for mixcnt = 1:kcnt                     % estimate pi*kernel(data,h)
    clust_ys = (labels == ks(mixcnt));  % datapoints in this mixture
    clust_no = (labels ~= ks(mixcnt));  % datapoints in this mixture - not!
    pik = sum(clust_ys)/n;
    if (sum(clust_ys) > 1) && not(MatrixProblem(cov(data(clust_ys,:))))
        % get densities for datapoints in this cluster
        [f,H,betas] = MPEProdKernelPDF(data(clust_ys,:),data(clust_ys,:),smth);
        posteriors(clust_ys,mixcnt) = pik*f;
        % get densities for datapoints not in this cluster
        if (sum(clust_no) >= 1)
            posteriors(clust_no,mixcnt) = pik*MPEProdKernelPDF(data(clust_ys,:),...
                data(clust_no,:),smth,H,betas);
        end
    end
end                 % mixtures loop
loglike  = sum(log(sum(posteriors,2))); % log-likelihood is sum(log(pi*kernel(data,h)))

penalty = 3*(kcnt*3*p + (kcnt-1));  % (k-1) mixing proportions + mu, h, beta for all dimensions
AIC = -2*loglike + penalty;