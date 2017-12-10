function [posteriors,newpi,newmu,newsigma] = KMixKernel_EM(data,labels,EMparams,smth,savename,init,pltflg)
% [posteriors, mixing proportions, mean vectors, covar matrices] = KMixKernel_EM(
% data, labels, EM parameters, covariance smoother, figure name, plot flag, init. method)
%  This implements the standard EM algorithm for a mixture of kernel density
%  estimators.  It returns the posterior probabilities of group belonging 
%  and estimates of the mean vector, covariance matrix, and mixing proportion
%  per cluster.  Pass in a smoothing code for CovSmooth if you want to smooth
%  the covariance matrix; otherwise, pass in smth.
%
%  Where
%  data --- (nxp) matrix of data
%  labels --- n-vector of mixtures assignments in range of 1:k
%  EM Parameters --- [convergence criteria, max iterations, bandwidth estimator type]
%  covariance smoother --- alpha code to pass covsmooth
%  figure name --- name to save figure, include full path; if you don't
%     want the plot saved as a file, pass in 0.
%  init. method --- string of initialization method (GARM,GKM,...)
%  plot flag --- 0 = show only at end; 1 = update on-the-fly
%  posteriors --- (nxk) matrix of belonging probabilities
%  mixing proportions --- (1xk) vector of mixing proportions
%  mean vectors --- (1,p,k) matrix of mean vector per cluster
%  covar matrices --- (p,p,k) matrix of covariance matrix per cluster
%
%  JAH 20070215
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT

[n,p] = size(data);
k = max(labels);    % requires at least 1 datapoint in all clusters from 1:max(labels)

if (nargin ~= 7) || (length(labels) ~= n) || (length(EMparams) < 3) || (not(isequal([1:k],unique(labels))))
    % wrong # of args, dimensional mismatch between data and labels, not all k included
    fprintf('KMixKernel_EM: INVALID USAGE-Please read the following instructions!\n'), help KMixKernel_EM, return
end

% extract EM parameters & initialize stuff
convgcrit = EMparams(1); maxiter = EMparams(2); bandtype = EMparams(3);
loglike = []; itercnt = 0; mixpropors = zeros(1,k); fhga = figure;
for mixcnt = 1:k
    mixpropors(mixcnt) = sum(labels == mixcnt)/n;
end                 % mixtures loop

if convgcrit < 1
    llstr = sprintf('%%0.%0.0ff',max(3,ceil(sqrt(log(1/convgcrit)))));
else
    llstr = '%0.3f';
end

while itercnt <= maxiter
    % E-step - estimation of posterior probabilities
    posteriors = zeros(n,k);
    for mixcnt = 1:k    % estimate pi*kernel(data,h)
        clust_ys = (labels == mixcnt);  % datapoints in this mixture
        clust_no = (labels ~= mixcnt);  % datapoints in this mixture - not!
        % get densities for datapoints in this cluster
        if sum(clust_ys) > 1
            [f,h] = MVKDE_Gauss(data(clust_ys,:),data(clust_ys,:),bandtype,smth,1);
            posteriors(clust_ys,mixcnt) = mixpropors(mixcnt)*f;
        end
        % get densities for datapoints not in this cluster
        if (sum(clust_no) >= 1) && (sum(clust_ys) > 1)
            posteriors(clust_no,mixcnt) = mixpropors(mixcnt)*...
                MVKDE_Gauss(data(clust_ys,:),data(clust_no,:),h,smth,1);
        end
    end         % mixtures loop
    sumpost = sum(posteriors,2);
    loglike  = [loglike,sum(log(sumpost))];         % log-likelihood is sum(log(pi*kernel(data,h)))
    posteriors = posteriors./repmat(sumpost,1,k);   % posterior probability of group membership
    [val,labels] = max(posteriors,[],2);            % get the current mixture assignments
    
    % check for convergence if enough iterations done
    if (itercnt > 5) && (abs(loglike(end) - loglike([end-1])) <= convgcrit)
        % loglikelihood converged        
        disp(table2str({'Iteration';'Log-likelihood'},[[0:itercnt]',loglike'],{'%0.0f';llstr},0))
        disp(sprintf('Log-likelihood converged in iteration %0.0f',itercnt))
        break;
    end
    
    % SHOW PROGRESS PLOT
    if (pltflg == 1)
        figure(fhga)
        plot([0:itercnt],loglike,'bo-');xlabel('Iteration');ylabel('Log-Likelihood');
        title(['EM(',init,') Progress']); drawnow
    end
    
    % M-step - estimation of the pi's, sigma's, and h's
    mixpropors = mean(posteriors);  % reestimate mixing proportions
    % don't need to estimate sigma (since only used possibly for h) or h (both computed in kernel function)
    itercnt = itercnt + 1;
end             % iterations loop

if itercnt > maxiter
    disp(sprintf('Maximum Number of Iterations %0.0f Achieved',maxiter))
    posteriors = nan(n,k);
    newpi = nan(1,k);
    newmu = nan(1,p,k);
    newsigma = nan(p,p,k);
    % finish up with the progress plot
    figure(fhga)
    plot([0:itercnt],loglike,'bo-'); xlabel('Iteration'); ylabel('Log-Likelihood');
    title(['EM(',init,') Progress - FAILED TO CONVERGE']); drawnow
else
    % compute final parameter estimates
    newsigma = zeros(p,p,k); newpi = zeros(1,k); newmu = zeros(1,p,k);
    for mixcnt = 1:k
        ind = (labels == mixcnt);
        newpi(mixcnt) = sum(ind)/n;
        newmu(:,:,mixcnt) = mean(data(ind,:),1);
        if sum(ind) > 1    % ensure we store cov = zeros(p) for singleton clusters
            newsigma(:,:,mixcnt) = CovSmooth(data(ind,:),smth,1,1,sum(ind));
        end    
    end                 % mixtures loop
    % finish up with the progress plot
    figure(fhga)
    plot([0:itercnt],loglike,'bo-'); xlabel('Iteration'); ylabel('Log-Likelihood');
    title(['EM(',init,') Progress']); drawnow
end
if not(isequal(savename,0))
    hgsave([savename,'_EM']); close(fhga)
end
