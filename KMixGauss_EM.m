function [posteriors,newpi,newmu,newsigma] = KMixGauss_EM(data,labels,EMparams,smth,savename,init,pltflg)
% [posteriors, mixing proportions, mean vectors, covar matrices] = KMixGauss_EM(
% data, labels, EM parameters, covariance smoother, figure name, init. method, plot flag)
%  This implements the standard EM algorithm for Gaussian mixtures.  It
%  returns the posterior probabilities of group belonging and estimates
%  of the mean vector and covariance matrix per cluster.  Pass in a
%  smoothing code for CovSmooth for cases in which a cluster's covariance
%  matrix is ill-conditioned.
%
%  Where
%  data --- (nxp) matrix of data
%  labels --- n-vector of mixtures assignments in range of 1:k
%  EM Parameters --- [convergence criteria, max iterations]
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
%  JAH 20061026
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT

[n,p] = size(data);
k = max(labels);    % requires at least 1 datapoint in all clusters from 1:max(labels)

if (nargin ~= 7) || (length(labels) ~= n) || (length(EMparams) < 2) || (not(isequal([1:k],unique(labels))))
    % wrong # of args, dimensional mismatch between data and labels, not all k included
    fprintf('KMixGauss_EM: INVALID USAGE-Please read the following instructions!\n'), help KMixGauss_EM, return
end

% extract EM parameters
convgcrit = EMparams(1); maxiter = EMparams(2);

% compute initial mean, covariance, and mixing proportion estimates from labels
meanvecs = zeros(1,p,k); covrmats = zeros(p,p,k); mixpropors = zeros(1,k);
for mixcnt = 1:k
    ind = (labels == mixcnt);
    mixpropors(mixcnt) = sum(ind)/n;
    meanvecs(:,:,mixcnt) = mean(data(ind,:),1);
    % estimate the covariance matrix: if ill-conditioned, use the smoother
    covrmats(:,:,mixcnt) = CovSmooth(data(ind,:),smth,1,1,sum(ind));
    % if still ill-conditioned, make all zeros
    if MatrixProblem(covrmats(:,:,mixcnt)) ~= 0; covrmats(:,:,mixcnt) = 0; end;
end                 % mixtures loop
posteriors = zeros(n,k); loglike = []; itercnt = 0; grps = reshape([1:(k*p)],p,k)';
fhga = figure;

if convgcrit < 1
    llstr = sprintf('%%0.%0.0ff',max(3,ceil(sqrt(log(1/convgcrit)))));
else
    llstr = '%0.3f';
end

while itercnt <= maxiter
    % E-step - estimation of posterior probabilities
    for mixcnt = 1:k    % estimate pi*pdf(data,muk,sigmak)
        if isequal(covrmats(:,:,mixcnt), zeros(p))
            % just let the posterior probability be 0
            posteriors(:,mixcnt) = 0;
        else
            posteriors(:,mixcnt) = mixpropors(mixcnt)*mvnpdf(data,meanvecs(:,:,mixcnt),covrmats(:,:,mixcnt));
        end
    end         % mixtures loop
    meancorr = sum(posteriors,2);                    % just reusing meancorr var - this has nothing to do with xi - muk
    loglike  = [loglike,sum(log(meancorr))];         % log-likelihood is sum(log(pi*pdf(data,muk,sigmak)))    
    posteriors = posteriors./repmat(meancorr,1,k);   % posterior probability of group membership
    
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
        plot([0:itercnt],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
        title(['EM(',init,') Progress']); drawnow
    end
    
    % M-step - estimation of the pi's, mu's, and sigma's
    mixpropors = mean(posteriors);  % reestimate mixing proportions
    
    xp = zeros(n,k*p);              % xi*pik - used for mean and covariance estimation
    for mixcnt = 1:k                % reestimate mean vectors & covariance matrices ...
        if mixpropors(mixcnt) == 0  % only if not 0 posterior probability
            meanvecs(:,:,mixcnt) = 0; covrmats(:,:,mixcnt) = 0;
        else
            nk = ceil(n*mixpropors(mixcnt));
            % estimate means
            xp(:,[grps(mixcnt,:)]) = data.*repmat(posteriors(:,mixcnt),1,p);
            meanvecs(:,:,mixcnt) = sum(xp(:,[grps(mixcnt,:)]))/nk;
            % estimate covariances
            meancorr = data - repmat(meanvecs(:,:,mixcnt),n,1); tmp = 0;
            for datcnt = 1:n        % this sucks, have to loop through datapoints :-(
                tmp = tmp + posteriors(datcnt,mixcnt)*(meancorr(datcnt,:)'*meancorr(datcnt,:));
            end     % data loop
            % estimate the covariance matrix: if ill-conditioned, use the smoother
            covrmats(:,:,mixcnt) = CovSmooth(tmp/nk,smth,0,1,nk);
            % if still ill-conditioned, make all zeros
            if MatrixProblem(covrmats(:,:,mixcnt)) ~= 0; covrmats(:,:,mixcnt) = 0; end;
        end
    end         % mixtures loop
    
    itercnt = itercnt + 1;
end             % iterations loop

if itercnt > maxiter
    disp(sprintf('Maximum Number of Iterations %0.0f Achieved',maxiter))
    posteriors = nan(n,k);
    newpi = nan(1,k);
    newmu = nan(size(meanvecs));
    newsigma = nan(size(covrmats));
    % finish up with the progress plot
    figure(fhga)
    plot([0:(itercnt-1)],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress: FAILED TO CONVERGE']);
else
    newpi = mixpropors; newmu = meanvecs; newsigma = covrmats;
    % finish up with the progress plot
    figure(fhga)
    plot([0:itercnt],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress']);
    
end
if not(isequal(savename,0))
    hgsave(fhga,[savename,'_EM']); close(fhga)
end
