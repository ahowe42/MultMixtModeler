function [posteriors,newpi,newmu,newsigma,newbetas] = KMixPEPK_EM(data,labels,EMparams,smth,savename,init,pltflg)
% [posteriors, mixing proportions, mean vectors, covar matrices, betas] = 
% KMixPEPK_EM(data, labels, EM parameters, covariance smoother, figure name, init. method, plot flag)
%  This implements the EM algorithm for a mixture model composed of a
%  mixture of a product kernel based on the power exponential kernel.  It 
%  returns the posterior probabilities of group belonging and estimates of
%  the parameters per cluster.  Pass in a smoothing code for CovSmooth for 
%  cases in which a cluster's covariance matrix is ill-conditioned.
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
%  betas --- (1,p,k) matrix of betas vector per cluster
%
%  JAH 20081005

[n,p] = size(data);
k = max(labels);    % requires at least 1 datapoint in all clusters from 1:max(labels)

if (nargin ~= 7) || (length(labels) ~= n) || (length(EMparams) < 2) || (not(isequal([1:k],unique(labels))))
    % wrong # of args, dimensional mismatch between data and labels, not all k included
    fprintf('KMixPEPK_EM: INVALID USAGE-Please read the following instructions!\n'), help KMixPEPK_EM, return
end

% extract EM parameters
convgcrit = EMparams(1); maxiter = EMparams(2);

% compute initial parameter estimates from labels
meanvecs = zeros(1,p,k); covrmats = zeros(p,p,k); mixpropors = zeros(1,k);
betavecs = zeros(1,p,k); hvecs = zeros(1,p,k);
for mixcnt = 1:k
    ind = (labels == mixcnt);
    if sum(ind) < 2      % no variance
        mixpropors(mixcnt) = 0;
        % don't need to init anything else (or at all, really), since all currently 0
    else
        mixpropors(mixcnt) = sum(ind)/n;
        meanvecs(:,:,mixcnt) = mean(data(ind,:),1);
        covrmats(:,:,mixcnt) = diag(diag(cov(data(ind,:))));
        hvecs(:,:,mixcnt) = (4/(sum(ind)*(p+2)))^(1/(p+4))*diag(covrmats(:,:,mixcnt));
        meancorr = data(ind,:) - repmat(meanvecs(:,:,mixcnt),sum(ind),1);
        for pcnt = 1:p
            dsq = sum(meancorr(:,pcnt).^4)/(sum(ind)*covrmats(pcnt,pcnt,mixcnt)^2);
            if dsq == 0
                disp('here')
            end
            betaest = @(b) abs(log(((gamma(1./(2*b)).*gamma(5./(2*b)))./(gamma(3./(2*b)).^2)))-log(dsq));
            [betavecs(1,pcnt,mixcnt),fval,flag] = fminbnd(betaest,0.1,10);
            if (flag ~= 1) || (fval == NaN) || (fval == Inf) || (betavecs(1,pcnt,mixcnt) == NaN) || (betavecs(1,pcnt,mixcnt) <=0)
                disp('Fminbnd failed to converge')
                betavecs(1,pcnt,mixcnt) = 0;
            end
        end     % dimensions loop
        % if still ill-conditioned, make all zeros
        if MatrixProblem(covrmats(:,:,mixcnt)) ~= 0; covrmats(:,:,mixcnt) = 0; end;
    end
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
    for mixcnt = 1:k    % estimate pi*pdf(data)
        if isequal(covrmats(:,:,mixcnt), zeros(p))
            % just let the posterior probability be 0
            posteriors(:,mixcnt) = 0;
        else
            clust_ys = (labels == mixcnt);  % datapoints in this mixture
            clust_no = (labels ~= mixcnt);  % datapoints in this mixture - not!
            % get densities for datapoints in this cluster
            if sum(clust_ys) > 1
                posteriors(clust_ys,mixcnt) = mixpropors(mixcnt)*MPEProdKernelPDF(...
                    data(clust_ys,:),data(clust_ys,:),smth,hvecs(:,:,mixcnt),betavecs(:,:,mixcnt));
            end
            % get densities for datapoints not in this cluster
            if (sum(clust_no) >= 1) && (sum(clust_ys) > 1)
                posteriors(clust_no,mixcnt) = mixpropors(mixcnt)*MPEProdKernelPDF(...
                    data(clust_ys,:),data(clust_no,:),smth,hvecs(:,:,mixcnt),betavecs(:,:,mixcnt));
            end
        end
    end         % mixtures loop
    meancorr = sum(posteriors,2);                    % just reusing meancorr var - this has nothing to do with xi - muk
    loglike  = [loglike,sum(log(meancorr))];         % log-likelihood is sum(log(pi*pdf(data)))    
    posteriors = posteriors./repmat(meancorr,1,k);   % posterior probability of group membership
    [val,labels] = max(posteriors,[],2);             % get the current mixture assignments
    
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
    
    % M-step - estimation of the mixture parameters
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
            % estimate the rest
            meancorr = data - repmat(meanvecs(:,:,mixcnt),n,1);
            for pcnt = 1:p
                % variances
                covrmats(pcnt,pcnt,mixcnt)=sum(posteriors(:,mixcnt).*meancorr(:,pcnt).^2)/nk;
                % betas
                dsq = sum(posteriors(:,mixcnt).*meancorr(:,pcnt).^4)/(nk*covrmats(pcnt,pcnt,mixcnt)^2);
                betaest = @(b) abs(log(((gamma(1./(2*b)).*gamma(5./(2*b)))./(gamma(3./(2*b)).^2)))-log(dsq));
                [betavecs(1,pcnt,mixcnt),fval,flag] = fminbnd(betaest,0.1,10);
                if (flag ~= 1) || (fval == NaN) || (fval == Inf) || (betavecs(1,pcnt,mixcnt) == NaN) || (betavecs(1,pcnt,mixcnt) <=0)
                    disp('Fminbnd failed to converge')
                    betavecs(1,pcnt,mixcnt) = 0;
                end
            end     % dimensions loop
            % if still ill-conditioned, make all zeros
            if MatrixProblem(covrmats(:,:,mixcnt)) ~= 0; covrmats(:,:,mixcnt) = 0; end;
            % now get H
            hvecs(:,:,mixcnt) = (4/(n*mixpropors(mixcnt)*(p+2)))^(1/(p+4))*diag(covrmats(:,:,mixcnt));
        end
    end             % mixtures loop    
    itercnt = itercnt + 1;
end                 % iterations loop

if itercnt > maxiter
    disp(sprintf('Maximum Number of Iterations %0.0f Achieved',maxiter))
    posteriors = nan(n,k);
    newpi = nan(1,k);
    newmu = nan(size(meanvecs));
    newsigma = nan(size(covrmats));
    newbetas = nan(size(betavecs));
    % finish up with the progress plot
    figure(fhga)
    plot([0:(itercnt-1)],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress: FAILED TO CONVERGE']);
else
    newpi = mixpropors; newmu = meanvecs; newsigma = covrmats; newbetas = betavecs;
    % finish up with the progress plot
    figure(fhga)
    plot([0:itercnt],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress']);    
end
if not(isequal(savename,0))
    hgsave(fhga,[savename,'_EM']); close(fhga)
end
