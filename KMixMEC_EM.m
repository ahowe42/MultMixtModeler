function [posteriors,newpi,newmu,newsigma,newparms] = KMixMEC_EM(data,labels,EMparams,ECtype,smth,savename,init,pltflg)
% [posteriors, mixing proportions, mean vectors, covar matrices, type-specific params] = 
% KMixMEC_EM(data, labels, EM parameters, EC subclass, covariance smoother, figure name, init. method, plot flag)
%  This implements the standard EM algorithm for multivariate elliptically
%  contoured distribution mixtures.  It returns the posterior probabilities 
%  of group belonging and estimates of the parameters per cluster.  Pass in a
%  smoothing code for CovSmooth for cases in which a cluster's covariance
%  matrix is ill-conditioned.
%
%  Where
%  data --- (nxp) matrix of data
%  labels --- n-vector of mixtures assignments in range of 1:k
%  EM Parameters --- [convergence criteria, max iterations]
%  EC subclass ---  'KT' (Kotz) or 'PVII' (Pearson Type VII)
%  covariance smoother --- alpha code to pass covsmooth
%  figure name --- name to save figure, include full path; if you don't
%     want the plot saved as a file, pass in 0.
%  init. method --- string of initialization method (GARM,GKM,...)
%  plot flag --- 0 = show only at end; 1 = update on-the-fly
%  posteriors --- (nxk) matrix of belonging probabilities
%  mixing proportions --- (1xk) vector of mixing proportions
%  mean vectors --- (1,p,k) matrix of mean vector per cluster
%  covar matrices --- (p,p,k) matrix of covariance matrix per cluster
%  type-specific parms --- (1,2/3,k) vector of [N,r,beta] (KT) or [N,nu] PVII
%
%  JAH 20080925

[n,p] = size(data);
k = max(labels);    % requires at least 1 datapoint in all clusters from 1:max(labels)

if (nargin ~= 8) || (length(labels) ~= n) || (length(EMparams) < 2) || ...
        (not(isequal([1:k],unique(labels))))
    % wrong # of args, dimensional mismatch between data and labels, not all k included
    fprintf('KMixMEC_EM: INVALID USAGE-Please read the following instructions!\n'), help KMixMEC_EM, return
end

% extract EM parameters
convgcrit = EMparams(1); maxiter = EMparams(2)
% compute initial mean, covariance, and mixing proportion estimates from labels
ECparms = zeros(1,2+isequal(ECtype,'KT'),k);
meanvecs = zeros(1,p,k); covrmats = zeros(p,p,k); mixpropors = zeros(1,k);

for mixcnt = 1:k
    ind = (labels == mixcnt);
    mixpropors(mixcnt) = sum(ind)/n;
    [meanvecs(:,:,mixcnt),covrmats(:,:,mixcnt),ECparms(:,:,mixcnt)] = MECEstimate(data(ind,:),ECtype,smth,1);
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
    oldpost = posteriors;
    % E-step - estimation of posterior probabilities
    for mixcnt = 1:k    % estimate pi*pdf(data,muk,sigmak)
        if isequal(covrmats(:,:,mixcnt), zeros(p))
            % just let the posterior probability be 0
            posteriors(:,mixcnt) = 0;
        else
            [jnk,dens] = MECPDF(data,ECtype,meanvecs(:,:,mixcnt),covrmats(:,:,mixcnt),ECparms(:,:,mixcnt));
            posteriors(:,mixcnt) = mixpropors(mixcnt)*dens;
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
    % use MAP rule to get groups for estimation of beta & nu
    [v,labs] = max(posteriors,[],2);
    
    xp = zeros(n,k*p);              % xi*pik - used for mean and covariance estimation
    for mixcnt = 1:k                % reestimate mean vectors & covariance matrices ...
        if mixpropors(mixcnt) < eps % only if not 0 posterior probability
            meanvecs(:,:,mixcnt) = zeros(1,p);
            covrmats(:,:,mixcnt) = zeros(p);
        else            
            ind = (labs == mixcnt); nn = sum(ind);
            nk = ceil(n*mixpropors(mixcnt));
            % estimate means
            xp(:,[grps(mixcnt,:)]) = data.*repmat(posteriors(:,mixcnt),1,p);
            meanvecs(:,:,mixcnt) = mean(xp(:,[grps(mixcnt,:)]))/mixpropors(mixcnt);        
            % estimate covariances & generator parameters
            meancorr = data - repmat(meanvecs(:,:,mixcnt),n,1);
            W = 0;
            for datcnt = 1:n        % this sucks, have to loop through datapoints :-(
                W = W + posteriors(datcnt,mixcnt)*(meancorr(datcnt,:)'*meancorr(datcnt,:));
            end     % data loop
            if isequal(ECtype,'KT')
                if (nn < p)
                    ECparms(1,3,mixcnt) = 0.75;
                else
                    % estimate beta first
                    dsq = 0;
                    siginv = inv(CovSmooth(data(ind,:),smth,1,1,sum(ind)));
                    if siginv(1,1) == inf
                        disp('debug')
                        [st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1)]);
                    end
                    for datcnt = 1:n        % this sucks, have to loop through datapoints :-(
                        dsq = dsq + ind(datcnt)*(meancorr(datcnt,:)*siginv*meancorr(datcnt,:)')^2;
                    end     % data loop
                    dsq = dsq/nn;
                    betaest = @(b) abs(log((((p^2)*gamma(p./(2*b)).*gamma((p+4)./(2*b)))./(gamma((p+2)./(2*b)).^2)))-log(dsq));
                    [ECparms(1,3,mixcnt),fval,flag] = fminbnd(betaest,0.1,10);
                    % if no solution found, set halfway between Laplace & Gaussian
                    if (flag ~= 1) || (fval == NaN) || (fval == Inf) || (ECparms(1,3,mixcnt) == NaN) || (ECparms(1,3,mixcnt) <=0)
                        disp('Fminbnd failed to converge')
                        ECparms(1,3,mixcnt) = 0.75;
                    end
                end
                % now get N and r
                ECparms(1,1,mixcnt) = 1; parms(1,2,mixcnt) = 0.5;
                % finally, the MLE for Sigma
                covrmats(:,:,mixcnt) = CovSmooth(p*((nk*p)/ECparms(1,3,mixcnt))^(-1/ECparms(1,3,mixcnt))*W,smth,0,1,nk);
                tmp = inv(covrmats(:,:,mixcnt));
                if (tmp(1,1) == inf) || MatrixProblem(covrmats(:,:,mixcnt))
                    disp('debug')
                    [st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1)]);
                end
            elseif isequal(ECtype,'PVII')
                if (nn < p)
                    ECparms(1,2,mixcnt) = 5;
                else
                    % estimate nu first                    
                    meancorr = CovSmooth(data(ind,:),smth,1,1,nn);    % just reusing meancorr var
                    a = 3*sum(diag(meancorr).^2);
                    b = sum(sum(data(ind,:).^4))/nn;
                    ECparms(1,2,mixcnt) = min(100,2*(a-2*b)/(a-b));
                    % sometimes, the MoM estimator returns negative values, so we'll do a numerical search
                    if ECparms(1,2,mixcnt) < 0
                        t = data(ind,:) - repmat(meanvecs(:,:,mixcnt),nn,1);
                        titt = diag(t*inv(meancorr)*t');
                        % remove -(n/2)*log(det(samp_covr)), since it is constant w.r.t. nu
                        mvpviiloglike = @(nu) -(nn*log(gamma((p+nu)/2)/((pi*nu)^(p/2)*gamma(nu/2)))-(p+nu)*sum(log(1+titt/nu))/2);
                        [ECparms(1,2,mixcnt),maxll,flag,output] = fminbnd(mvpviiloglike,3,100);
                        % if no solution found, set to t with 5 degrees of freedom
                        if (flag ~= 1) || (maxll == NaN) || (maxll == Inf) || (ECparms(1,2,mixcnt) == NaN) || (ECparms(1,2,mixcnt) <=0)
                            disp('Fminbnd failed to converge')                
                            ECparms(1,2,mixcnt) = 5;
                        end                    
                    end
                end
                % now get N
                ECparms(1,1,mixcnt) = (p+ECparms(1,2,mixcnt))/2;
                % finally, the MLE for Sigma - same as for Gaussian
                covrmats(:,:,mixcnt) = CovSmooth(W/nk,smth,0,1,nk);
            end
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
    newparms = nan(size(ECparms));
    % finish up with the progress plot
    figure(fhga)
    plot([0:(itercnt-1)],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress: FAILED TO CONVERGE']);
else
    newpi = mixpropors; newmu = meanvecs; newsigma = covrmats; newparms = ECparms;
    % finish up with the progress plot
    figure(fhga)
    plot([0:itercnt],loglike,'bo-'); xlabel('Iteration');ylabel('Log-Likelihood');
    title(['EM(',init,') Progress']);
    
end
if not(isequal(savename,0))
    hgsave(fhga,[savename,'_EM']); close(fhga)
end
