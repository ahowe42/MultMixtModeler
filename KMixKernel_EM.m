function [posteriors,newpi,newmu,newsigma] = KMixKernel_EM(data,labels,EMparams,bandtype,smth,savename,init,pltflg)
%{
 [posteriors, mixing proportions, mean vectors, covar matrices] = KMixKernel_EM(
  data, labels, EM parameters, bandwidth estimator, covariance smoother, figure name, plot flag, init. method)
   This implements the standard EM algorithm for a mixture of kernel density
   estimators.  It returns the posterior probabilities of group belonging 
   and estimates of the mean vector, covariance matrix, and mixing proportion
   per cluster.  Pass in a smoothing code for CovSmooth if you want to smooth
   the covariance matrix; otherwise, pass in smth.
 
   Where
   data --- (nxp) matrix of data
   labels --- n-vector of mixtures assignments in range of 1:k
   EM Parameters --- [convergence criteria, max iterations]
   bandwidth estimator --- code for MVKDE_Gauss for bandwidth matrix creation
   covariance smoother --- alpha code to pass covsmooth
   figure name --- name to save figure, include full path; if you don't
      want the plot saved as a file, pass in 0.
   init. method --- string of initialization method (GARM,GKM,...)
   plot flag --- 0 = show only at end; 1 = update on-the-fly
   posteriors --- (nxk) matrix of belonging probabilities
   mixing proportions --- (1xk) vector of mixing proportions
   mean vectors --- (1,p,k) matrix of mean vector per cluster
   covar matrices --- (p,p,k) matrix of covariance matrix per cluster
 
Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe
%}

[n,p] = size(data);
k = max(labels);    % requires at least 1 datapoint in all clusters from 1:max(labels)

if (nargin ~= 8) || (length(labels) ~= n) || (length(EMparams) < 2) || (not(isequal([1:k],unique(labels))))
    % wrong # of args, dimensional mismatch between data and labels, not all k included
    fprintf('KMixKernel_EM: INVALID USAGE-Please read the following instructions!\n'), help KMixKernel_EM, return
end

% extract EM parameters & initialize stuff
convgcrit = EMparams(1); maxiter = EMparams(2);
loglike = []; itercnt = 0; mixpropors = zeros(1,k); fhga = figure;
Hmats = zeros(p,p,k);
for mixcnt = 1:k
    mixpropors(mixcnt) = sum(labels == mixcnt)/n;
end                 % mixtures loop

if convgcrit < 1
    llstr = sprintf('%%0.%0.0ff',max(3,ceil(sqrt(log(1/convgcrit)))));
else
    llstr = '%0.3f';
end

grps = reshape([1:(k*p)],p,k)';
while itercnt <= maxiter
    % E-step - estimation of posterior probabilities
    posteriors = zeros(n,k);
    for mixcnt = 1:k    % estimate pi*kernel(data,h)
        clust_ys = (labels == mixcnt);  % datapoints in this mixture
        clust_no = (labels ~= mixcnt);  % datapoints in this mixture - not!
        % get densities for datapoints in this cluster
        if sum(clust_ys) > 1
            % if first time, estimate H, and return it, otherwise, just use mine
            if (itercnt == 0)
                [f,Hmats(:,:,mixcnt)] = MVKDE_Gauss(data(clust_ys,:),data(clust_ys,:),bandtype,smth,1);
            else
                [f,tmp] = MVKDE_Gauss(data(clust_ys,:),data(clust_ys,:),Hmats(:,:,mixcnt),smth,1);
            end
            posteriors(clust_ys,mixcnt) = mixpropors(mixcnt)*f;
        end
        % get densities for datapoints not in this cluster
        if (sum(clust_no) >= 1) && (sum(clust_ys) > 1)
            posteriors(clust_no,mixcnt) = mixpropors(mixcnt)*...
                MVKDE_Gauss(data(clust_ys,:),data(clust_no,:),bandtype,smth,1);
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
    
    xp = zeros(n,k*p);              % xi*pik - used for mean and covariance estimation
    for mixcnt = 1:k                % reestimate covariance matrices ...
        if mixpropors(mixcnt) == 0  % only if not 0 posterior probability
            Hmats(:,:,mixcnt) = 0;
        else
            nk = ceil(n*mixpropors(mixcnt));
            % estimate means
            xp(:,[grps(mixcnt,:)]) = data.*repmat(posteriors(:,mixcnt),1,p);
            % estimate covariances
            meancorr = data - repmat(sum(xp(:,[grps(mixcnt,:)]))/nk,n,1); W = 0;
            for datcnt = 1:n        % this sucks, have to loop through datapoints :-(
                W = W + posteriors(datcnt,mixcnt)*(meancorr(datcnt,:)'*meancorr(datcnt,:));
            end     % data loop
            % estimate the covariance matrix: if ill-conditioned, use the smoother
            Hmats(:,:,mixcnt) = CovSmooth(W/nk,smth,0,1,nk);
            % if still ill-conditioned, make all zeros
            if MatrixProblem(Hmats(:,:,mixcnt)) ~= 0; Hmats(:,:,mixcnt) = 0; end;
            % now do H
            switch bandtype
                case 1  % 1 = [((pn)^-1)*tr(W)]*I = tr(W/(n*pi*p)*I
                    Hmats(:,:,mixcnt) = trace(Hmats(:,:,mixcnt)/p)*eye(p);
                case 2  % 2 = diag(W)/n = diag(W/(n*pi))
                    Hmats(:,:,mixcnt) = diag(Hmats(:,:,mixcnt));
                case 3  % 3 = W/n = W/(n*pi) - no change needed
                    %Hmats(:,:,mixcnt) = W/n;
                case 4  % 4 = n^(-1/(p + 4))]*sqrt(sigma_regul)
                    Hmats(:,:,mixcnt) = ((4/(p + 2))^(2/(p + 4)))*(nk^(-2/(p + 4)))*Hmats(:,:,mixcnt);
                case 5  % 5 = [(4/(p+2))^(1/(p+4))*((n*pi)^(-1/(p+4)))]*diag(sigma_regul)
                    Hmats(:,:,mixcnt) = diag((4/(p+2))^(1/(p+4))*(nk^(-1/(p+4)))*diag(Hmats(:,:,mixcnt)));
                case 6  % 6 = [(4/(2*p+1))^(1/(p+4))*(n^(-1/(p+4)))]*C1(sigma_regul)
                    Hmats(:,:,mixcnt) = (4/(2*p+1))^(1/(p+4))*(nk^(-1/(p+4)))*EntComp(Hmats(:,:,mixcnt),1)*eye(p);
            end
        end
    end         % mixtures loop
    
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
    plot([0:(itercnt-1)],loglike,'bo-'); xlabel('Iteration'); ylabel('Log-Likelihood');
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
    hgsave(fhga,[savename,'_EM']); close(fhga)
end

%{
JAH 20070215, checked for octave 3.4.3

Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe

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