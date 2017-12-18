function AIC = KMixPEPK_AIC(data,labels,smth,ifprob)
%{
  AIC = KMixPEPK_AIC(data, labels, covariance smoother,problem only)
   Compute AIC for a mixture model based on a multivariate product kernel
   with the PE kernel.

   Where
   data --- (nxp) matrix of data
   labels --- n-vector of mixtures assignments in range of 1:k
   covariance smoother --- alpha code to pass covsmooth
   problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
   AIC --- value of AIC for a mixture of kernels model
 
  Copyright (C) 2008 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}

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

%{
JAH 20081006, checked for octave 3.4.3

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