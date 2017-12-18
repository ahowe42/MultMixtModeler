function ICOMP = KMixGauss_ICOMPPERF(data, labels, smth, ifprob)
%{
  ICOMP = KMixGauss_ICOMPPERF(data matrix, mixture assignments, ...
  covariance smoother, problem only)
   This will compute ICOMP(Performance) for a mixture of normals model.
   Note that it requires a global vector 'ys' holding the actual class
   labels for the sigma_hat = (1/n)sum((y-yhat)^2) computation.  If first
   element of ys is NaN, the function will return Inf;

   Where
   data matrix --- (nxp) data matrix
   mixture assignments --- n-valued vector of assignments 1:k
   covariance smoother --- alpha code to pass covsmooth
   problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
   ICOMP --- ICOMP(IFIM) for gaussian mixture model

  Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}

[n,p] = size(data); global ys;

if (length(labels) ~= n) || (nargin ~= 4)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMixGauss_ICOMPPERF: INVALID USAGE-Please read the following instructions!\n'), help KMixGauss_ICOMPPERF, return
end

if isnan(ys(1)); ICOMP = Inf; return; end;

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
mix_covrmats = zeros(p,p,kcnt);

for mixcnt = 1:kcnt
    ind = (labels == ks(mixcnt));
    if (sum(ind) <= 1); continue; end;  % empty or singleton cluster - no variance    
    % estimate the parameters    
    mix_covrmats(:,:,mixcnt) = CovSmooth(data(ind,:),smth,1,ifprob,sum(ind));
    if MatrixProblem(mix_covrmats(:,:,mixcnt)) ~= 0
        mix_covrmats(:,:,mixcnt) = 0;
    end
end     % mixtures loop

sigmasq_hat = sum((ys-labels).^2)/n;
if sigmasq_hat == 0
    % perfect match, set loglikelihood to 0?
    loglike = 0;
else
    loglike = (n*log(2*pi) + n*log(sigmasq_hat) + n)/(-2);
end
penalty = 2*EntComp(sum(mix_covrmats,3),0);
% finish ICOMP
ICOMP = -2*loglike + penalty;

%{
JAH 20070130, checked for octave 3.4.3

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