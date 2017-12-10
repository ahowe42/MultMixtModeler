function TotRMDist = KMix_TotInClustRMDist(data, labels, smth, regul_scale)
% total distance = KMix_TotInClustRMDist(data matrix, mixture assignments, ...
%  covariance smoother, regularizing scale)
%  This will compute the total regularized mahalanobis distance across all
%  mixtures in a mixture model.
%
%  Where
%  data matrix --- (nxp) data matrix
%  mixture assignments --- n-valued vector of assignments 1:k
%  covariance smoother --- alpha code to pass covsmooth
%  Regularizing Scale --- if Inf, det(sigma)^c is replaced with C1(sigma)
%  total distance --- Total within cluster RM distance for mixture model

[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 4)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMix_TotInClustRMDist: INVALID USAGE-Please read the following instructions!\n'), help KMix_TotInClustRMDist, return
end

ks = unique(labels); kcnt = length(ks); % must do this in case some k missing (shouldn't)
novars = ones(1,kcnt);
mix_meanvecs = zeros(1,p,kcnt); mix_covrmats = zeros(p,p,kcnt);
rm_scales = zeros(1,kcnt);

% compute the cluster statistics
for mixcnt = 1:kcnt
    ind = (labels == ks(mixcnt));
    if (sum(ind) <= 1); continue; end;  % empty or singleton cluster - no variance
    % estimate the parameters
    mix_meanvecs(:,:,mixcnt) = mean(data(ind,:),1);     % cluster means
    mix_covrmats(:,:,mixcnt) = CovSmooth(data(ind,:),smth,1,0,sum(ind));
    if not(isequal(mix_covrmats(:,:,mixcnt),zeros(p)))
        novars(mixcnt) = 0; % indicate that this cluster does not have no variance
    else
        continue;
    end
    % compute the scaling portion of the RM
    if regul_scale == Inf
        rm_scales(mixcnt) = EntComp(mix_covrmats(:,:,mixcnt),0);
    else
        rm_scales(mixcnt) = (det(mix_covrmats(:,:,mixcnt)))^regul_scale;
    end
    % invert the covariance matrix
    mix_covrmats(:,:,mixcnt) = inv(mix_covrmats(:,:,mixcnt));
end     % mixtures loop

% compute the rm distances
TotRMDist = 0;
for datcnt = 1:n
    ind = find(labels(datcnt) == ks);
    if (novars(ind) == 1)
        % cluster with no variance - distance is 0
    else
        meancorr = data(datcnt,:) - mix_meanvecs(:,:,ind);
        TotRMDist = TotRMDist + rm_scales(ind)*...
            (meancorr*mix_covrmats(:,:,ind)*(meancorr'));
    end
end     % datapoints loop

% JAH 20060105, adapted for octave 3.4.3 20120324
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
