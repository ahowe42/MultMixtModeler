function SBC = KMixGauss_SBC(data, labels, smth, ifprob)
% SBC = KMixGauss_SBC(data matrix, mixture assignments, ...
% covariance smoother, problem only)
%  This will compute SBC for a mixture of normals model.
%
%  Where
%  data matrix --- (nxp) data matrix
%  mixture assignments --- n-valued vector of assignments 1:k
%  covariance smoother --- alpha code to pass covsmooth
%  problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
%  SBC --- value of SBC for gaussian mixture model
%
%  JAH 200670106
%  Copyright Prof. Hamparsum Bozdogan & J. Andrew Howe
%  All rights reserved, see LICENSE.TXT

[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 4)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMixGauss_SBC: INVALID USAGE-Please read the following instructions!\n'), help KMixGauss_SBC, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
mix_meanvecs = zeros(1,p,kcnt); mix_covrmats = zeros(p,p,kcnt);
mix_propors = zeros(1,kcnt); novars = ones(1,kcnt);

for mixcnt = 1:kcnt
    ind = (labels == ks(mixcnt));
    if (sum(ind) <= 1); continue; end;  % empty or singleton cluster - no variance    
    % estimate the parameters
    mix_propors(mixcnt) = sum(ind)/n;                   % cluster proportion
    mix_meanvecs(:,:,mixcnt) = mean(data(ind,:),1);     % cluster means
    mix_covrmats(:,:,mixcnt) = CovSmooth(data(ind,:),smth,1,ifprob,sum(ind));
    if MatrixProblem(mix_covrmats(:,:,mixcnt)) == 0
        novars(mixcnt) = 0; % indicate that this cluster has estimable variance
    end
end     % mixtures loop

densities = zeros(n,kcnt);
for mixcnt = 1:kcnt    % estimate pi*pdf(data,muk,sigmak)
    if (novars(mixcnt) == 0)
        densities(:,mixcnt) = mix_propors(mixcnt)*mvnpdf(data,mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt));
    end    
end     % mixtures loop
densities = sum(densities,2);
%densities(densities == 0) = eps;
loglike = sum(log(densities));

penalty = (kcnt*p + kcnt*p*(p+1)/2 + (kcnt-1))*log(n);
SBC = -2*loglike + penalty;
