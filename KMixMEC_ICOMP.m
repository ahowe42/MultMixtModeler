function ICOMP = KMixMEC_ICOMP(data,labels,smth,type,ifprob)
% ICOMP = KMixMEC_ICOMP(data matrix, mixture assignments, covariance smoother,...
% subtype, problem only)
%  This will compute ICOMP for a mixture of elliptically-contoured model.
%
%  Where
%  data matrix --- (nxp) data matrix
%  mixture assignments --- n-valued vector of assignments 1:k
%  covariance smoother --- alpha code to pass covsmooth
%  subtype --- 'KT' (Kotz) or 'PVII' (Pearson Type VII)
%  problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
%  ICOMP --- ICOMP(IFIM) for the EC mixture model
%
%  JAH 20081009

[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 5)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMixMEC_ICOMP: INVALID USAGE-Please read the following instructions!\n'), help KMixMEC_ICOMP, return
end

ks = unique(labels); kcnt = length(ks); % must do this because some k can be missing
mix_meanvecs = zeros(1,p,kcnt); mix_covrmats = zeros(p,p,kcnt);
mix_propors = zeros(1,kcnt); novars = ones(1,kcnt);
Q1 = zeros(p,0.5*p*(p+1));
IFIMp = 0.5*p*(p+3);
traIFIMs = 0; detIFIMs = 1;
%mix_ifims = zeros(IFIMp,IFIMp,kcnt);

if isequal(type,'KT')
    mix_ECparms = zeros(1,3,kcnt); % 3 parms for Kotz
elseif isequal(type,'PVII')
    mix_ECparms = zeros(1,2,kcnt); % 3 parms for Pearson VII
else
    fprintf('KMixMEC_ICOMP: TYPE INVALID-Please read the following instructions!\n'), help KMixMEC_ICOMP, return
end

Dp = DupMatrix(p); Dpt = Dp';
% estimate the parameters
for mixcnt = 1:kcnt
    ind = (labels == ks(mixcnt));
    if (sum(ind) <= 1); continue; end;  % empty or singleton cluster - no variance    
    mix_propors(mixcnt) = sum(ind)/n;                   % cluster proportion
    [mu,sigma,parms] = MECEstimate(data(ind,:),type,smth,ifprob);
    mix_meanvecs(:,:,mixcnt) = mu;                      % cluster means
    mix_covrmats(:,:,mixcnt) = sigma;                   % cluster covariances
    mix_ECparms(:,:,mixcnt) = parms;                    % cluster generator parameters
    if MatrixProblem(mix_covrmats(:,:,mixcnt)) == 0
        novars(mixcnt) = 0; % indicate that this cluster has estimable variance
    else
        continue;
    end
    is = inv(sigma); vis = is(:);
    % compute the IFIM for this group
    if isequal(type,'KT')
        Q4 = Dpt*(n/2*kron(is,is)-(1-mix_ECparms(1,1,mixcnt)-(mix_ECparms(1,3,mixcnt)-1)...
            *(mix_ECparms(1,3,mixcnt)+n*p*0.5-1))*vis*vis'/p^2)*Dp;
        Q2 = (((mix_ECparms(1,1,mixcnt)+n*p*0.5-1)/(mix_ECparms(1,2,mixcnt)*...
            mix_ECparms(1,3,mixcnt)))^(1/mix_ECparms(1,3,mixcnt)))*sigma/(n*n*p);
    else    % PVII
        Q4 = Dpt*(n/2*kron(is,is)-((mix_ECparms(1,1,mixcnt)*n^2)/(2*mix_ECparms(1,1,mixcnt)-p-n*p)^2)*vis*vis')*Dp;
        Q2 = (mix_ECparms(1,2,mixcnt)/(2*mix_ECparms(1,1,mixcnt)-n*p))*sigma/n;
    end
    traIFIMs = traIFIMs + trace(Q2)+trace(inv(Q4));
    detIFIMs = detIFIMs*det([Q2,Q1;Q1',inv(Q4)]);
%    mix_ifims(:,:,mixcnt) = [Q2,Q1;Q1',inv(Q4)];
end                 % mixtures loop

% compute the maximized likelihood for the entire mixture model
densities = zeros(n,kcnt);
for mixcnt = 1:kcnt    % estimate pi*pdf(data,thetak)
    if (novars(mixcnt) == 0)
        [loglike,densities(:,mixcnt)] = MECPDF(data,type,mix_meanvecs(:,:,mixcnt),...
            mix_covrmats(:,:,mixcnt),mix_ECparms(:,:,mixcnt));
    end    
end     % mixtures loop
densities = sum(densities,2);
densities(densities == 0) = eps;
loglike = sum(log(densities));

%% now build the IFIM :-( and compute the complexity
%IFIMP = kcnt*(1+IFIMp);
%IFIM = zeros(IFIMP);
%IFIM([1:kcnt],[1:kcnt]) = diag(1./mix_propors);
%% get ranges for inserting blocks
%rc = [kcnt:IFIMp:IFIMP]'; rc = [rc([1:(end-1)])+1,rc([2:end])];
%% now insert the blocks
%for mixcnt = 1:kcnt
%    IFIM([rc(mixcnt,1):rc(mixcnt,2)],[rc(mixcnt,1):rc(mixcnt,2)]) = mix_ifims(:,:,mixcnt);
%end                 % mixtures loop
%% now we can do the complexity - but IFIM is not square?
%penalty = 2*EntComp(IFIM,1);

% don't need to build it, since block diagonal :-)
traIFIMs = traIFIMs+sum(1./mix_propors);
detIFIMs = detIFIMs*prod(1./mix_propors);
rnk = kcnt*(1+IFIMp);
C1 = 0.5*rnk*log(traIFIMs/rnk)-0.5*log(detIFIMs);
% number params per mixture: p means, p(p+1)/2 variances/covariances, 1 beta/nu, 1 mix proportion
penalty = (kcnt*p + kcnt*p*(p+1)/2 + kcnt*1 + (kcnt-1)) + log(n)*C1;

% finish ICOMP
ICOMP = -2*loglike + penalty;
