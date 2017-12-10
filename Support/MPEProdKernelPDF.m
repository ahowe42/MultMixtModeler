function [ED,H,betas] = MPEProdKernelPDF(data,dataeval,smth,H,betas)
% [densities,H,betas] = MPEProdKernelPDF(data,data eval,smoother,bandwidth vector, betas vector)
%  Compute the density estimates for some data using the product kernel with
%  the PE kernel function.
%
%  Where
%  data --- (nxp) data matrix
%  data eval --- (mxp) matrix of data at which to evaluate densities
%  smoother --- alpha code for covsmooth ('MLE' for none)
%  bandwidth vector --- optional (1xp) vector of bandwidth values to use
%  betas vector --- optional (1xp) vector of kurtosis parameters to use
%  densities --- (mx1) vector of probability densities
%  H --- (1xp) vector of bandwidth computed using Silerman's rule of Thumb (or just input passed back out)
%  betas --- (1xp) vector of beta estimates (or just input passed back out)
%
%  JAH 20081003

[n,p] = size(data); [ne,pe] = size(dataeval);

if (isscalar(data)) || (isscalar(dataeval)) || ((nargin~=3) && (nargin~=5)) || (p~=pe)
    fprintf('MPEProdKernelPDF: INVALID USAGE-Please read the following instructions!\n'), help MPEProdKernelPDF , return
end

if (nargin ~= 5)
    % get the covariance matrix, smoothed, if necessary
    sigmas = diag(CovSmooth(data,smth,1,1,n))';
    H = ((4/(n*(p+2)))^(1/(p+4)))*sigmas;

    % must loop through the dimensions and estimate the Betas
    % First: get the squared regularized mahalanobis distances (square of kernel of Gaussian)
    stand = ((data - repmat(mean(data),n,1)).^2)./repmat(sigmas,n,1);
    avgregmahalsq = mean(stand.^2);
    if sum(isinf(avgregmahalsq))>0
        disp('debug')
        [st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1)]);
    end
    if not(isequal(size(avgregmahalsq),[1,p]))
        disp('debug')
        [st,i] = dbstack; eval(['dbstop in ',st(i).name,' at ',num2str(st(i).line+1)]);
    end
    % Second: do numerical search to get betas
    betas = zeros(1,p);
    for pcnt = 1:p
        if avgregmahalsq(pcnt) > 0
            betaest = @(b) abs(log(((gamma(1/(2*b)).*gamma(5/(2*b)))/(gamma(3./(2*b)).^2)))-log(avgregmahalsq(pcnt)));
            [betas(pcnt),fval,flag] = fminbnd(betaest,0.1,10);
            if (flag ~= 1) || (fval == NaN) || (fval == Inf) || (betas(pcnt) == NaN) || (betas(pcnt) <=0)
                disp(sprintf('Fminbnd failed to converge: Dbarsq = %0.4f',avgregmahalsq(pcnt)))
                betas(pcnt) = 1;
            end
        else
            disp('Dbarsq = 0, some dimension has all same values')
            betas(pcnt) = 1;
        end        
    end                 % loop through dimensions :-(
end

% ok, now that we have H and betas, do the real work
Hs = repmat(H,n,1);             % make a matrix with H row for each data row
Bs = repmat(betas,n,1);         % make a matrix with beta row for each data row
K = @(t) exp(-0.5*t.^Bs)./(gamma(1+1./(2*Bs)).*(2.^(1+1./(2*Bs))));
ED = zeros(ne,1);
for ncnt = 1:ne
    t = ((data - repmat(dataeval(ncnt,:),n,1)).^2)./Hs; % compute the kernel of the kernel, hahaha
    ED(ncnt) = sum(prod(K(t),2));
end                 % datapoints loop
ED = ED/(n*prod(H));

% SLOWER LOOPING WAY
%K = @(t,B) exp(-0.5*t.^B)/(gamma(1+1/(2*B)).*(2^(1+1/(2*B))));
%ED = zeros(ne,1);
%for ncnt = 1:ne
%    t=1;
%    for pcnt = 1:p
%        t = t.*K(((data(:,pcnt) - data(ncnt,pcnt)).^2)/H(pcnt),betas(pcnt));
%    end             % variables loop
%    ED(ncnt) = sum(t);
%end                 % datapoints loop
%ED = ED/(ne*prod(H));