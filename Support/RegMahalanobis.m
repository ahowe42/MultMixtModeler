function regmahaldist = RegMahalanobis(datapoint, meanvec, covarmat, smth, scal, n)
% regularized mahalanobis distance = RegMahalanobis(data point, mean vector, ...
% covariance matrix, smooth code, scale param, samp size)
%  This computes and returns the regularized mahalanobis distance given by 
%  D(x,mu) = ((det(sigma_reg))^c)*(x-mu)*'inv(sigma_reg)*(x-mu) or
%  C1(sigma_reg)*(x-mu)*'inv(sigma_reg)*(x-mu).
%
%  Where
%  data point ---(1xp) vector or scalar observation of interest
%  mean vector --- (1xp) vector of means of entire dataset
%  covariance matrix --- (pxp) covariance matrix of entire dataset
%  smooth meth --- String code for CovSmooth; don't use MAXENT or MXETEB. 
%     If ridge regularization is required, pass RDGREG and the value of 
%     alpha instead of the sample size.  If none is required, pass MLE.
%  scale param --- Scalar value used for power of determinant of sigma pass
%     Inf to scale with C1(sigma) instead.
%  samp size --- Scalar sample size (used for regularizing), or alpha if
%     smooth meth == RDGREG
%  regularized mahalanobis distance --- Distance from data point to mean.
%
%  Example: xmvn = randn(100,4); RegMahalanobis(xmvn(1,:),mean(xmvn),cov(xmvn),'MLE/EB',Inf,100)
%
% See Also EntComp, CovSmooth, MatrixProblem

smth = upper(smth);

if (nargin ~= 6)
    % not all 6 args
    fprintf('RegMahalanobis: INVALID USAGE-Please read the following instructions!\n'), help RegMahalanobis, return
end

if (~isscalar(covarmat))
    % smooth / regularize the covariance matrix
    if isequal(smth,'RDGREG')
        covarmat = covarmat + eye(size(covarmat,1))*n;
    elseif not(isequal(smth,'MLE'))
        covarmat = CovSmooth(covarmat,smth,0,0,n);
    end
end
invcovarmat = inv(covarmat);      % invert the covariance matrix

% compute the scaling portion of the RM
if scal == Inf
    preval = EntComp(covarmat,1);   % this was, erroneously, EntComp(covarmat,0) 20070404 JAH
else
    preval = (det(covarmat))^scal;
end

% compute the result
xmm = datapoint - meanvec;
regmahaldist = preval*(xmm*invcovarmat*(xmm'));

% JAH 20061007, adapted for octave 3.4.3 20120315
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
