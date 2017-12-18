function [kdens,hmat,nmparam] = MVKDE_Gauss(datamatrix,dataeval,htype,smth,ifprob)
%{
  [density estimates, bandwidth matrix, parameter count] = MVKDE_Gauss(data matrix,...
  eval data, bandwidth estimator, covariance smoother, problem only)
   Perform multivariate kernel density estimation using the Gaussian kernel
   function.  This function allows the use of a completely general
   bandwidth matrix (diagonal or not) that just must be positive definite.
   If an incorrect bandwidth estimator code is entered, 3 will be used.  If
   the bandwidth matrix is inestimable inf's will be returned.

   Where
   data matrix --- (nxp) matrix of data, can be univariate
   eval data --- (mxp) data for which density estimates should be evaluated
   bandwidth estimator --- numeric code indicating which bandwidth matrix estimator to use:
    ESTIMATOR                   SHAPE
    1 = [((p)^-1)*tr(sigma)]*I  spherical
    2 = diag(sigma)             ellipsoidal
    3 = W/n                     linear kernel
    4 = ((4/(p+2))^(2/(p+4)))*(n^(-2/(p+4)))*sigma
    5 = [(4/(p+2))^(1/(p+4))*(n^(-1/(p+4)))]*sqrt(diag(sigma))
    6 = [(4/(p+2))^(1/(p+4))*(n^(-1/(p+4)))]*C1(sigma)*I
    you can also pass in the pxp bandwidth matrix instead of the estimator type
   covariance smoother --- alpha code to pass covsmooth
   problem only --- 1 = instruct CovSmooth to only smooth if problem, 0 = always
   density estimates --- (mx1) vector of density estimates
   bandwidth matrix --- (pxp) matrix specifying amounts and directions of smoothing used
   parameter count --- number nonrestricted elements of H matrix

   dat = genrndmvnorm(100,2,[-5,5],[2,0.5;0.5,1]); MVKDE_Gauss(dat,[0,0],4,'MLE',1)

   See Also KDE, BivarKDE, ProdKDE, KDE_HComp, KDE_ICOMP, MPEProdKernelPDF.
   
  Copyright (C) 2007 J. Andrew Howe; see below
%}
   
   
[n,p] = size(datamatrix); [evan,evap] = size(dataeval);

if (isscalar(datamatrix)) || (p ~= evap)
    % data scalar, data eval not same dim as datamatrix
    fprintf('MVKDE_Gauss: INVALID USAGE-Please read the following instructions!\n'), help MVKDE_Gauss, return
end

if (sum(size(htype) == p) == 2) && ( p~= 1)
    % bandwidth matrix passed, not estimator type
    hmat = htype; nmparam = NaN;
else
    % invalid htype, assign to 3
    if (sum(htype == [1:6]) ~= 1); htype = 3; end;
    % estimate the covariance matrix: if ill-conditioned, use the smoother
    sigma = CovSmooth(datamatrix,smth,1,ifprob,n);
    W = (n-1)*sigma;  % unscale covariance to get sum((x-mu)'*(x-mu))
    switch htype
        case 1
            hmat = trace(sigma)*eye(p)/p;
            nmparam = 1;
        case 2
            hmat = diag(diag(sigma));
            nmparam = p;
        case 3
            hmat = W/n;
            nmparam = p*(p+1)/2;
        case 4
            hmat = ((4/(p + 2))^(2/(p + 4)))*(n^(-2/(p + 4)))*sigma;
            nmparam = p*(p+1)/2;
        case 5
            hmat = diag((4/(p+2))^(1/(p+4))*(n^(-1/(p+4)))*sqrt(diag(sigma)));
            nmparam = p;
        case 6
            hmat = diag((4/(p+2))^(1/(p+4))*(n^(-1/(p+4))))*EntComp(sigma,1)*eye(p);
            nmparam = 1;
    end
end

if MatrixProblem(hmat) ~= 0
    %disp('Bad Covariance Matrix')
    %[st,i] = dbstack; eval(['dbstop in ''MVKDE_Gauss.m'' at ',num2str(st(i).line+1)]);
    hmat = inf(p); nmparam = inf; kdens = inf(evan,1);
    return
end

negsrtdetH = det(hmat)^(-0.5);      % reciprocal square root determinant of H
invH = inv(hmat);                   % inverse of H
kdens = zeros(evan,n);

% loop and calc the estimated densities
for ncnt = 1:evan
    % if we repeat the datapoint n times, we can, in one step, compute the
    % density contribution of the entire datamatrix to this datapoint - it
    % is in the diagonals of diff*invH*diff'.
    %diff = repmat(dataeval(ncnt,:),n,1) - datamatrix;
    %kdens(ncnt,:) = negsrtdetH*diag(exp(-0.5*diff*invH*diff'));
    % this is equivalent to:
    % JAH 20080115 - though it seems the matrix method would be preferred,
    % for high p, it seems to be substantially slower
    for dcnt = 1:n
        diff = dataeval(ncnt,:) - datamatrix(dcnt,:);
        kdens(ncnt,dcnt) = negsrtdetH*exp(-0.5*diff*invH*diff');
    end
end                 % datapoints to evaluate loop
kdens = sum(kdens,2)/(n*(2*pi)^(p/2));

%{
JAH 20070219, adapted for octave 3.4.3 20120315

Copyright (C) 2007 J. Andrew Howe

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
