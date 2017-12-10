function smcov = CovSmooth(data_matrix,meth,swtch,ifprob,numobs,regul_amnt)
%{
  covariance = CovSmooth(data/covariance matrix, smoothing method code,
  data input switch, problem only, number observations, regularization amount)

  Compute a smooth or robust estimator of a covariance matrix, either from
  a covariance matrix or the original data.  If covsmooth is set to only
  smooth in case of a bad matrix , and the matrix is not bad, it just 
  returns the MLE. If no parameters are passed (covs = CovSmooth()), this
  returns a cell array with all the smoothing codes.

  Where
  data/covariance matrix --- either data or covariance
  smoothing method code --- string as shown below
  'MLE' = maximum likelihood estimator
  'MAXENT' = maximum entropy - Fiebig 1982 - DON'T USE IF SWITCH = 0
  'MLE/EB' = maximum likelihood / empirical bayes
  'MXETEB' = maximum entropy / empirical bayes - DON'T USE IF SWITCH = 0
  'STIRDG' = stipulated ridge - Shurygin 1983
  'STIDAG' = stipulated diagonal - Shurygin 1983
  'CONSUM' = convex sum - Chen 1976
  'LEDWLF' = ledoit wolf - Ledoit & Wolf 2003
  'THOMAZ' = Thomaz eigen stabilization Thomaz 2004
  'RDGREG' = ridge regularization, requires 6th parameter (alpha)
  data input switch --- 1 = data, 0 = covariance
  problem only --- 1 = If not MatrixProblem, don't smooth; 0 = always
  number observations --- salar number of observations in dataset
  regularization amount --- for sigma +alphaI, alpha (1x1)
  covariance --- smoothed covariance matrix, if appropriate

  Example: CovSmooth(rand(4),'CONSUM',0,1,100)

  See Also MatrixProblem

  Copyright (C) 2006 J. Andrew Howe; see below
%}

% if no parameters, pass out list of choices
if nargin == 0
    smcov = {'MLE';'MAXENT';'MLE/EB';'MXETEB';'STIRDG';'STIDAG';'CONSUM';...
        'LEDWLF';'THOMAZ';'RDGREG'};
    return;
end

meth = upper(meth); [d1,d2] = size(data_matrix);

if ((swtch == 0) && d1 ~= d2) || ((isequal(meth,'I')) && nargin ~= 6)
    % input a vector; if switch says covariance, not square, ridge regul but no amount
    fprintf('CovSmooth: INVALID USAGE-Please read the following instructions!\n'), help CovSmooth, return
end

% if number rows is 1 . . .
if (d1 == 1) && (swtch == 1)
    smcov = zeros(d2); return;      % and input is data, pass back 0
elseif (d1 == 1) && (swtch == 0)
    smcov = data_matrix; return;    % and input is cov (probably scalar), pass back unaltered cov
end

% compute covariance matrix if data passed
if swtch == 1
    cov_mat = cov(data_matrix,1); [numobs,p] = size(data_matrix);
else
    cov_mat = data_matrix; p = length(data_matrix);
end

% ensure there is actually some complexity, before working
if (isequal(cov_mat,zeros(p)) || isscalar(cov_mat))
    smcov = cov_mat; return;
end

% only proceed if told to always smooth, or told to smooth if there's a problem, and there is
if (ifprob == 1) && (MatrixProblem(cov_mat) == 0); smcov = cov_mat; return; end;

cmtrace = trace(cov_mat); eyep = eye(p);

switch(meth)
    case {'A','MAXENT'}
        [XX,IX] = sort(data_matrix);
        XX = [XX(1,:) ; XX ; XX(numobs,:)];

        XX = (XX(1:numobs+1,:) + XX(2:numobs+2,:))./2;
        DX = XX(2:numobs+1,:) - XX(1:numobs,:);
        DX = sum(DX.^2);
        DX = DX./(12*numobs);

        EX = (XX(2,:) - XX(1,:)).^2 + (XX(numobs+1,:) - XX(numobs,:)).^2;
        EX = EX./(6*numobs);

        D1 = diag(DX+EX);

        XX = ( XX(1:numobs,:) + XX(2:numobs+1,:) )./2;

        % Reorder
        for j = 1:p
            XX(:,j) = XX(IX(:,j),j); 
        end

        XX = XX - ones(numobs,1)*mean(XX);
        XX = (XX'*XX)./numobs;

        smcov = XX + D1;
    case {'B','MLE/EB'}
        smcov = cov_mat + ((p-1)/(numobs*trace(inv(cov_mat))))*eyep;
    case {'C','MXETEB'}    % this is a composite of B on A, so recursively get A, then do B on it
        cov_mat = CovSmooth(data_matrix,'A',swtch,ifprob,numobs);
        smcov = cov_mat + ((p-1)/(numobs*trace(cov_mat)))*eyep;
    case {'D','STIRDG'}
        smcov = cov_mat + p*(p-1)/(2*numobs*trace(inv(cov_mat)))*eyep;
    case {'E','STIDAG'}
        %ste = diag(diag(cov_mat).^(-0.5)); r = ste'*cov_mat*ste;
        r = StoR(cov_mat);  % updated to use StoR JAH 20080108
        pi = p*(p-1)/(2*numobs*(trace(inv(r))-p));
        smcov = (1 - pi)*cov_mat + pi*diag(diag(cov_mat));
    case {'F','CONSUM'}
        beta = (cmtrace^2)/trace(cov_mat^2);
        m = round(mean([0,2*(p*(1+beta) - 2)/(p-beta)]));
        smcov = (numobs/(numobs+m))*cov_mat + (1-(numobs/(numobs+m)))*(cmtrace/p)*eyep;
    case {'G','LEDWLF'}
        smcov = ledoit(data_matrix);    % this is not my code
    case {'H','THOMAZ'}
        [vecs,vals] = eig(cov_mat);
        vecs = real(vecs); vals = real(diag(vals));
        lambdastar = diag(max([vals,repmat(mean(vals),p,1)],[],2));
        smcov = (vecs*lambdastar*vecs');
    case {'I','RDGREG'}
        smcov = cov_mat + regul_amnt*eyep;
    otherwise
        smcov = cov_mat;
end

%{
JAH 20060214, Brant Quinton April 2010, JAH adapted for octave 3.4.3 20120305

Copyright (C) 2006 J. Andrew Howe

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
