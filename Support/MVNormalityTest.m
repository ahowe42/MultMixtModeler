function [kurt,skew] = MVNormalityTest(data,silent,alpha)
%{
  [kurt results, skew results] = MVNormalityTest(data, silent, alpha)

  Performs Mardia's tests for multivariate normality on a data sample.
  The tests for skewness and kurtosis will both be performed.  
  The hypothesis test for kurtosis is two-sided, with the test stat being
  distributed as N(0,1) under the null.  If the sample value is less than
  the theoretical, the data is platykurtic (flat with light tails, widely
  dispersed - think Uniform).  If the sample value is greater than the 
  theoretical, the data is leptokurtic (peaked with heavy tails - think 
  Laplace).  The skewness test stat is assumed to be chi-squared(p(p+1)(p+2)/6)
  distributed under the null, and is a one-sided test.

  Where
  data  --- (nxp) matrix of data to evaluate
  silent --- 1 = don't print results to screen or create plot; 0 = do
  alpha --- optional scalar significance level (default = 0.05)
  kurt results --- (5x1) vector with theoretical value, sample value, test
     statistic, critical value, p-value
  skew results --- (5x1) vector with theoretical value, sample value, test
     statistic, critical value, p-value

  Example: [K,S] = MVNormalityTest(genrndmvpexp(1000,2,[0;0],[1,0;0,1],1),0);
  Example: [K,S] = MVNormalityTest(genrndmvpexp(1000,2,[0;0],[1,0;0,1],1.5),0);
  Example: [K,S] = MVNormalityTest(genrndmvpexp(1000,2,[0;0],[1,0;0,1],0.75),0);

  Copyright (C) 2007 J. Andrew Howe
%}
  
if isscalar(data) || ((nargin ~= 2) && (nargin ~= 3))
    % data should be a vector or matrix, wrong # of args
    fprintf('MVNormalityTest: INVALID USAGE-Please read the following instructions!\n'), help MVNormalityTest, return
end

% prepare
lin = [repmat('+-',1,15),'+'];
[n,p] = size(data);
samp_mean = mean(data);
samp_covr = cov(data);
%samp_covr = covsmooth(data,'MLE/EB',1,1,n);
Sinv = inv(samp_covr);
mc = (data - repmat(samp_mean,n,1));
samp_kurt = 0; samp_skew = 0;

if (nargin == 2); alpha = 0.05; end;
Zval = norminv(1 - alpha/2);
dof = p*(p + 1)*(p + 2)/6;
chi2val = chi2inv(1-alpha,dof);

% test for kurtosis
for datcnt = 1:n
    samp_kurt = samp_kurt + (mc(datcnt,:)*Sinv*mc(datcnt,:)').^2;
end                 % datapoints loop
samp_kurt = samp_kurt/n;
theo_kurt = p*(p + 2);
test_stat = (samp_kurt - theo_kurt)/sqrt(8*p*(p + 2)/n);
pval = 2*(1 - normcdf(abs(test_stat)));
kurt = [theo_kurt;samp_kurt;test_stat;Zval;pval];
% display result
if not(silent)
    disp(lin)
    disp(sprintf('KURTOSIS\nH0: Data ~Np(mu,sigma)\nTheoretical = %0.2f\nSample = %0.2f\nTest Stat = %0.2f\n%0.1f%% Accept Range = [%0.2f,%0.2f]\np-value = %0.5f',theo_kurt,samp_kurt,test_stat,100*(1-alpha),-Zval,Zval,pval));
    if abs(test_stat) > Zval
        disp('Reject H0: Data is not Gaussian')
    else
        disp('Can''t Reject H0: Data could be Gaussian')
    end
end

% test for skewness
for d1cnt = 1:n
    mcd1 = mc(d1cnt,:);
    for d2cnt = 1:n
        samp_skew = samp_skew + (mcd1*Sinv*mc(d2cnt,:)').^3;
    end             % data loop
end                 % data loop
samp_skew = samp_skew/(n^2);
theo_skew = 0;
test_stat = (n/6)*samp_skew;
pval = 1 - chi2cdf(test_stat,dof);
skew = [theo_skew;samp_skew;test_stat;chi2val;pval];

% display result and do nice accept / reject plots
if not(silent)
    disp(lin)
    disp(sprintf('SKEWNESS\nH0: Data ~Np(mu,sigma)\nTheoretical = %0.2f\nSample = %0.2f\nTest Stat = %0.2f\n%0.1f%% Accept Range = [0,%0.2f]\np-value = %0.5f',theo_skew,samp_skew,test_stat,100*(1-alpha),chi2val,pval));
    if abs(test_stat) > chi2val
        disp('Reject H0: Data is not Gaussian')
    else
        disp('Can''t Reject H0: Data could be Gaussian')
    end
    disp(lin)
    % KURT
    FH = figure('units','inches'); loc = get(FH,'position');
    x = sort([linspace(-4,4,100),kurt(3)]);     % make sure test stat eval'd, if too big/small
    subplot(1,2,1),plot(x,normpdf(x),'b-'),title('Test for Kurtosis')
    % draw lines for accept region
    hold on
    plot(repmat(-Zval,10,1),linspace(0,normpdf(-Zval),10),'r--',...
        repmat(Zval,10,1),linspace(0,normpdf(-Zval),10),'r--')
    % now the test stat
    plot(repmat(kurt(3),10,1),linspace(0,normpdf(kurt(3)),10),'g--')
    %text(kurt(3),normpdf(kurt(3)),{sprintf('*Test Stat: %0.2f',kurt(3));sprintf('p-value: %0.2f',kurt(5))},'Rotation',0);
	% JAH 20120305 change text multiline using sprintf
	text(kurt(3),normpdf(kurt(3)),sprintf('*Test Stat: %0.2f\np-value: %0.2f',kurt([3,5])),'Rotation',0);    
	hold off
    % SKEW
    x = sort([linspace(0,chi2inv(0.99,dof),100),skew(3)]);     % make sure test stat eval'd, if too big/small
    subplot(1,2,2),plot(x,chi2pdf(x,dof),'b-'),title('Test for Skewness')
    % draw lines for accept region
    hold on
    plot(repmat(chi2val,10,1),linspace(0,chi2pdf(chi2val,dof),10),'r--')
    % now the test stat
    plot(repmat(skew(3),10,1),linspace(0,chi2pdf(skew(3),dof),10),'g--')
    %text(skew(3),chi2pdf(skew(3),dof),{sprintf('*Test Stat: %0.2f',skew(3));sprintf('p-value: %0.2f',skew(5))},'Rotation',0);
	% same thing here JAH 20120305
    text(skew(3),chi2pdf(skew(3),dof),sprintf('*Test Stat: %0.2f\np-value: %0.2f',skew([3,5])),'Rotation',0);
    hold off
    set(FH,'position',[loc([1,2]),6,2.5])
end

%{
JAH 20070516, adapted for octave 3.4.3 20120305

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
