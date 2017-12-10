function X = genrndmvpexp(n,p,Mu,Sigma,Beta,plotflg)
%{
  rands = genrndmvpexp(n,p,Mu,Sigma,Beta,scatter plot)

  Returns a sample matrix from multivariate Power Exponential distribution
  with a given mean vector, covariance matrix, and kurtosis parameter. For
  repeatability, the states for both randn and rand must be set.

  Where
  n --- scalar number of sample sets to be generated
  p --- scalar number of variables
  Mu --- (px1) mean vector
  Sigma --- (pxp) covariance matrix
  Beta --- scalar kurtosis parameter of PE distribution
  scatter plot --- 1 = make scatter plot, 0 (default) = don't
  rands --- (nxp) matrix of random data

  Example: pexp = genrndmvpexp(500,2,[0;0],[1,0;0,1],2,1);

  See Also GENRNDMIXMVPEXP.
  
  Copyright (C) 2006 J. Andrew Howe; see below
%}

if (nargin < 5) || (nargin > 6) || (not(isscalar(p))) || (not(isscalar(n))) || (not(isscalar(Beta))) ...
        || (p < 1) || (n < 1) || (sum(size(Mu) == [p,1]) ~= 2) || ((sum(size(Sigma) == [p,p]) ~= 2))
    % not 5/6 args, nonscalar p/n/beta, nonpositive p or n, missized mu or sigma
    help genrndmvpexp, return
end

if nargin == 5; plotflg = 0; end;

% Generate n-by-p random points uniformly distributed on the surface of a
% hypersphere of p dimension
a = randn(n,p);                 %a = normrnd(0,1,n,p); changed JAH 20080204
%a_t = a';                      % comment out JAH 20080204
a1 = sqrt(sum(a.^2,2));         %a1 = sum(a_t.^2)'; changed JAH 20080204
a2 = a1*ones(1,p);              %a2 = sqrt(a1*ones(1,p)); changed JAH 20080204
u = a./a2;

% Generate Gamma(1/2,p/2beta)
g = gamrnd(p/2/Beta,2,n,1);
R = g.^(1/2/Beta);

% Cholesky Factorization of SIGMA
A = chol(Sigma);

% Generate the PE matrix
X = (Mu*ones(1,n) + (ones(p,1)*R').*(A'*u'))';

if plotflg; scatter(X(:,1),X(:,2),[],[],'.'),xlabel('x1'),ylabel('x2'),title('Bivariate Power Exponential'); end;

%{
JAH 20061230, 20080204, 20120218

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
