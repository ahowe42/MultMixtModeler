function meanm = MixMeanInit(data,ks)
% initialized means = MixMeanInit(data matrix, number groups)
%  Determine the initial estimated mean vectors for data drawn from a
%  mixture distributions using Dr. Hamparsum Bozdogan's mean initialization
%  technique.
%
%  Where
%  data matrix --- (nxp) matrix of data to analyze
%  number groups --- (scalar) number of groups or classes (k) supposed in data
%  initialized means --- (kxp) matrix of initialized means, with the mean
%     vector for group i in row i.

if (nargin ~= 2)
    % wrong number arguments
    fprintf('MixMeanInit: INVALID USAGE-Please read the following instructions!\n'), help MixMeanInit, return
end

minvec = min(data);     % minimum values
maxvec = max(data);     % maximum values
p = size(data,2);       % number dimensions
meanm = zeros(p,ks);
kms = ks - 1;

for ns = 1:1:p
   xbm = maxvec(ns)*eye(ks);
   xbm(1,1) = 0.5*(minvec(ns) + xbm(1,1));
   xbm(2:ks,1) = minvec(ns)*ones(kms,1);
   for is = 2:1:ks
      ism = is - 1;
      xbm(is,1) = 0.5*(xbm(is,1) + xbm(ism,1));
      for js = 2:1:ism
         xbm(is,js) = 0.5*(xbm(ism,(js - 1)) + xbm(ism,js));
      end
      xbm(is,is) = 0.5*(xbm(ism,ism) + xbm(is,is));
   end
   meanm(ns,:) = xbm(ks,:);
end
meanm = meanm';

% Written By:  Ike Patterson On:  25 May 96, updated JAH 20071008, adapted for octave 3.4.3 20120310
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
