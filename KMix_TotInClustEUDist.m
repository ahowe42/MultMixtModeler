function total_var = KMix_TotInClustEUDist(data, labels)
% total distance = KMix_TotInClustEUDist(data matrix, mixture assignments)
%  This will compute the total within cluster Euclidian variance for a mixture model.
%
%  Where
%  data matrix --- (nxp) data matrix
%  mixture assignments --- n-valued vector of assignments 1:k
%  total distance --- total within cluster variance

[n,p] = size(data);

if (length(labels) ~= n) || (nargin ~= 2)
    % mismatched dimensions, wrong # of arguments
    fprintf('KMix_TotInClustEUDist: INVALID USAGE-Please read the following instructions!\n'), help KMix_TotInClustEUDist, return
end

ks = unique(labels); kcnt = length(ks); % must do this in case some k missing (shouldn't)
total_var = zeros(1,kcnt); means = zeros(1,p,kcnt);

% compute the cluster means
for mixcnt = 1:kcnt
    means(:,:,mixcnt) = mean(data(labels == ks(mixcnt),:),1);
end     % mixtures loop

% get the total variance for each cluster
for datcnt = 1:n
    thisclst = find(labels(datcnt) == ks);
    total_var(thisclst) = total_var(thisclst) + sum((data(datcnt,:) - means(:,:,thisclst)).^2);
end     % data loop
total_var = sum(total_var);

% JAH 20061215, adapted for octave 3.4.3 20120324
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
