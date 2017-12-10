function tabout = DispMeanCovar(mu, sigma, sizestr)
% string table = DispMeanCovar(mean vector, covariance matrix, size string)
%  This function is simple - take in a mean vector and covariance matrix,
%  and display both in a nice formatted string.
%
%  Where
%  mean vector  --- (1xp) mean vector
%  covariance matrix --- (pxp) covariance matrix
%  size string --- optional character string like %0.2f used for formatting
%     all values, if omitted, I will figure out the proper string
%  string table --- string holding nicely formatted parameters
%
%  Example: disp(DispMeanCovar([-3,9,-7],[2,0.5,-0.5;0.5,1,0.25;-0.5,0.25,0.5]))
%
% See Also: StrPad, table2str, ICSubTable, MatrixtoStr, MakeLaTeXTable.

p = length(mu);
ps = size(sigma);
if ((nargin ~= 2) && (nargin ~= 3)) || (ps(1) ~= ps(2)) || (p ~= ps(1))
    % wrong number args, covar not square, size mismatch between mean and covar
    fprintf('DispMeanCovar: INVALID USAGE-Please read the following instructions!\n'), help DispMeanCovar, return
end

if nargin ~= 3
    % get the proper size string for sprintf
    allvec = unique([mu(:);sigma(:)]);          % vector of unique values
    allvecS = sprintf('%0.0f',max(allvec));     % string the largest
    % now get the length of the largest and add 1 if any negatives present
    ndigits = length(allvecS) + 1*(sign(min(allvec)) == -1);
    % now get the longest decimal
    ndecs = 0;
    for valcnt = 1:length(allvec)
       ps = strfind(num2str(allvec(valcnt)),'.');
       if not(isempty(ps)) && ps > ndecs
           ndecs = length(num2str(allvec(valcnt))) - ps;
       end
    end
    sizestr = sprintf('%%0%0.0f.%0.0ff',ndigits+ndecs,ndecs);
end

lenp = length(sprintf('%0.0f',p));
% make a table of both params - leave room on front for 1:p
tbcovr = table2str(cellstr(num2str([1:p]')),[sigma;mu],{sizestr},2 + lenp);
% separate the mean vector and pipe it
tbmean = tbcovr((p + 4),:); tbmean([3,end]) = '|';
% separate the covar matrix
tbcovr = [tbcovr([2:(p + 3)],:);tbcovr(end,:)];
% put a box around the covar matrix
tbcovr(2,[3,4]) = '-'; tbcovr(end,[3,4]) = '-'; tbcovr([3:(end-1)],3) = '|';
% insert the 1:p on the front
tbcovr([3:(end-1)],[1:lenp]) = num2str([1:p]');
% pipe the covar matrix
tbcovr([3:(end-1)],end) = '|';
% finish the box for the mean vector
tbmean = [tbcovr(end,:);tbmean;tbcovr(end,:)];
% finalize - put it all together
maxlen = size(tbmean,2);
if (maxlen < 11 ) || (maxlen < 17)
    % :-( must resize tables to fit titles
    tabout = StrPad('Mean Vector',max([maxlen,17]),'R');
    for cnt = 1:3   % resize mean vector
        tabout = [tabout;StrPad(tbmean(cnt,:),17,'L')];
    end
    tabout = [tabout;'Covariance Matrix'];
    for cnt = 1:(p+3)   % resize mean vector
        tabout = [tabout;StrPad(tbcovr(cnt,:),17,'L')];
    end    
else
    tabout = [StrPad('Mean Vector',max([maxlen,17]),'R');tbmean;StrPad('Covariance Matrix',maxlen,'R');tbcovr];
end

% JAH 20070810, adapted for octave 3.4.3 20120315
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered