function tabout = MatrixtoStr(matrx, sizestr)
% string table = MatrixtoStr(data matrix, size string)
%  This function is simple - take any numeric data matrix, and display in 
%  a nice formatted string.
%
%  Where
%  data matrix  --- (nxp) data matrix
%  size string --- optional character string like %0.2f used for formatting
%     all values, if omitted, I will figure out the proper string
%  string table --- string holding nicely formatted data table
%
%  Example: disp(MatrixtoStr([2,0.5,-0.5;0.5,1,0.25;-0.5,0.25,0.5]))
%
% See Also table2str, ICSubTable, DispMeanCovar, StrPad, MakeLaTeXTable, 

if ((nargin ~= 1) && (nargin ~= 2))
    % wrong number args
    fprintf('MatrixtoStr: INVALID USAGE-Please read the following instructions!\n'), help MatrixtoStr, return
end

[n,p] = size(matrx);

if nargin ~= 2
    if islogical(matrx)
        % "sign" is undefined for logicals, and besides, the correct
        % sizestr is obvious for logicals.
        sizestr = '%1.0f';
    else
        % get the proper size string for sprintf
        allvec = unique(matrx(:));                  % vector of unique values
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
end

lenp = length(sprintf('%0.0f',p)); lenn = length(sprintf('%0.0f',n));
% make a table of data - leave room on front for 1:n
tabout = table2str(cellstr(num2str([1:p]')),matrx,{sizestr},2 + lenn);
% drop first line
tabout = tabout([2:end],:);
% put a box around the matrix
tabout([2,end],[3,4]) = '-'; tabout([3:(end-1)],[3,end]) = '|';
% insert the 1:n on the front
tabout([3:(end-1)],[1:lenn]) = num2str([1:n]');

% JAH 20070824, adapted for octave 3.4.3 20120315
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
