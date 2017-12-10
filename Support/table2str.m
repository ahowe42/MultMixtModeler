function tablestr = table2str(colheads,data,format,colpad,rowheads)
%{
  string table = table2str(col heads, data, formats, col pad, row heads)

  This function takes in a table of numbers, and returns a nicely
  formatted character string table that can be displayed with disp().

  Where
  col heads --- (px1) cell array of column headings
  data --- (nxp) matrix of data
  formats --- (px1) cell array of sprintf style column formats; if cell{1}
   is input, every column is formatted identically.
  col pad --- (1x1) or (1x2) number of blank columns on the far left (+) 
     and/or far right (-), can pass an element for each side.
  row heads --- (nx1) optional cell array of row headings
  string table --- string matrix of formatted table

  See Also StrPad, ICSubTable, DispMeanCovar, MatrixtoStr, MakeLaTeXTable.
  
  Copyright (C) 2005 J. Andrew Howe; see below
%}

[datr,datc] = size(data);
if (length(colheads) ~= datc) || ((nargin ~= 4) && (nargin ~= 5)) || ((nargin == 5) && (length(rowheads) ~= datr))
    % col head count should equal data column count, should be 4 or 5 inputs, row head count should equal data row count
    fprintf('table2str: INVALID USAGE-Please read the following instructions!\n'), help table2str, return
end

% if format just 1, replicate
if (sum(size(format) == 1) == 2)
    format = repmat(format,datc,1);
end

lens = zeros(datr+1,datc);
% loop through and get the length(colheads(i))
for cnt = 1:datc
    lens(1,cnt) = 1+size(colheads{cnt},2);
end
% loop through and get the length(num2str(data(i,j)))
for rcnt = 1:datr
    for ccnt = 1:datc
        lens(1+rcnt,ccnt) = 1+size(num2str(data(rcnt,ccnt),format{ccnt}),2);
    end    
end
maxes = max(lens); bars = repmat('-',[1,sum(maxes)]);
% make output table
tablestr = '';
for cnt = 1:datc
    extra = maxes(cnt)-size(colheads{cnt},2);
    if mod(extra,2) == 0    % even number - divide extra space on either side of col hEAD equally        
        tablestr = [tablestr, repmat(' ',[1,extra/2]), colheads{cnt}, repmat(' ',[1,extra/2])];
    else    % odd number - divide extra space on either side of data with extra on right
        tablestr = [tablestr, repmat(' ',[1,floor(extra/2)]), colheads{cnt}, repmat(' ',[1,ceil(extra/2)])];
    end
end
tablestr = [bars; tablestr; bars];
for rcnt = 1:datr
    row = '';
    for ccnt = 1:datc
        row = [row,StrPad(num2str(data(rcnt,ccnt),format{ccnt}),maxes(ccnt),'S')];
    end
    tablestr = [tablestr; row];
end
tablestr = [tablestr; bars];

% add the row headings, if input
if nargin == 5
    % get the length of the longest header
    maxes = 0;
    for cnt = 1:datr
        if length(rowheads{cnt}) > maxes; maxes = length(rowheads{cnt}); end;
    end
    maxes = maxes + 1;
    bars = repmat('-',1,maxes);
    rhds = [bars;blanks(maxes);bars];
    for cnt = 1:datr
        rhds = [rhds;StrPad(rowheads{cnt},maxes,'R')];
    end
    tablestr = [[rhds;bars],tablestr];
end

% colpad everything if requested
if not(isequal(colpad,0))
    for icnt = 1:length(colpad)
        if colpad(icnt) > 0
            tablestr = [repmat(' ',[datr+4,colpad(icnt)]), tablestr];
        elseif colpad(icnt) < 0
            tablestr = [tablestr,repmat(' ',[datr+4,-colpad(icnt)])];
        end        
    end             % colpad loop
end

%{
JAH 20051207, adapted for octave 3.4.3 20120305

Copyright (C) 2005 J. Andrew Howe

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