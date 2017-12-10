function outstr = StrPad(inpstr,len,LRS,padchar)
%{
  padded string = StrPad(input string, desired length, pad type, pad char)

  Take a string and pad it with spaces (or a desired character.  Padding
  can be done on the left, right, or split so the input string is
  centered.  If the length is already >= the desired length, nothing is
  done.

  Where
  input string --- string (no matrices) to be padded
  desired length --- maximum number of characters required in the output string
  pad type --- L:left, R:right, S:split
  pad char --- (optional) character to be used for padding
  padded string --- input string with padding added (if appropriate)

  Example: a = StrPad('Andrew Howe',15,'S','*')

 See Also table2str, ICSubTable, DispMeanCovar, MatrixtoStr, MakeLaTeXTable.
 
 Copyright (C) 2006 J. Andrew Howe
%}

LRS = upper(LRS); [row,col] = size(inpstr);
if (len < 1) || ((sum(LRS == 'LRS') == 0)) || (not(ischar(inpstr))) || ((nargin == 4) && (sum(size(padchar) == [1,1]) ~= 2))
    % invalid length, invalid padding location, input not string, or invalid pad char instruction
    fprintf('StrPad: INVALID USAGE-Please read the following instructions!\n'), help StrPad, return
end

% no padding character
if nargin ~= 4; padchar = ' '; end;

if col >= len    % if already longer than (or as long as) len, just return it
    outstr = inpstr; return;
else
    outstr = [];
end

for cntr = 1:row
    ln = length(inpstr(cntr,:));
    if LRS == 'L'        % left pad
        outstr = [outstr;[repmat(padchar,1,len-ln),inpstr(cntr,:)]];
    elseif LRS == 'R'    % right pad
        outstr = [outstr;[inpstr(cntr,:),repmat(padchar,1,len-ln)]];
    elseif LRS == 'S'    % split padding evenly
        d = len-ln;
        if mod(d,2) == 0    % can split 50-50
            outstr = [outstr;[repmat(padchar,1,d/2),inpstr(cntr,:),repmat(padchar,1,d/2)]];
        else                % can't split 50-50
            outstr = [outstr;[repmat(padchar,1,floor(d/2)),inpstr(cntr,:),repmat(padchar,1,ceil(d/2))]];
        end
    end    
end

%{
JAH 20060408, adapted for octave 3.4.3 20120305
 
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