%{
  Usage: genrndmixmvpexp
  This script will generate random numbers from a mixture of three bivariate
  power exponential distributions.  The returned data is stored in Data
  such that the first column is the group labels, and the other columns are
  the data.  Parameters used are stored in Mu1, Mu2, Mu3, Sigma1, Sigma2,
  Sigma3, Beta1, Beta2, Beta3, n1, n2, n3.  A formatted string with these
  parameters is stored in paratit.

  Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}
  
clear, close all, clc

Data = []; k = 3; p = 2;

% class 1
prompt = {'Mu 1:','Mu 2:','Sigma 1(>0):','Sigma 2(>0):','Sigma 3/4(>0):','Beta:','N:'};
def = {'0','0','2','5','0','2','85'};
answer = inputdlg(prompt,'Parameters for set #1',1,def,'off');
Mu1 = [str2num(char(answer(1,1))), str2num(char(answer(2,1)))];
Sigma1 = [str2num(char(answer(3,1))),str2num(char(answer(5,1)));...
    str2num(char(answer(5,1))),str2num(char(answer(4,1)))];
Beta1 = str2num(char(answer(6,1)));
n1 = str2num(char(answer(7,1)));
X = genrndmvpexp(n1,p,Mu1',Sigma1,Beta1); Data=[Data;[1*ones(n1,1),X]];

% class 2
prompt = {'Mu 1:','Mu 2:','Sigma 1(>0):','Sigma 2(>0):','Sigma 3/4(>0):','Beta:','N:'};
def = {'2','-2','1.5','3','0.5','5','85'};
answer = inputdlg(prompt,'Parameters for set #2',1,def,'off');
Mu2 = [str2num(char(answer(1,1))), str2num(char(answer(2,1)))];
Sigma2 = [str2num(char(answer(3,1))),str2num(char(answer(5,1)));...
    str2num(char(answer(5,1))),str2num(char(answer(4,1)))];
Beta2 = str2num(char(answer(6,1)));
n2 = str2num(char(answer(7,1)));
X = genrndmvpexp(n2,p,Mu2',Sigma2,Beta2); Data=[Data;[2*ones(n2,1),X]];

% class 3
prompt = {'Mu 1:','Mu 2:','Sigma 1(>0):','Sigma 2(>0):','Sigma 3/4(>0):','Beta:','N:'};
def = {'4','4','0.5','2','0.25','1','85'};
answer = inputdlg(prompt,'Parameters for set #3',1,def,'off');
Mu3 = [str2num(char(answer(1,1))) str2num(char(answer(2,1)))];
Sigma3 = [str2num(char(answer(3,1))),str2num(char(answer(5,1)));...
    str2num(char(answer(5,1))),str2num(char(answer(4,1)))];
Beta3 = str2num(char(answer(6,1)));
n3 = str2num(char(answer(7,1)));
X = genrndmvpexp(n3,p,Mu3',Sigma3,Beta3); Data=[Data;[3*ones(n3,1),X]];

% plot them
paratit = sprintf('Mu = [%1.1f,%1.1f; %1.1f,%1.1f; %1.1f,%1.1f]\nSigma = [%1.1f,%1.1f;%1.1f,%1.1f; %1.1f,%1.1f;%1.1f,%1.1f; %1.1f,%1.1f;%1.1f,%1.1f]\nBeta = [%1.1f; %1.1f; %1.1f]\nN = [%1.1f,%1.1f,%1.1f]',...
    Mu1,Mu2,Mu3,Sigma1,Sigma2,Sigma3,Beta1,Beta2,Beta3,n1,n2,n3); % title for plots
PlotBivarClusts(Data(:,[2,3]),Data(:,1));
title(sprintf('Mixture of %d Bivariate Power Exponential Random Variables\n%s',k, paratit)), xlabel('x1'), ylabel('x2')

% display
str = 'Class %1.0f Parms:Mu=[%1.2f,%1.2f],Sigma=[%1.2f,%1.2f; %1.2f,%1.2f],Beta=%1.2f,N=%1.1f';
disp(sprintf(str,1,Mu1,Sigma1,Beta1,n1))
disp(sprintf(str,2,Mu2,Sigma2,Beta2,n2))
disp(sprintf(str,3,Mu3,Sigma3,Beta3,n3))

%{
JAH 20061230, adapted for octave 3.4.3 20120312

Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe

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