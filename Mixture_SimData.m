%{
  Usage: Mixture_SimData
  After Mixture_LoadData is run, if data_type = -1, this script can be used
  to simulate data from a known mixture of power exponentials.  Simulation
  parameter variables must be set up as in Mixture_LoadData; n (number
  observations) comes directly from Mixture_GUI.  If you want to simulate
  from a mixture of Gaussians, set betas = 1.

  Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}
Data = []; labels = [];
ns = round(MMM.n*mixpropors); n = sum(ns);  % NOTE: n comes from the GUI
for kcnt = 1:K
    Data = [Data;genrndmvpexp(ns(kcnt),p,meanvecs(:,:,kcnt)',covrmats(:,:,kcnt),betas(kcnt))];        
    labels = [labels;ones(ns(kcnt),1)*kcnt];
end
clear kcnt ns

%{
JAH 20070122, 20120224, adapted for octave 3.4.3

Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe

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