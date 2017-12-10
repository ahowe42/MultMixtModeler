% usage: Mixture_SimData
% After Mixture_LoadData is run, if data_type = -1, this script can be used
% to simulate data from a known mixture of power exponentials.  Simulation
% parameter variables must be set up as in Mixture_LoadData; n (number
% observations) comes directly from Mixture_GUI.  If you want to simulate
% from a mixture of Gaussians, set betas = 1.

Data = []; labels = [];
ns = round(MMM.n*mixpropors); n = sum(ns);  % NOTE: n comes from the GUI
for kcnt = 1:K
    Data = [Data;genrndmvpexp(ns(kcnt),p,meanvecs(:,:,kcnt)',covrmats(:,:,kcnt),betas(kcnt))];        
    labels = [labels;ones(ns(kcnt),1)*kcnt];
end
clear kcnt ns

% JAH 20070122, 20120224
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
