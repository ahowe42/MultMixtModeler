% usage: Mixture_LoadData
% This will load the data for M3 in data_file.  If simulated data is 
% selected, this will load the simulation parameters, but will not do the 
% simulation.  There are 3 possible values for data_type in M3: -1=simulated, 
% 0=known with labels, 1 = unknown (no labels).  Required input file formats:
% simulated:
%   row 1:        k
%   row 2:        tab separated mean vector
%   row 3-p+2:    tab separated covariance matrix
%   row p+3:      beta tab mixing proportion
%   must have (k-1)*(p+2) additional rows in the same format
% known:
%   [num_obs, num_dims+1]-sized, tab-delimited file
%   with labels (1,2,...) in FIRST column
% unknown:
%   [num_obs, num_dims]-sized, tab-delimited file

% data_mixnorm_?: 1=mixed/overlap,2=spherical/nonoverlap,3=ellipsoidal/overlap,4=spherical/overlap
% JamesData_?: 1=ellipsoidal/nonoverlap,2=ellipsoidal/overlap,3=mixed/overlap

datinput = dlmread([MMM.data_path,filesep,MMM.data_file]);

switch MMM.data_type
    case -1     % simulated data must be regenerated for each run
        K = datinput(1,1); p = size(datinput,2);
        meanvecs = zeros(1,p,K); covrmats = zeros(p,p,K);
        mixpropors = zeros(1,K); betas = mixpropors;

        betarows = [2:(p+2):(K*(p+2)+1)]; meanrows = [3:(p+2):(K*(p+2)+1)];
        covrrows = [4:(p+2):(K*(p+2)+1)]; covrrows = [covrrows;covrrows+(p-1)];
        for kcnt = 1:K
            meanvecs(:,:,kcnt) = datinput(meanrows(kcnt),:);
            covrmats(:,:,kcnt) = datinput([covrrows(1,kcnt):covrrows(2,kcnt)],:);
            betas(kcnt) = datinput(betarows(kcnt),1);
            mixpropors(kcnt) = datinput(betarows(kcnt),2);
        end
        mixpropors = mixpropors/sum(mixpropors);
    case 0      % data of known structure (with labels)
        labels = datinput(:,1); Data = datinput(:,[2:end]);
        K = length(unique(labels)); [n,p] = size(Data);
        % "defragment" data - sort so clusters contiguous
        [labels,ind] = sort(labels); Data = Data(ind,:);
        tab = tabulate(labels); mixpropors = tab(tab(:,2) ~= 0,3)/100;
    case 1      % data of unknown structure
        Data = datinput; [n,p] = size(Data);
end
clear kcnt *rows datinput ind tab

% JAH 20070109, 20120224
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
