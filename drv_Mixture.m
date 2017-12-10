% usage: THIS MUST BE CALLED FROM Mixture_noGUI
% This script will perform multiple simulations/replications of the core M3
% scripts - Mixture_GARM, Mixture_GKM, Mixture_KM.

% close all figures
close ALL;

clc, nw = clock;
warning('off','MATLAB:divideByZero'); warning('off','MATLAB:dispatcher:InexactMatch');

% make the data - range specific output directory, and the MULTI-RUN prefix
df = MMM.data_file([1:strfind(MMM.data_file,'.')-1]);
dir = [MMM.PStype,filesep,df,sprintf('_%0.0f',MMM.Kmax)];
if exist([mydir,filesep,'output',filesep,dir],'dir') ~= 7
    mkdir([mydir,filesep,'output'],dir);
end
dir = [filesep,'output',filesep,dir,filesep];

if isempty(GAmat)
    % BEGIN A NEW GARM/GKM MULTI-RUN
    GAscores = []; GAfil = []; GAchroms = [];
    GAmat = [sprintf(['MIXTMR_',df,'_%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f.mat'],nw)];
    GAmat = [dir,GAmat]; save([mydir,GAmat]);
else
    % AUGMENT A GARM/GKM MULTI-RUN
    GAmax = MMM.GAmax;
    load(GAmat);
    MMM.GAmax = GAmax;
    clear GAmax tmp; save([mydir,GAmat]);
end

Mixture_LoadData
while GAcnt <= MMM.GAmax
    if MMM.data_type == -1; Mixture_SimData; end;
    % update progress
    clc;
    MMM.lblCurrRun = sprintf('%1d of %1d',GAcnt,MMM.GAmax);
    Mixture_DispStatus
    % do the work
    eval(['Mixture_',MMM.init_type]);
    clear rnd_stat SCORES_CERRS f t stt tottim best_clustassigns
    % close all figures
    if MMM.GAmax > 1; close ALL; end;
    % don't want to overwrite the status variables
    tmp.lblCurrRun = MMM.lblCurrRun; tmp.lblCurrIC = MMM.lblCurrIC;
    tmp.lblCurrk = MMM.lblCurrk; tmp.lblTimRep = MMM.lblTimRep;
    tmp.lblTimK = MMM.lblTimK; tmp.lblTimInit = MMM.lblTimInit;
    tmp.lblTimIC = MMM.lblTimIC; tmp.lblTimTot = MMM.lblTimTot;
    % save the stuff
    load([mydir,GAmat])
    MMM.lblCurrRun = tmp.lblCurrRun; MMM.lblCurrIC = tmp.lblCurrIC;
    MMM.lblCurrk = tmp.lblCurrk; MMM.lblTimRep = tmp.lblTimRep;
    MMM.lblTimK = tmp.lblTimK; MMM.lblTimInit = tmp.lblTimInit;
    MMM.lblTimIC = tmp.lblTimIC; MMM.lblTimTot = tmp.lblTimTot;    
    GAchroms = [GAchroms;[ones(size(best_chroms,1),1)*GAcnt,best_chroms]];
    GAscores = [GAscores;[ones(size(best_scores,1),1)*GAcnt,best_scores]];
    GAfil = [GAfil;StrPad(savename([(length(mydir)+1):length(savename)]),100,'R',' ')];
    clear best* savename tmp
    GAcnt = GAcnt + 1; iccnt = 1; save([mydir,GAmat]);
end
warning('on','MATLAB:divideByZero'); warning('on','MATLAB:dispatcher:InexactMatch');
Mixture_MRSummary

% JAH 20061219, adapted for octave 3.4.3 20120324
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
