% usage: MultRegGaSub_noGUI
% This will perform multiple simulations/replications of the 
% multivariate mixture model, with multiple choices of model type.

%%% STILL NEED TO SETUP FOR THE DIFFERENT KINDS OF SUBSET ANALYSIS, AS WELL
%%% AS THE INFLUENCE ANALYSIS !!!
% Mixture_AllSubsAnal
% Mixture_GASubsAnal
% Mixture_OutDetect

clear

mydir = dbstack; mydir = mydir(end);
mydir = which(mydir.file);
tmp = strfind(mydir,filesep);
mydir = mydir([1:(tmp(end)-1)]);
% ensure output directory exists
if exist([mydir,filesep,'output'],'dir') ~= 7
    mkdir([mydir,filesep,'output']);
end

% SETUP THE DROPDOWN-TYPE CHOICE FIELDS
%MMM.regul_funcs = CovSmooth();          % covariance regularization
%MMM.init_types = {'GARM','GKM','KM'};   % mixture initialization
%MMM.em_type = {'GEM','EM'};             % mixture optimization

% FIRST SEE IF WE'RE AUGMENTING AN EXISTING RUN
GAmat = ''; GAcnt = 1;
[fn,pn] = uigetfile('*.mat', 'Please Select Experiment Results to Augment to Load (OPTIONAL).',[mydir,filesep,'output']);
if pn ~= 0                  % load the results & set relevant parameters
    load([pn,filesep,fn]);
    GAmat = [pn,fn];
    GAcnt = size(GAfil,1)+1;
    answer = inputdlg({'Total Number Replications'},'Multivariate Mixture Modeler',1,{int2str(GAcnt)},'off');
    if isequal(answer,{})
        % user clicked cancel, so let's cancel
        return;
    else
        MMM.GAmax = str2double(answer{1});
    end
else                        % setup everything from scratch
    % INITIALIZE MODEL PARAMETERS
    % first get the common model settings
    MMM = load('MixtureParams.mat');
    % set the Ks evaluated
    if MMM.JustK
        MMM.Knum = 1;                   % just Kmax
        MMM.Krange = MMM.Kmax;
    else
        MMM.Knum = MMM.Kmax;            % 1 to Kmax
        MMM.Krange = [1:MMM.Kmax];
    end
    % now, load any model-specific parameters the user may wish
    [fn,pn] = uigetfile('*.mat', 'Please Select a Set of Model-Specific Parameters to Load (OPTIONAL).',[mydir,filesep]);
    if pn == 0
        % user canceled, so just set to Gaussian
        tmp = load('GaussianParams.mat');
        % loop through all the field in this PS params file, and add them to the existing struct
        flds = fieldnames(tmp);
        for i = 1:numel(flds)
          MMM = setfield(MMM, flds{i}, tmp.(flds{i}));
        end
    else
        % need to add items in this parameter file into the MMM object
        tmp = load([pn,filesep,fn]);
        % loop through all the field in this PS params file, and add them to the existing struct
        flds = fieldnames(tmp);
        for i = 1:numel(flds)
          MMM = setfield(MMM, flds{i}, tmp.(flds{i}));
        end
        % need to copy ECtype to htype for EC
        if isequal(MMM.PStype,'EC')
            MMM.htype = MMM.ECtype;
        elseif isequal(MMM.PStype,'Gaussian') || isequal(MMM.PStype,'PEKern')
            MMM.htype = NaN;
        end;
    end

    % INITIALIZE THE DATA PARAMETERS
    % setup parms data_type: There are 3 possible values for data_type in M3: -1=simulated, 
    % 0=known with labels, 1 = unknown (no labels)
    [MMM.data_file,MMM.data_path] = uigetfile('*.m', 'Load Data File for Modeling',[mydir,filesep]);
    if MMM.data_path == 0
        % user clicked cancel, so let's cancel
        return;        
    else
        qans = questdlg('What Type of Data is This?','Multivariate Mixture Modeler','Simulated','Known','Unknown','Simulated');
        switch qans
            case 'Simulated'
                MMM.data_type = -1;
                % get the number observations to simulate  
                answer = inputdlg({'Number Simulated Observations'},'Multivariate Mixture Modeler',1,{int2str(150)},'off');
                if isequal(answer,{})
                    return;
                else
                    MMM.n = str2double(answer{1});
                end                
            case 'Known'
                MMM.data_type = 0;
            case 'Unknown'
                MMM.data_type = 1;
            otherwise
                return;
        end
    end
    % check if data_type seems to match what is in the file before doing anything
    datinput = dlmread([MMM.data_path,filesep,MMM.data_file]);    % load the datafile
    if (((size(find(datinput(1,:)),2) == 1) + (MMM.data_type == -1)) == 1)
        % assuming data is not univariate, first row should only have just 1 nonzero
        % element IF data is simulated, otherwise, it's probably known or unknown
        qans = questdlg('Please Check - is the Data Structure Correct?','Potential Data Error');
        if not(isequal(qans,'Yes')); return; end;
    end

    % clear unnecessary shit
    clear datinput flds fn i pn qans tmp;
end    

% FINALLY RUN THE MODEL
% initialize modeling progress strings
MMM.lblCurrRun = 'NONE';
MMM.lblCurrIC = 'NONE';
MMM.lblCurrk= 'NONE';
MMM.lblTimRep = 'Replic: 0.0';
MMM.lblTimK = 'k: 0.0';
MMM.lblTimInit = [MMM.init_type,': 0.0'];
MMM.lblTimIC = 'IC: 0.0';
MMM.lblTimTot = 'Total (m): 0.0';
% set the start time
MMM.totstt = clock;

% temporary debug stuff - set some parameters so it'll run faster
%MMM.Knum = 3; MMM.Krange = [1,2,3]; MMM.num_generns = 10; MMM.popul_size = 10; MMM.GAmax = 3;
% temporary debug stuff - control the tests
%return
%MMM.init_type = 'GARM'; MMM.optim_type = 'GEM'; MMM.optim_func = 'KMixPEPK_GEM'; drv_Mixture
%MMM.init_type = 'GARM'; MMM.optim_type = 'EM'; MMM.optim_func = 'KMixPEPK_EM'; drv_Mixture
%MMM.init_type = 'GKM'; MMM.optim_type = 'GEM'; MMM.optim_func = 'KMixPEPK_GEM'; drv_Mixture
%MMM.init_type = 'GKM'; MMM.optim_type = 'EM'; MMM.optim_func = 'KMixPEPK_EM'; drv_Mixture
%MMM.init_type = 'KM'; MMM.optim_type = 'EM'; MMM.optim_func = 'KMixPEPK_EM'; drv_Mixture
% let's do this thing!
drv_Mixture