%{
  Usage: THIS MUST BE CALLED FROM drv_Mixture
  This script will perform one run of mixture modeling using the K-means
  algorithm to initialize the EM algorithm.  NOTE: this can only be used
  with the EM algorithm

  Copyright (C) 2007 Prof. Hamparsum Bozdogan & J. Andrew Howe
%}

drvstt = clock; rnd_stat = sum(drvstt*1000000); rand('state',rnd_stat);
stt = sprintf('%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f',drvstt);
savename = [mydir,dir,MMM.init_type,'+',MMM.optim_type,'_',stt]; diary([savename,'.out']);
disp(repmat('#',1,50)), Mixture_DispParms
disp(sprintf('Random State: %0.0f',rnd_stat))
disp(['Save files: ',savename([(length(mydir)+1):length(savename)])]), disp(repmat('#',1,50))

maxP2disp = 8;      % if p > this, don't display estimated mu, sigma, h
% SCORES_CERRS = [score,%_correct,actual_k,attempt_k]
SCORES_CERRS = ones(MMM.Knum,4,10)*Inf;
SCORES_CERRS(:,2,:) = NaN;                  % initialize %_correct for data_type = 1
SCORES_CERRS(:,4,:) = repmat(MMM.Krange',1,10); % and attempt_k
best_clustassigns = zeros(n,MMM.Knum,10); global ys; ys = NaN;  % for ICOMP_PERF with k =1;

for kcnt = 1:MMM.Knum
    kstt = clock;
    MMM.lblCurrk = sprintf('%0.0f of %0.0f',MMM.Krange([kcnt,MMM.Knum]));
    MMM.lblCurrIC = MMM.init_type;
    Mixture_DispStatus
    best_score = Inf;
    if MMM.Krange(kcnt) == 1    % do 1 cluster
        if (MMM.data_type ~= 1)
            [val,ind] = max(mixpropors);
            best_clustassigns(:,1,:) = ones(1,n,10)*ind;
        else
            best_clustassigns(:,1,:) = ones(1,n,10);
        end
        mix_propors = 1;
        if isequal(MMM.PStype,'EC')
            mix_meanvecs = zeros(1,p,1); mix_covrmats = zeros(p,p,1);
            mix_ECparms = zeros(1,2+isequal(MMM.htype,'KT'),1);
            [mix_meanvecs(:,:,1),mix_covrmats(:,:,1),mix_ECparms(:,:,1)]= MECEstimate(Data,MMM.htype,MMM.regul_func,0);
        else
            mix_meanvecs(:,:,1) = mean(Data,1);
            mix_covrmats(:,:,1) = CovSmooth(Data,MMM.regul_func,1,1,n);
        end
        % display final parameters
        disp(sprintf('Cluster %0.0f',1))
        disp(['pi = ',num2str(100*mix_propors(1),'%0.2f%%')])
        if p <= maxP2disp
            disp(DispMeanCovar(mix_meanvecs(:,:,1),mix_covrmats(:,:,1),'%0.3f'))
        end
        switch MMM.PStype
            case 'Gaussian' % nothing special to do
            case 'Kernel'
                if p <= maxP2disp
                    [tmp,H] = MVKDE_Gauss(Data,zeros(1,p),MMM.htype,MMM.regul_func,0);                
                    disp('Kernel Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'))
                end
            case 'EC'
                if isequal(MMM.htype,'KT')
                    disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',MMM.htype,mix_ECparms(:,:,1)))
                else
                    disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',MMM.htype,mix_ECparms(:,:,1)))
                end
            case 'PEKern'
                [tmp,H,Bs] = MPEProdKernelPDF(Data,zeros(1,p),MMM.regul_func);
                disp('PE Kernel Bandwidths'),disp(mat2str(H,2))
                disp('PE Kernel Betas'),disp(mat2str(Bs,2))
        end        
        % compute Percent Correct & objective value, do confusion matrix
        SCORES_CERRS(1,3,:) = 1;
        for infcnt = 1:10
            if not(isempty(MMM.InfCrit{infcnt}))
                switch MMM.PStype
                    case {'Gaussian','PEKern'}
                        SCORES_CERRS(1,1,infcnt) = feval(MMM.InfCrit{infcnt},...
                            Data,best_clustassigns(:,1,infcnt),MMM.regul_func,0);
                    case 'Kernel'
                        SCORES_CERRS(1,1,infcnt) = feval(MMM.InfCrit{infcnt},...
                            Data,best_clustassigns(:,1,infcnt),MMM.regul_func,MMM.htype,0);
                    case 'EC'
                        SCORES_CERRS(1,1,infcnt) = feval(MMM.InfCrit{infcnt},...
                            Data,best_clustassigns(:,1,infcnt),MMM.regul_func,MMM.htype,0);
                end
            end
        end         % information criteria loop
        if MMM.data_type ~= 1   % not unknown
            SCORES_CERRS(1,2,:) = 100*sum(best_clustassigns(:,1,1) == labels)/n;
            disp(sprintf('Percent Correct = %0.2f%%\n',SCORES_CERRS(1,2,1)))
        end
        continue;
    end
    
    k = MMM.Krange(kcnt); tic
    % Begin K-means Algorithm (KM)
    disp(sprintf('KM Algorithm Entered for k = %0.0f',k))

    % use Bozdogan's initialization scheme
    initmeanvecs = MixMeanInit(Data,k);
    initmeanvecs = reshape(initmeanvecs',1,p,k);

    % initialize more "things"
    mix_covrmats = zeros(p,p,k); mix_propors = zeros(1,k); itercnt = 0;
    mix_meanvecs = initmeanvecs(:,:,[1:MMM.Krange(kcnt)]); totdists = [];
    mix_ECparms = zeros(1,2+isequal(MMM.htype,'KT'),k);
    population = zeros(1,n);
    
    while(itercnt <= MMM.maxiter)
        % ASSIGN DATAPOINTS TO CLOSEST MIXTURE
        for datcnt = 1:n
            dist_from_clusts = zeros(k,1);
            % compute euclidian distances from each centroid
            for mixcnt = 1:k
                dist_from_clusts(mixcnt) = sum((Data(datcnt,:) - mix_meanvecs(:,:,mixcnt)).^2);
            end     % mixtures loop
            % assign datapoint to nearest cluster
            [val,ind] = min(dist_from_clusts);
            population(datcnt) = ind;
        end         % datapoints loop
        
        % RECOMPUTE THE CENTROIDS & NEW TOTAL WITHIN-CLUSTER DISTANCE
        for mixcnt = 1:k
            ind = (population == mixcnt);
            if sum(ind) ~= 0
                mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
            end
        end         % mixtures loop
        totdists = [totdists,KMix_TotInClustEUDist(Data,population)];
        
        % DO WE HAVE CONVERGENCE?
        if (itercnt > 5) && (abs(totdists(end) - totdists([end-1])) <= MMM.convgcrit)
            % total within-cluster distance converged        
            disp(table2str({'Iteration';'Total Disance'},[[0:itercnt]',totdists'],{'%0.0f';'%0.3f'},0))
            disp(sprintf('Total Within-Cluster Distance converged in iteration %0.0f',itercnt))
            break;
        end
        itercnt = itercnt + 1;
    end
    
    % since it's ok for KM to drop clusters, let's make sure the EM
    % algorithm gets labels from 1:some k
    if (sum(unique(population)) ~= sum([1:k]))
        tab = tabulate(population); tab = tab(tab(:,2) ~= 0,:);
        for mixcnt = 1:k
            ind = find(tab(:,1) == mixcnt);
            if not(isempty(ind))
                population([population == mixcnt]) = ind;
            end
        end         % mixtures loop
    end
    
    ys = population; % here for ICOMP_PERF functions

    KM_toc = toc;
    disp(sprintf('KM Completed in \n\t%1.0f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',KM_toc./[1,60,3600]));
    MMM.lblTimInit = [MMM.init_type,sprintf(': %0.2f',KM_toc/60)];
    MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
    Mixture_DispStatus
    % End K-means Algorithm (KM)

    % Begin Expectation-Maximization EM
    tic;
    objec_func = 'Log-Likelihood';
    MMM.lblCurrIC = objec_func;
    Mixture_DispStatus
    disp(sprintf([MMM.optim_type,' Algorithm Entered for k = %0.0f: ',objec_func],k))

    % actually do the optimization    
    if isequal(MMM.PStype, 'EC') || isequal(MMM.PStype, 'Kernel')
        [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
            Data,population,[MMM.convgcrit,MMM.maxiter],MMM.htype,MMM.regul_func,...
            [savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
    else
        [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
            Data,population,[MMM.convgcrit,MMM.maxiter],MMM.regul_func,...
            [savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
    end
    emconvg = not(isnan(mix_propors(1)));

    % properly align cluster labels
    if (emconvg == 1) && (MMM.data_type ~= 1)
        [val,bins] = max(posteriors,[],2);              % assign to clusters based on maximum posterior
        kselu = unique(bins); ksell = length(kselu);    % account for the possibility of missing k
        cm = ConfMatrix(labels,bins,1);                 % rows actual, columns estimated
        cm = cm([1:(end-1)],[1:(end-1)]); cm = cm([1:(K + 1)],:);
        % determine most frequent matches
        maxs = ones(3,k)*Inf;
        for mixcnt = 1:ksell
            if sum(bins == kselu(mixcnt)) ~= 0
                [val,ind] = max(cm([2:end],cm(1,:) == kselu(mixcnt)));
                maxs(:,mixcnt) = [ind;val;kselu(mixcnt)];
            end
        end         % mixtures loop
        % assign based on max matches
        for mixcnt = 1:K
            val = max(maxs(2,maxs(1,:) == mixcnt));
            if not(isempty(val))
                ind = maxs(3,and(maxs(1,:) == mixcnt, maxs(2,:) == val));
                bins(bins == ind(1)) = 100 + mixcnt;
            end
        end         % mixtures loop
        % see if we need to do more clean up
        for mixcnt = 1:K
            cm = (sum(bins == mixcnt) > 0);
            if cm; break; end;
        end         % mixtures loop
        if cm
            bins = bins - 100;
            % ensure 1:k accounted for (unless something missing, which is ok)
            for mixcnt = 2:k
                ind = (bins == mixcnt);
                if sum(ind) == 0
                    % find a number that's not supposed to be there and reassign it
                    tab = tabulate(bins); tab = tab(:,1); tab = tab(not(ISIN(tab,[1:k],1)));
                    if not(isempty(tab)); bins(bins == tab(1)) = mixcnt; end;
                end % mixtures loop
            end
        else
            % just take off the 100
            ind = find(bins > 100); bins(ind) = bins(ind) - 100;
        end
        best_clustassigns(:,kcnt,1) = bins;
    elseif  emconvg == 1
        [val,best_clustassigns(:,kcnt,1)] = max(posteriors,[],2);
    else
        best_clustassigns(:,kcnt,1) = NaN; % em did not converge
    end
    kselu = unique(best_clustassigns(:,kcnt,1)); ksell = length(kselu);    % account for the possibility of missing k

    % recompute (in case assignments changed) and display estimates
    disp(sprintf(['\nResults for k = %0.0f: ',objec_func],k))
    if (emconvg == 1) && (MMM.data_type ~= 1)    
        for mixcnt = 1:ksell
            ind = (best_clustassigns(:,kcnt,1) == kselu(mixcnt));
            if sum(ind) ~= 0
                mix_propors(mixcnt) = sum(ind)/n;                
                if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                    mix_covrmats(:,:,mixcnt) = 0;
                else
                    if isequal(MMM.PStype,'EC')
                        [mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),...
                            mix_ECparms(:,:,mixcnt)]= MECEstimate(Data(ind,:),...
                            MMM.htype,MMM.regul_func,0);                        
                    else
                        mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                        mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),MMM.regul_func,1,1,sum(ind));
                    end
                end
            else
                mix_propors(mixcnt) = 0;
                mix_meanvecs(:,:,mixcnt) = nan(1,p);
                mix_covrmats(:,:,mixcnt) = nan(p,p);
            end
        end         % mixtures loop
    end

    if (emconvg == 1)
        for mixcnt = 1:ksell
            disp(sprintf('Cluster %0.0f',kselu(mixcnt)))
            disp(['pi = ',num2str(100*mix_propors(mixcnt),'%0.2f%%')])
            if p <= maxP2disp
                disp(DispMeanCovar(mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'%0.3f'))
            end            
            switch MMM.PStype
                case 'Gaussian' % nothing special to do
                case 'Kernel'
                    if p <= maxP2disp
                        [tmp,H] = MVKDE_Gauss(Data([best_clustassigns(:,kcnt,infcnt) ...
                            == kselu(mixcnt)],:),zeros(1,p),MMM.htype,MMM.regul_func,0);
                        disp('Kernel Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'));
                    end
                case 'EC'
                    if isequal(MMM.htype,'KT')
                        disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',MMM.htype,mix_ECparms(:,:,mixcnt)))
                    else
                        disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',MMM.htype,mix_ECparms(:,:,mixcnt)))
                    end
                case 'PEKern'
                    [tmp,H,Bs] = MPEProdKernelPDF(Data([best_clustassigns(:,kcnt,infcnt)...
                        == kselu(mixcnt)],:),zeros(1,p),MMM.regul_func);
                    disp('PE Kernel Bandwidths'),disp(mat2str(H,2))
                    disp('PE Kernel Betas'),disp(mat2str(Bs,2))
            end
        end         % mixtures loop
        % compute Percent Correct & do confusion matrix
        SCORES_CERRS(kcnt,3,:) = ksell;
        if MMM.data_type ~= 1   % not unknown
            SCORES_CERRS(kcnt,2,:) = 100*sum(best_clustassigns(:,kcnt,1) == labels)/n;
            confmat = ConfMatrix(labels,best_clustassigns(:,kcnt,1),0); % JAH 20120227
            disp(sprintf('Percent Correct = %0.2f%%\n',SCORES_CERRS(kcnt,2,1)))
        end
    else
        disp([MMM.optim_func,' FAILED TO CONVERGE - NO ESTIMATES!'])
        SCORES_CERRS(kcnt,1,:) = Inf;
        SCORES_CERRS(kcnt,3,:) = NaN;
        SCORES_CERRS(kcnt,2,:) = NaN;
    end

    % compute information criteria, since this only does EM
    if (emconvg == 1)
        for cnt = 1:10
            if not(isempty(MMM.InfCrit{cnt}))
                switch MMM.PStype
                    case {'Gaussian','PEKern'}
                        SCORES_CERRS(kcnt,1,cnt) = feval(MMM.InfCrit{cnt},...
                            Data,best_clustassigns(:,kcnt,1)',MMM.regul_func,0);
                    case 'Kernel'
                        SCORES_CERRS(kcnt,1,cnt) = feval(MMM.InfCrit{cnt},...
                            Data,best_clustassigns(:,kcnt,1)',MMM.regul_func,MMM.htype,0);
                    case 'EC'
                        SCORES_CERRS(1,1,cnt) = feval(MMM.InfCrit{cnt},...
                            Data,best_clustassigns(:,kcnt,1)',MMM.regul_func,MMM.htype,0);                        
                end
            end
        end     % information criteria loop
    end

    % plot mixtures, if bivariate
    if (p == 2) && (emconvg == 1)
        fhclst = figure;
        % create the plot title in steps
        plottit = {['Structure Discovered By ',MMM.optim_type,'(',MMM.init_type,')'];...
            sprintf('k attempted/fit: %0.0f/%0.0f',MMM.Krange(kcnt),ksell)};
        if MMM.data_type == 1   % unknown
            PlotBivarClusts(Data,best_clustassigns(:,kcnt,1)),...
                title(plottit,'interpreter','none')
        else
            set(fhclst,'Position',[10,10,775,425]);
            subplot(1,2,1),PlotBivarClusts(Data,labels),title('Actual Structure')
            subplot(1,2,2),PlotBivarClusts(Data,best_clustassigns(:,kcnt,1)),...
                title({char(plottit);['Percent Correct: ',...
                num2str(SCORES_CERRS(kcnt,2,1),'%0.2f%%')]},'interpreter','none')
        end
        figure(fhclst), drawnow
        % save the figure as both a publication .eps and editable .fig/.ofig
        hgsave(fhclst,[savename,'_',num2str(k),'_',objec_func]);
        print(fhclst,[savename,'_',num2str(k),'_',objec_func,'.eps'],'-depsc');
    end

    MMM.lblTimIC = sprintf('IC: %0.2f',toc/60); 
    % End Expectation-Maximization EM algorithm
    
    MMM.lblTimK = sprintf('k: %0.2f',etime(clock,kstt)/60);
    MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
    Mixture_DispStatus
end                 % ks loop

% DISPLAY FINAL RESULTS BY INFORMATION CRITERIA
disp('FINAL RESULTS BY INFORMATION CRITERIA')
lin = '================================';
best_chroms = []; best_scores = [];
for infcnt = 1:10
    if isempty(MMM.InfCrit{infcnt}); continue; end;
    objec_func = MMM.InfCrit{infcnt};
    
    best_table = SCORES_CERRS(:,[4,1,2,3],infcnt);
    [best_score,best_k] = min(best_table(:,2));
    
    if (best_score == Inf)
        % Inf most likely means the optimization method was EM and it
        % didn't converge
        disp('Best Score is Infinity - It Seems That the Optimization Method Failed!!!')
        best_chroms = [best_chroms;[infcnt,best_clustassigns(:,best_k,infcnt)']];
        best_scores = [best_scores;[infcnt,best_table(best_k,:)]];
        continue
    end

    % since EM: for best k, copy first best_clustassigns into all others
    best_clustassigns(:,best_k,infcnt) = best_clustassigns(:,best_k,1);
    
    kselu = unique(best_clustassigns(:,best_k,infcnt)); ksell = length(kselu);    % account for the possibility of missing k
    % display time
    disp(' '), disp(lin)
    tab = table2str({'k';objec_func;'Percent Correct';'Number Clusters'},...
        best_table,{'%0.0f';'%0.2f';'%0.3f';'%0.0f'},1);
    tab(best_k + 3,1) = '*';
    disp(tab),disp(['* indicates best set of mixtures by minimizing ',objec_func])
    disp('Best Estimates From Best k:')
    % compute best estimates
    for mixcnt = 1:ksell
        ind = (best_clustassigns(:,best_k,infcnt) == kselu(mixcnt));
        if sum(ind) ~= 0
            mix_propors(mixcnt) = sum(ind)/n;
            if sum(ind) == 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0; ...
                    mix_ECparms(:,:,mixcnt) = 0;
            else
                if isequal(MMM.PStype,'EC')
                    [mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),...
                        mix_ECparms(:,:,mixcnt)]= MECEstimate(Data(ind,:),...
                        MMM.htype,MMM.regul_func,1);
                else
                    mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                    mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),...
                        MMM.regul_func,1,1,sum(ind));
                end
            end
        else
            mix_propors(mixcnt) = 0;
            mix_meanvecs(:,:,mixcnt) = nan(1,p);
            mix_covrmats(:,:,mixcnt) = nan(p,p);
        end
    end             % mixtures loop
    % display best estimates
    for mixcnt = 1:ksell
        disp(sprintf('Cluster %0.0f',kselu(mixcnt)))
        disp(sprintf('Mixing Proportion = %0.2f%%',100*mix_propors(mixcnt)))
        if (p <= maxP2disp) && not(isequal(MMM.PStype,'EC'))
            disp(DispMeanCovar(mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'%0.3f'))
        end
        switch MMM.PStype
            case 'Gaussian' % nothing special to do
            case 'Kernel'
                if p <= maxP2disp
                    [tmp,H] = MVKDE_Gauss(Data([best_clustassigns(:,best_k,infcnt) ...
                        == kselu(mixcnt)],:),zeros(1,p),MMM.htype,MMM.regul_func,0);
                    disp('Kernal Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'))
                end
            case 'EC'
                [mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),mix_ECparms(:,:,mixcnt)] = MECEstimate(Data([best_clustassigns(:,best_k,infcnt) ...
                    == kselu(mixcnt)],:),MMM.htype,MMM.regul_func,0);
                disp(DispMeanCovar(mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'%0.3f'))
                if isequal(MMM.htype,'KT')
                    disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',MMM.htype,mix_ECparms(:,:,mixcnt)))
                else
                    disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',MMM.htype,mix_ECparms(:,:,mixcnt)))
                end
            case 'PEKern'
                if sum(best_clustassigns(:,best_k,infcnt) == kselu(mixcnt)) > 1
                    [tmp,H,Bs] = MPEProdKernelPDF(Data([best_clustassigns(:,best_k,infcnt) ...
                        == kselu(mixcnt)],:),zeros(1,p),MMM.regul_func);
                else
                    H = zeros(1,p); Bs = H;
                end
                disp('PE Kernel Bandwidths'),disp(mat2str(H,2))
                disp('PE Kernel Betas'),disp(mat2str(Bs,2))
        end
    end             % mixtures loop
    if MMM.data_type ~= 1
        disp(lin), disp('Confusion Matrix from Best Estimates')
        ConfMatrix(labels,best_clustassigns(:,best_k,infcnt),0); % JAH 20120227
    end
    disp(lin)
    
    % now done with display, save stuff for the driver
    best_chroms = [best_chroms;[infcnt,best_clustassigns(:,best_k,infcnt)']];
    best_scores = [best_scores;[infcnt,best_table(best_k,:)]];
end                 % information criteria loop
    
if MMM.data_type == -1  % simulated
    disp(lin), disp('Actual Parameters')
    for mixcnt = 1:K
        disp(sprintf('Cluster %0.0f',mixcnt))
        disp(sprintf('Mixing Proportion = %0.2f%%',100*mixpropors(mixcnt)))
        % no reason to test for p <= maxP2disp, since I'll never do a sim that big
        disp(DispMeanCovar(meanvecs(:,:,mixcnt),covrmats(:,:,mixcnt),'%0.3f'))
        disp(sprintf('MPE Kurtosis Parameter = %0.2f',betas(mixcnt)))
    end
    disp(lin)
end

tottim = etime(clock,drvstt);
disp(sprintf('Mixture_KM Completed in \n\t%1.0f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',tottim./[1,60,3600]));
MMM.lblTimRep = sprintf('Replic: %0.1f',tottim/60);
MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
Mixture_DispStatus
diary off

clear fh* H* X bins col dat ind tab posteriors tmp rng kselu kcnt best_chrom
clear termcount mix_* minval mixcnt mutcnt dist_from_clusts Knum maxs cm
clear mutprobs new_pop num_mut pop_fitness population minind best_k kstt
clear std_fitness val sumks xoverpoint stdindex lin k tab drvstt val Krange
clear ksell popcnt datcnt gencnt cnt GKM_toc best_table emconvg objec_func
clear genscores confmat popul_size_save infcnt best_score;
clear global ys
save([savename,'.mat'])


%{
JAH 20070309 (copied from Mixture_GKM), adapted for octave 3.4.3

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