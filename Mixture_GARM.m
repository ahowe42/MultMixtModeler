%{
  Usage: THIS MUST BE CALLED FROM drv_Mixture
  This script will perform one run of mixture modeling using the genetic
  algorithm with regularized mahalanobis distance to initialize either the 
  EM or genetic EM algorithm.

  Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe
%}
  
drvstt = clock; rnd_stat = sum(drvstt*1000000); rand('state',rnd_stat);
stt = sprintf('%4.0f%02.0f%02.0f_%02.0f%02.0f%02.0f',drvstt);
savename = [mydir,dir,MMM.init_type,'+',MMM.optim_type,'_',stt]; diary([savename,'.out']);
disp(repmat('#',1,50)), Mixture_DispParms
disp(sprintf('Random State: %0.0f',rnd_stat))
disp(['Save files: ',savename([(length(mydir)+1):length(savename)])]), disp(repmat('#',1,50))

% initialize overall "things"
maxP2disp = 8;      % if p > this, don't display estimated mu, sigma, h
popul_size_save = MMM.popul_size; I = zeros(p);
% SCORES_CERRS = [score,%_correct,actual_k,attempt_k]
SCORES_CERRS = ones(MMM.Knum,4,10)*Inf;
SCORES_CERRS(:,2,:) = NaN;                  % initialize %_correct for data_type = 1
SCORES_CERRS(:,4,:) = repmat(MMM.Krange',1,10); % and attempt_k
best_clustassigns = zeros(n,MMM.Knum,10); global ys; ys = NaN;  % for ICOMP_PERF with k = 1;

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
            [mix_meanvecs(:,:,1),mix_covrmats(:,:,1),mix_ECparms(:,:,1)]= ...
                MECEstimate(Data,MMM.htype,MMM.regul_func,0);
        else
            mix_meanvecs(:,:,1) = mean(Data,1);
            mix_covrmats(:,:,1) = CovSmooth(Data,MMM.regul_func,1,0,n);
        end
        % display final parameters
        disp(sprintf('Cluster %0.0f',1))
        disp(['pi = ',num2str(100*mix_propors(1),'%0.2f%%')])
        if p <= maxP2disp
            disp(DispMeanCovar(mix_meanvecs(:,:,1),mix_covrmats(:,:,1),'%0.3f'))
        end
        switch MMM.PStype
            case 'Gaussian' % nothing special
            case 'Kernel'
                if p <= maxP2disp
                    [tmp,H] = MVKDE_Gauss(Data,zeros(1,p),MMM.htype,MMM.regul_func,0);                
                    disp('Kernel Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'))
                end
            case 'EC'
                if isequal(MMM.htype,'KT')
                    disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',...
                        MMM.htype,mix_ECparms(:,:,1)))
                else
                    disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',MMM.htype,...
                        mix_ECparms(:,:,1)))
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
                    case {'Kernel','EC'}
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
    % Begin Genetic Algorithm on Regularized Mahalanobis Distance (GARM)
    disp(sprintf('GARM Algorithm Entered for k = %0.0f',k))
    population = unidrnd(k,MMM.popul_size,n);   % initialize chromosomes
    tos = fix(n*[1:k]/k); frs = [1,1+tos([1:(k-1)])];
    for cnt = 1:k; population(1,[frs(cnt):tos(cnt)]) = cnt; end    

    % initialize more "things"
    mix_meanvecs = zeros(1,p,k); mix_covrmats = zeros(p,p,k); sumks = sum([1:k]);
    mix_propors = zeros(1,k); regmahal_from_clusts = mix_propors;
    genscores = zeros(MMM.num_generns,2); termcount = 0;
    mix_ECparms = zeros(1,2+isequal(MMM.htype,'KT'),k);

    for gencnt = 1:MMM.num_generns
        num_mut = ceil(MMM.prob_mutate*MMM.popul_size); % number of chromosomes in population to mutate    

        % CHECK FOR ILLEGAL CHROMOSOMES & COMPUTE OBJECTIVE FUNCTION VALUES
        pop_fitness = ones(MMM.popul_size,1)*Inf;
        for popcnt = 1:MMM.popul_size
            % ensure each chromosome has at least p members of each cluster
            popclusts = unique(population(popcnt,:));
            if (sum(popclusts) ~= sumks)
                for mixcnt = 1:k
                    % this mixture missing, so psuedorandomly assign p datapoints to it
                    % to ensure a random assignment doesn't get overwritten, do so in 1/ks
                    if (sum(popclusts == mixcnt) == 0)
                        rng = min(n,ceil(n*[(mixcnt - 1):mixcnt]/k)) + [1,0];
                        population(popcnt,unidrnd(rng(2) - rng(1),p,1) + rng(1)) = mixcnt;
                    end
                end % mixtures loop
            end        
            % evaluate the total within-cluster regularized mahalanobis distance
            pop_fitness(popcnt) = KMix_TotInClustRMDist(Data,population(popcnt,:),MMM.regul_func,MMM.regul_scale);
        end         % chromosomes loop
        [minval,minind]= min(pop_fitness);
        genscores(gencnt,:) = [min(pop_fitness),mean(pop_fitness(not(isnan(pop_fitness))))];
        if (gencnt == MMM.num_generns); break; end;       % don't bother with the next generation

        % EARLY TERMINATION ALLOWED?
        if genscores(gencnt,1) < best_score
            best_score = genscores(gencnt,1);
            best_chrom = population(minind,:);
            termcount = 1;
        elseif (genscores(gencnt,1) > best_score) && not(MMM.elitism)
            % if elitism is on, we can still do early termination with this
            termcount = termcount + 1;
        elseif genscores(gencnt,1) == best_score
            termcount = termcount + 1;
        end
        if termcount >= MMM.nochange_terminate
            disp(['Early Termination On Generation ',num2str(gencnt),' of ',num2str(MMM.num_generns)]);
            genscores = genscores([1:gencnt],:);
            break;
        end

        % SELECTION OF NEXT GENERATION
        % sort scores appropriately
        [val, stdindex] = sort(pop_fitness);
        % prepare bins for roulette - bigger bins at the beginning with lower scores
        bins =  cumsum([MMM.popul_size:-1:1]/(MMM.popul_size*(MMM.popul_size + 1)/2))';
        % roulette selection - find into which bin the random falls
        new_pop = sum(repmat(rand(1,MMM.popul_size),MMM.popul_size,1) >= repmat(bins,1,MMM.popul_size))+1;
        new_pop = population(stdindex(new_pop),:);

        % CHROMOSOME MUTATION
        popmut = unidrnd(MMM.popul_size,num_mut,1);
        for mutcnt = 1:num_mut
            % first we need to compute the mean vectors and covariance matrices
            for mixcnt = 1:k
                ind = new_pop(popmut(mutcnt),:) == mixcnt;
                if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                    mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0;
                else
                    mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);                
                    mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),MMM.regul_func,1,0,sum(ind));
                end                
            end     % mixtures loop

            ind = find(rand(n,1) <= MMM.prob_mutate);  % randomly choose datapoints to mutate
            for datcnt = 1:length(ind)
                dat = Data(ind(datcnt),:);
                % compute regularized mahalanobis distances from each centroid
                for mixcnt = 1:k                
                    if isequal(mix_covrmats(:,:,mixcnt),I)
                        % singleton or no variation cluster, distance is 0
                        regmahal_from_clusts(mixcnt) = 0;
                    else
                        regmahal_from_clusts(mixcnt) = RegMahalanobis(dat,...
                            mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'MLE',MMM.regul_scale,0);
                    end
                end % mixtures loop
                % compute mutation probabilities
                regmahal_from_clusts = max(regmahal_from_clusts) - regmahal_from_clusts;
                mutprobs = regmahal_from_clusts / sum(regmahal_from_clusts);

                % choose into what mixture (if any) to mutate
                mutprobs = max(find((mutprobs - rand(1,k)) > 0));
                if (~isempty(mutprobs)) && (new_pop(popmut(mutcnt),ind(datcnt)) ~= mutprobs)
                    val = new_pop(popmut(mutcnt),ind(datcnt));      % current label
                    new_pop(popmut(mutcnt),ind(datcnt)) = mutprobs; % reassign label
                    % recompute parameters for both affected clusters
                    tmp = (new_pop(popmut(mutcnt),:) == val);
                    mix_propors(val) = sum(tmp)/n;
                    if mix_propors(val)*n <= 1
                        mix_covrmats(:,:,val) = 0; mix_meanvecs(:,:,val) = 0;
                    else
                        mix_meanvecs(:,:,val) = mean(Data(tmp,:),1);
                        mix_covrmats(:,:,val) = CovSmooth(Data(tmp,:),MMM.regul_func,1,0,sum(tmp));
                    end
                    tmp = (new_pop(popmut(mutcnt),:) == mutprobs);
                    mix_propors(mutprobs) = sum(tmp)/n;
                    if mix_propors(mutprobs)*n <= 1
                        mix_covrmats(:,:,mutprobs) = 0; mix_meanvecs(:,:,mutprobs) = 0;
                    else
                        mix_meanvecs(:,:,mutprobs) = mean(Data(tmp,:),1);
                        mix_covrmats(:,:,mutprobs) = CovSmooth(Data(tmp,:),MMM.regul_func,1,0,sum(tmp));
                    end
                end
            end     % datapoints loop
        end         % chromosomes to mutate loop

        % CROSSOVER OPERATION ON NEW POPULATION
        new_pop = new_pop(randperm(MMM.popul_size),:);  % randomly permute rows
        for popcnt = 2:2:MMM.popul_size
            if rand <= MMM.prob_xover
                xoverpoint = unidrnd(n - 2) + 1;  % ensure xover point not on ends
                tmp = new_pop(popcnt - 1,:);
                % trade right-most portions
                new_pop(popcnt-1,[(xoverpoint + 1):n]) = new_pop(popcnt,[(xoverpoint + 1):n]);
                new_pop(popcnt,[(xoverpoint + 1):n]) = tmp([(xoverpoint + 1):n]);
            end
        end         % chromosomes loop

        % CHECK FOR ILLEGAL CHROMOSOMES
        for popcnt = 1:MMM.popul_size
            % ensure each chromosome has at least 1 member of each cluster
            popclusts = unique(new_pop(popcnt,:));
            if (sum(popclusts) ~= sumks)
                for mixcnt = 1:k
                    % this mixture missing, so psuedorandomly assign p datapoints to it
                    % to ensure a random assignment doesn't get overwritten, do so in 1/ks
                    if (sum(popclusts == mixcnt) == 0)                    
                        rng = min(n,ceil(n*[(mixcnt - 1):mixcnt]/k)) + [1,0];
                        new_pop(popcnt,unidrnd(rng(2) - rng(1),p,1) + rng(1)) = mixcnt;
                    end
                end % mixtures loop
            end
        end         % chromosomes loop

        % MAHALANOBIS OPERATION
        % choose chromosomes without replacement - must do this to avoid reevaluation of empty clusters
        popmut = randperm(MMM.popul_size); popmut = popmut([1:num_mut]);
        for mutcnt = 1:num_mut
            % first we need to compute the mean vectors and covariance matrices (again . . .)
            for mixcnt = 1:k
                ind = new_pop(popmut(mutcnt),:) == mixcnt;                
                if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                    mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0;
                else
                    mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                    mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),MMM.regul_func,1,0,sum(ind));
                 end
            end     % mixtures loop

            for datcnt = 1:n
                % compute regularized mahalanobis distances from each centroid
                for mixcnt = 1:k
                    if isequal(mix_covrmats(:,:,mixcnt),I)
                        % singleton or no variance cluster, distance is 0
                        regmahal_from_clusts(mixcnt) = 0;
                    else
                        regmahal_from_clusts(mixcnt) = RegMahalanobis(Data(datcnt,:),...
                            mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'MLE',MMM.regul_scale,0);
                    end
                end % mixtures loop            
                % assign datapoint to nearest cluster
                [val,ind] = min(regmahal_from_clusts);
                if (new_pop(popmut(mutcnt),datcnt) ~= ind)
                    val = new_pop(popmut(mutcnt),datcnt);   % current label
                    new_pop(popmut(mutcnt),datcnt) = ind;   % reassign
                    % recompute parameters for both affected clusters
                    tmp = (new_pop(popmut(mutcnt),:) == val);
                    mix_propors(val) = sum(tmp)/n;
                    if mix_propors(val)*n <= 1
                        mix_covrmats(:,:,val) = 0; mix_meanvecs(:,:,val) = 0;
                    else
                        mix_meanvecs(:,:,val) = mean(Data(tmp,:),1);
                        mix_covrmats(:,:,val) = CovSmooth(Data(tmp,:),MMM.regul_func,1,0,sum(tmp));
                     end
                    tmp = (new_pop(popmut(mutcnt),:) == ind);
                    mix_propors(ind) = sum(tmp)/n;
                    if mix_propors(ind)*n <= 1
                        mix_covrmats(:,:,ind) = 0; mix_meanvecs(:,:,ind) = 0;
                    else
                        mix_meanvecs(:,:,ind) = mean(Data(tmp,:),1);
                        mix_covrmats(:,:,ind) = CovSmooth(Data(tmp,:),MMM.regul_func,1,0,sum(tmp));
                    end
                end
            end     % datapoints loop
        end         % chromosomes to mahalanobis operate

        % CONVEY BEST INDIVIDUAL INTO NEW POPULATION
        if MMM.elitism
            population = [new_pop;population(stdindex(1),:)];
        else
            population = new_pop;
        end
        MMM.popul_size = size(population,1);
        if rem(gencnt,2) == 1
            disp(sprintf('GARM Generation %0.0f of %0.0f: Best Score = %0.4f, %0.0f', gencnt,MMM.num_generns,best_score,termcount));
        end
    end             % generations loop
    MMM.popul_size = popul_size_save;
    ys = best_chrom; % here for ICOMP_PERF functions

    % Begin Expectation-Maximization (EM / GEM)
    % prepare for the specific optimizer
    if isequal(MMM.optim_type,'EM')    % EM
        % nothing to do ...
    else                        % GEM
        if not(MMM.elitism) % if elitism off, no guarantee best solution in final pop, so add it
            pop_fitness = [pop_fitness;best_score]; population = [population;best_chrom];
        end
        % start GEM with roulette selection of final population chroms
        [val, stdindex] = sort(pop_fitness);           % sort scores appropriately
        % prepare bins for roulette - bigger bins at the beginning with lower scores
        bins =  cumsum(2*[popul_size_save:-1:1]/(popul_size_save*(popul_size_save + 1)))';
        % roulette selection - find into which bin the random falls
        new_pop = sum(repmat(rand(1,popul_size_save),popul_size_save,1) >= repmat(bins,1,popul_size_save)) + 1;
        new_pop = population(stdindex(new_pop),:);
    end
    GARM_toc = toc;
    disp(sprintf('GARM Completed in \n\t%1.0f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',GARM_toc./[1,60,3600]));
    MMM.lblTimInit = [MMM.init_type,sprintf(': %0.2f',GARM_toc/60)];
    MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
    Mixture_DispStatus
    % End Genetic Algorithm on Regularized Mahalanobis Distance (GARM)    

    for infcnt = 1:10
        if isempty(MMM.InfCrit{infcnt}); continue; end;
        tic;
        
        % disply progress info
        if isequal(MMM.optim_type,'EM'); objec_func = 'Log-Likelihood';
        else; objec_func = MMM.InfCrit{infcnt}; end;
        MMM.lblCurrIC = objec_func;
        Mixture_DispStatus
        disp(sprintf([MMM.optim_type,' Algorithm Entered for k = %0.0f: ',objec_func],k))
        
        % actually do the optimization
        if isequal(MMM.optim_type,'EM')
            if isequal(MMM.PStype, 'EC') || isequal(MMM.PStype, 'Kernel')
                [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
                    Data,best_chrom,[MMM.convgcrit,MMM.maxiter],MMM.htype,MMM.regul_func,...
                    [savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
            else
                [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
                    Data,best_chrom,[MMM.convgcrit,MMM.maxiter],MMM.regul_func,...
                    [savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
            end
            emconvg = not(isnan(mix_propors(1)));
        else
            if isequal(MMM.PStype, 'EC') || isequal(MMM.PStype, 'Kernel')
                [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
                    Data,[popul_size_save,MMM.num_generns,MMM.elitism,MMM.prob_xover,...
                    MMM.prob_mutate,MMM.nochange_terminate],MMM.htype,objec_func,1,new_pop,...
                    MMM.regul_func,[savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
            else
                [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = feval(MMM.optim_func,...
                    Data,[popul_size_save,MMM.num_generns,MMM.elitism,MMM.prob_xover,...
                    MMM.prob_mutate,MMM.nochange_terminate],objec_func,1,new_pop,...
                    MMM.regul_func,[savename,'_',num2str(k)],MMM.init_type,MMM.showplot);
            end
            emconvg = 1;    % just here to work with EM
        end
        
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
            end     % mixtures loop
            % assign based on max matches
            for mixcnt = 1:K
                val = max(maxs(2,maxs(1,:) == mixcnt));
                if not(isempty(val))
                    ind = maxs(3,and(maxs(1,:) == mixcnt, maxs(2,:) == val));
                    bins(bins == ind(1)) = 100 + mixcnt;
                end
            end     % mixtures loop
            % see if we need to do more clean up
            for mixcnt = 1:K
                cm = (sum(bins == mixcnt) > 0);
                if cm; break; end;
            end     % mixtures loop
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
            best_clustassigns(:,kcnt,infcnt) = bins;
        elseif  emconvg == 1
            [val,best_clustassigns(:,kcnt,infcnt)] = max(posteriors,[],2);
        else
            best_clustassigns(:,kcnt,infcnt) = NaN; % em did not converge
        end
        kselu = unique(best_clustassigns(:,kcnt,infcnt)); ksell = length(kselu);    % account for the possibility of missing k
        
        % recompute (in case assignments changed) and display estimates
        disp(sprintf(['\nResults for k = %0.0f: ',objec_func],k))
        if (emconvg == 1) && (MMM.data_type ~= 1)    
            for mixcnt = 1:ksell
                ind = (best_clustassigns(:,kcnt,infcnt) == kselu(mixcnt));
                if sum(ind) ~= 0
                    mix_propors(mixcnt) = sum(ind)/n;
                    if sum(ind) == 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                        mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0; mix_ECparms(:,:,mixcnt) = 0;
                    else
                        if isequal(MMM.PStype,'EC')
                            [mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),...
                                mix_ECparms(:,:,mixcnt)]= MECEstimate(Data(ind,:),...
                                MMM.htype,MMM.regul_func,0);
                        else
                            mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                            mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),...
                                MMM.regul_func,1,0,sum(ind));
                        end
                    end
                else
                    mix_propors(mixcnt) = 0;
                    mix_meanvecs(:,:,mixcnt) = NaN;
                    mix_covrmats(:,:,mixcnt) = NaN;
                    mix_ECparms(:,:,mixcnt) = NaN;
                end
            end     % mixtures loop
        end
        
        if (emconvg == 1)
            for mixcnt = 1:ksell
                disp(sprintf('Cluster %0.0f',kselu(mixcnt)))
                disp(['pi = ',num2str(100*mix_propors(mixcnt),'%0.2f%%')])
                if p <= maxP2disp
                    disp(DispMeanCovar(mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'%0.3f'))
                end
                switch MMM.PStype
                    case 'Gaussian' % nothing special
                    case 'Kernel'
                        if p <= maxP2disp
                            [tmp,H] = MVKDE_Gauss(Data([best_clustassigns(:,kcnt,infcnt) ...
                                == kselu(mixcnt)],:),zeros(1,p),MMM.htype,MMM.regul_func,0);
                            disp('Kernel Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'));
                        end
                    case 'EC'
                        if isequal(MMM.htype,'KT')
                            disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',...
                                MMM.htype,mix_ECparms(:,:,mixcnt)))
                        else
                            disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',...
                                MMM.htype,mix_ECparms(:,:,mixcnt)))
                        end
                    case 'PEPK'
                        [tmp,H,Bs] = MPEProdKernelPDF(Data([best_clustassigns(:,kcnt,infcnt) ...
                            == kselu(mixcnt)],:),zeros(1,p),MMM.regul_func);
                        disp('PE Kernel Bandwidths'),disp(mat2str(H,2))
                        disp('PE Kernel Betas'),disp(mat2str(Bs,2))                        
                end
            end         % mixtures loop
            % compute Percent Correct & objective value, do confusion matrix
            switch MMM.PStype
                case{'Gaussian','PEKern'}
                    SCORES_CERRS(kcnt,1,infcnt) = feval(MMM.InfCrit{infcnt},Data,...
                        best_clustassigns(:,kcnt,infcnt)',MMM.regul_func,0);
                case {'Kernel','EC'}
                    SCORES_CERRS(kcnt,1,infcnt) = feval(MMM.InfCrit{infcnt},Data,...
                        best_clustassigns(:,kcnt,infcnt)',MMM.regul_func,MMM.htype,0);
            end
            SCORES_CERRS(kcnt,3,infcnt) = ksell;
            if MMM.data_type ~= 1   % not unknown
                SCORES_CERRS(kcnt,2,infcnt) = 100*sum(best_clustassigns(:,kcnt,infcnt) == labels)/n;
                if emconvg == 1
                    confmat = ConfMatrix(labels,best_clustassigns(:,kcnt,infcnt),0); % JAH 20120227
                    disp(sprintf('Percent Correct = %0.2f%%\n',SCORES_CERRS(kcnt,2,infcnt)))
                end
            end
        else
            disp([MMM.optim_func,' FAILED TO CONVERGE - NO ESTIMATES!'])
            SCORES_CERRS(kcnt,1,:) = Inf;
            SCORES_CERRS(kcnt,3,:) = NaN;
            SCORES_CERRS(kcnt,2,:) = NaN;
        end
        
        % if EM, only need to do once, just copy results into all IC
        if (infcnt == 1) && isequal(MMM.optim_type,'EM') && (emconvg == 1)
            for cnt = 1:10
                if not(isempty(MMM.InfCrit{cnt}))
                    switch MMM.PStype
                        case{'Gaussian','PEKern'}
                            SCORES_CERRS(kcnt,1,cnt) = feval(MMM.InfCrit{cnt},...
                                Data,best_clustassigns(:,kcnt,1)',MMM.regul_func,0);
                        case {'Kernel','EC'}
                            SCORES_CERRS(kcnt,1,cnt) = feval(MMM.InfCrit{cnt},...
                                Data,best_clustassigns(:,kcnt,1)',MMM.regul_func,MMM.htype,0);
                    end
                end
            end     % information criteria loop
            SCORES_CERRS(kcnt,3,:) = SCORES_CERRS(kcnt,3,1);            
            if MMM.data_type ~= 1   % not unknown
                SCORES_CERRS(kcnt,2,:) = SCORES_CERRS(kcnt,2,1);
            end            
        end        
        
        % plot mixtures, if bivariate
        if (p == 2) && (emconvg == 1)
            fhclst = figure;
            % create the plot title in steps
            plottit = sprintf('Structure Discovered By %s(%s)\nk attempted/fit: %0.0f/%0.0f', MMM.optim_type,MMM.init_type,MMM.Krange(kcnt),ksell);
            if not(isequal(MMM.optim_type,'EM'))
                plottit = sprintf('%s\n%s Value: %0.4',plottit,objec_func,SCORES_CERRS(kcnt,1,infcnt)); % JAH 20120324
            end
            if MMM.data_type == 1   % unknown
                PlotBivarClusts(Data,best_clustassigns(:,kcnt,infcnt))
                title(plottit,'interpreter','none')
            else
                set(fhclst,'Position',[10,10,775,425]);
                subplot(1,2,1),PlotBivarClusts(Data,labels),title('Actual Structure')
                subplot(1,2,2),PlotBivarClusts(Data,best_clustassigns(:,kcnt,infcnt))
				title(sprintf('%s\nPercent Correct: %0.2f%%',plottit,SCORES_CERRS(kcnt,2,infcnt)),'interpreter','none'); % JAH 20120324
            end
            figure(fhclst), drawnow
            % save the figure as both a publication .eps and editable .fig/.ofig
            hgsave(fhclst,[savename,'_',num2str(k),'_',objec_func]);
			print(fhclst,[savename,'_',num2str(k),'_',objec_func,'.eps'],'-depsc');
        end
        MMM.lblTimIC = sprintf('IC: %0.2f',toc/60);
        MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
        Mixture_DispStatus
        if isequal(MMM.optim_type,'EM'); break; end;
    end             % information criteria loop
    % End Expectation-Maximization (EM / GEM) algorithm
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

    % if EM: for best k, copy first best_clustassigns into all others
    if isequal(MMM.optim_type,'EM')
        best_clustassigns(:,best_k,infcnt) = best_clustassigns(:,best_k,1);
    end
    
    kselu = unique(best_clustassigns(:,best_k,infcnt)); ksell = length(kselu);    % account for the possibility of missing k
    % display time
    disp(' '), disp(lin)
    tab = table2str({'k';objec_func;'Percent Correct';'Number Clusters'},...
        best_table,{'%0.0f';'%0.2f';'%0.3f';'%0.0f'},1);
    tab(best_k + 3,1) = '*';
    disp(tab),disp(['* indicates best set of mixtures by minimizing ',objec_func])
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
            mix_meanvecs(:,:,mixcnt) = NaN;
            mix_covrmats(:,:,mixcnt) = NaN;
            mix_ECparms(:,:,mixcnt) = NaN;
        end
    end             % mixtures loop
    % display best estimates
    disp('Best Estimates From Best k:')
    for mixcnt = 1:ksell
        disp(sprintf('Cluster %0.0f',kselu(mixcnt)))
        disp(sprintf('Mixing Proportion = %0.2f%%',100*mix_propors(mixcnt)))
        if p <= maxP2disp
            disp(DispMeanCovar(mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt),'%0.3f'))
        end
        switch MMM.PStype
            case 'Gaussian' % nothing special to do
            case 'Kernel'
                if p <= maxP2disp
                    [tmp,H] = MVKDE_Gauss(Data([best_clustassigns(:,best_k,infcnt) ...
                        == kselu(mixcnt)],:),zeros(1,p),MMM.htype,MMM.regul_func,1);
                    disp('Kernal Bandwidth Matrix'),disp(MatrixtoStr(H,'%0.3f'))
                end
            case 'EC'
                if isequal(MMM.htype,'KT')
                    disp(sprintf('%s Parameters\nN=%0.2f,r=%0.2f,beta=%0.2f',...
                        MMM.htype,mix_ECparms(:,:,mixcnt)))
                else
                    disp(sprintf('%s Parameters\nN=%0.2f,nu=%0.2f',MMM.htype,...
                        mix_ECparms(:,:,mixcnt)))
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
disp(sprintf('Mixture_GARM Completed in \n\t%1.0f Seconds\n\t%1.4f Minutes\n\t%1.4f Hours',tottim./[1,60,3600]));
MMM.lblTimRep = sprintf('Replic: %0.1f',tottim/60);
MMM.lblTimTot = sprintf('Total (m): %1.4f',etime(clock,MMM.totstt)/60);
Mixture_DispStatus
diary off

clear fh* H* X GARM_toc bins col dat ind tab posteriors tmp rng I cm kcnt
clear termcount cnt minval mixcnt mutcnt regmahal_from_clusts Knum kselu
clear mutprobs new_pop num_mut pop_fitness popclusts population minind popmut
clear std_fitness val sumks xoverpoint stdindex lin k tab drvstt val Krange
clear ksell popcnt datcnt gencnt infcnt best_score maxs objec_func best_k kstt
clear best_chrom best_table mix_* emconvg genscores confmat popul_size_save
clear global ys bs tos frs
save([savename,'.mat'])

%{
JAH GARM began on 20061007; EM began on 20061025; GEM began on 20061203, adapted for octave 3.4.3 20120324

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