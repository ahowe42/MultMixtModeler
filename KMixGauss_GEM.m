function [posteriors,mix_propors,mix_meanvecs,mix_covrmats] = KMixGauss_GEM(Data,GEMparams,objec_func,ic,population,smth,savename,init,pltflg)
%{ 
  [posteriors, mixing proportions, mean vectors, covar matrices] = KMixGauss_GEM(
  data, GA parameters, objec funct, inf_crit_flag, init popul, ...
  covar smoother, figure name, init. method, plot flag)
   This implements the genetic EM algorithm for Gaussian mixtures.  It
   returns the posterior probabilities of group belonging and estimates
   of the mean vector and covariance matrix per cluster.  Pass in a
   smoothing code for CovSmooth for cases in which a cluster's covariance
   matrix is ill-conditioned.

   Where
   data --- (nxp) matrix of data
   GA Parameters --- [population size, number generations, elitism, ...
        xover rate, mutate rate, early termination threshold]
   Objective function --- fitness function for GA
   inf_crit_flag --- 1 = information criteria so maximize neg obj func value,
        otherwise, just maximize obj func value.
   initial population --- (population size x length(data)) of 1:k
   covariance smoother --- alpha code to pass covsmooth
   figure name --- name to save figure, include full path; if you don't
      want the plot saved as a file, pass in 0.
   init. method --- string of initialization method (GARM,GKM,...)
   plot flag --- 0 = show only at end; 1 = update on-the-fly
   posteriors --- (nxk) matrix of belonging probabilities
   mixing proportions --- (1xk) vector of mixing proportions
   mean vectors --- (1,p,k) matrix of mean vector per cluster
   covar matrices --- (p,p,k) matrix of covariance matrix per cluster
   
  Copyright (C) 2006 Prof. Hamparsum Bozdogan & J. Andrew Howe; see below
%}

[n,p] = size(Data); k = max(population(:));

if (nargin ~= 9) || (length(GEMparams) < 6) || (sum(size(population) == [GEMparams(1),n]) ~= 2) || (not(isequal([1:k],reshape(unique(population(:)),1,k))))
    % wrong # of args, wrong # of GA params, initial population wrong size, not all k represented
    fprintf('KMixGauss_GEM: INVALID USAGE-Please read the following instructions!\n'), help KMixGauss_GEM, return
end

if ic == 1; ic = -1; else ic = 1; end;
% extract the GA parameters
popul_size = GEMparams(1); num_generns = GEMparams(2); elitism = GEMparams(3);
prob_xover = GEMparams(4); prob_mutate = GEMparams(5); nochange_terminate = GEMparams(6);
% initialize things
mix_meanvecs = zeros(1,p,k); mix_covrmats = zeros(p,p,k); sumks = sum([1:k]);
mix_propors = zeros(1,k); posteriors = mix_propors; fhga = figure; I = zeros(p);
genscores = zeros(num_generns,2); best_score = -Inf; termcount = 0; pltxyy = [];

% run the GA
for gencnt = 1:num_generns
    num_mut = ceil(prob_mutate*popul_size); % number of chromosomes in population to mutate    

    % COMPUTE OBJECTIVE FUNCTION VALUES
    pop_fitness = ones(popul_size,1)*Inf;
    for popcnt = 1:popul_size
        % evaluate the objective function for each chromosome
        % note we are maximizing the negative of the information criteria
        pop_fitness(popcnt) = ic*feval(objec_func,Data,population(popcnt,:),smth,1);
    end         % chromosomes loop
    [maxval,maxind]= max(pop_fitness);
    genscores(gencnt,:) = [max(pop_fitness),mean(pop_fitness(not(isnan(pop_fitness))))];

    % UPDATE DATA FOR PLOT
    pltxyy = [pltxyy;[gencnt,genscores(gencnt,1),genscores(gencnt,2)]];
    if (pltflg == 1)
        figure(fhga)
        [AX,H1,H2] = plotyy(pltxyy(:,1),pltxyy(:,2),pltxyy(:,1),pltxyy(:,3)); xlabel('Generation');
        title(['GEM(',init,') Progress: Objective function ',objec_func],'interpreter','none');
        set(get(AX(1),'Ylabel'),'String','Maximum Value (o)','color','b');
        set(H1,'color','b','marker','o'); set(AX(2),'ycolor','b');
        set(get(AX(2),'Ylabel'),'String','Average Value (*)','color','r');
        set(AX(2),'ycolor','r'); set(H2,'color','r','marker','*');
        drawnow
    end

    if (gencnt == num_generns); break; end;       % don't bother with the next generation

    % EARLY TERMINATION ALLOWED?
    if genscores(gencnt,1) > best_score
        best_score = genscores(gencnt,1);
        best_chrom = population(maxind,:);
        termcount = 1;
    elseif (genscores(gencnt,1) < best_score) && not(elitism)
        % if elitism is on, we can still do early termination with this
        termcount = termcount + 1;
    elseif genscores(gencnt,1) == best_score
        termcount = termcount + 1;
    end
    if termcount >= nochange_terminate
        disp(['Early Termination On Generation ',num2str(gencnt),' of ',num2str(num_generns)]);
        genscores = genscores([1:gencnt],:);
        break;
    end
    

    % SELECTION OF NEXT GENERATION
    % sort scores appropriately
    [std_fitness, stdindex] = sort(pop_fitness,'descend');
    % prepare bins for roulette
    bins =  cumsum([popul_size:-1:1]/(popul_size*(popul_size + 1)/2))';
    % roulette selection - find into which bin the random falls
    new_pop = sum(repmat(rand(1,popul_size),popul_size,1) >= repmat(bins,1,popul_size))+1;
    new_pop = population(stdindex(new_pop),:);

    % CHROMOSOME MUTATION
    popmut = unidrnd(popul_size,num_mut,1);
    for mutcnt = 1:num_mut
        % first we need to compute the mean vectors and covariance matrices
        for mixcnt = 1:k
            ind = new_pop(popmut(mutcnt),:) == mixcnt;
            mix_propors(mixcnt) = sum(ind)/n;            
            if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0;
            else
                mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                % estimate the covariance matrix: if ill-conditioned, use the smoother
                mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),smth,1,1,sum(ind));
                % if still ill-conditioned, make all zeros
                if MatrixProblem(mix_covrmats(:,:,mixcnt)) ~= 0; mix_covrmats(:,:,mixcnt) = 0; end;
            end
        end     % mixtures loop

        ind = find(rand(n,1) <= prob_mutate);  % randomly choose datapoints to mutate
        for datcnt = 1:length(ind)
            dat = Data(ind(datcnt),:);
            % compute posterior probability of cluster memberships
            for mixcnt = 1:k
                if isequal(mix_covrmats(:,:,mixcnt),I)
                    % singleton, let pi be 0
                    posteriors(mixcnt) = 0;
                else
                    posteriors(mixcnt) = mix_propors(mixcnt)*mvnpdf(dat,...
                        mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt));
                end
            end % mixtures loop
            posteriors = posteriors/sum(posteriors);
            % compute mutation probabilities
            posteriors = max(posteriors) - posteriors;
            
            % only change if not if all clusters equally likely
            if sum(posteriors) ~= 0
                mutprobs = posteriors / sum(posteriors);
                % choose into what mixture (if any) to mutate
                mutprobs = max(find((mutprobs - rand(1,k)) < 0));
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
                        % estimate the covariance matrix: if ill-conditioned, use the smoother
                        % JAH corrected from mixcnt to val 20080831
                        mix_covrmats(:,:,val) = CovSmooth(Data(tmp,:),smth,1,1,sum(tmp));
                        % if still ill-conditioned, make all zeros
                        if MatrixProblem(mix_covrmats(:,:,val)) ~= 0; mix_covrmats(:,:,val) = 0; end;
                    end
                    tmp = (new_pop(popmut(mutcnt),:) == mutprobs);
                    mix_propors(mutprobs) = sum(tmp)/n;
                    if mix_propors(mutprobs)*n <= 1
                        mix_covrmats(:,:,mutprobs) = 0; mix_meanvecs(:,:,mutprobs) = 0;
                    else
                        mix_meanvecs(:,:,mutprobs) = mean(Data(tmp,:),1);
                        % estimate the covariance matrix: if ill-conditioned, use the smoother
                        % JAH corrected from mixcnt to mutprobs 20080831
                        mix_covrmats(:,:,mutprobs) = CovSmooth(Data(tmp,:),smth,1,1,sum(tmp));
                        % if still ill-conditioned, make all zeros
                        if MatrixProblem(mix_covrmats(:,:,mutprobs)) ~= 0; mix_covrmats(:,:,mutprobs) = 0; end;
                    end
                end
            end
        end     % datapoints loop
    end         % chromosomes to mutate loop

    % CROSSOVER OPERATION ON NEW POPULATION
    new_pop = new_pop(randperm(popul_size),:);  % randomly permute rows
    for popcnt = 2:2:popul_size
        if rand <= prob_xover
            xoverpoint = unidrnd(n - 2) + 1;  % ensure xover point not on ends
            tmp = new_pop(popcnt - 1,:);
            % trade right-most portions
            new_pop(popcnt-1,[(xoverpoint + 1):n]) = new_pop(popcnt,[(xoverpoint + 1):n]);
            new_pop(popcnt,[(xoverpoint + 1):n]) = tmp([(xoverpoint + 1):n]);
        end
    end         % chromosomes loop

    % POSTERIOR OPERATION
    % choose chromosomes without replacement - must do this to avoid reevaluation of empty clusters
    popmut = randperm(popul_size); popmut = popmut([1:num_mut]);
    for mutcnt = 1:num_mut
        % first we need to compute the mixture parameters (again ...)
        for mixcnt = 1:k
            ind = new_pop(popmut(mutcnt),:) == mixcnt;
            mix_propors(mixcnt) = sum(ind)/n;
            if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
                mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0;
            else
                mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
                % estimate the covariance matrix: if ill-conditioned, use the smoother
                mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),smth,1,1,sum(ind));
                % if still ill-conditioned, make all zeros
                if MatrixProblem(mix_covrmats(:,:,mixcnt)) ~= 0; mix_covrmats(:,:,mixcnt) = 0; end;
            end
        end     % mixtures loop

        for datcnt = 1:n
            % compute posterior probability for each cluster
            for mixcnt = 1:k
                if isequal(mix_covrmats(:,:,mixcnt),I)
                    % singleton, let pi be 0
                    posteriors(mixcnt) = 0;
                else
                    posteriors(mixcnt) = mix_propors(mixcnt)*mvnpdf(Data(datcnt,:),...
                        mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt));
                end
            end % mixtures loop
            posteriors = posteriors/sum(posteriors);            
            % assign datapoint to nearest cluster, if not all equally likely
            if (sum(posteriors == posteriors(1)) ~= k)
                [val,ind] = max(posteriors);
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
                        % estimate the covariance matrix: if ill-conditioned, use the smoother
                        % JAH corrected from mixcnt to val 20080831
                        mix_covrmats(:,:,val) = CovSmooth(Data(tmp,:),smth,1,1,sum(tmp));
                        % if still ill-conditioned, make all zeros
                        if MatrixProblem(mix_covrmats(:,:,val)) ~= 0; mix_covrmats(:,:,val) = 0; end;
                    end
                    tmp = (new_pop(popmut(mutcnt),:) == ind);
                    mix_propors(ind) = sum(tmp)/n;
                    if mix_propors(ind)*n <= 1
                        mix_covrmats(:,:,ind) = 0; mix_meanvecs(:,:,ind) = 0;
                    else
                        mix_meanvecs(:,:,ind) = mean(Data(tmp,:),1);
                        % estimate the covariance matrix: if ill-conditioned, use the smoother
                        % JAH corrected from mixcnt to ind 20080831
                        mix_covrmats(:,:,ind) = CovSmooth(Data(tmp,:),smth,1,1,sum(tmp));
                        % if still ill-conditioned, make all zeros
                        if MatrixProblem(mix_covrmats(:,:,ind)) ~= 0; mix_covrmats(:,:,ind) = 0; end;
                    end
                end
            end
        end     % datapoints loop
    end         % chromosomes to posterior operate

    % CONVEY BEST INDIVIDUAL INTO NEW POPULATION
    if elitism
        population = [new_pop;population(stdindex(1),:)];
    else
        population = new_pop;
    end
    popul_size = size(population,1);
    if rem(gencnt,2) == 1
        disp(sprintf('GEM Generation %0.0f of %0.0f: Best Score = %0.4f, %0.0f',gencnt,num_generns,best_score,termcount));
    end
end             % generations loop
% finish up with the progress plot
figure(fhga)
[AX,H1,H2] = plotyy(pltxyy(:,1),pltxyy(:,2),pltxyy(:,1),pltxyy(:,3)); xlabel('Generation');
title(['GEM(',init,') Progress: Objective function ',objec_func],'interpreter','none');
set(get(AX(1),'Ylabel'),'String','Maximum Value (o)','color','b');
set(H1,'color','b','marker','o'); set(AX(2),'ycolor','b');
set(get(AX(2),'Ylabel'),'String','Average Value (*)','color','r');
set(AX(2),'ycolor','r'); set(H2,'color','r','marker','*');
if not(isequal(savename,0))
    hgsave(fhga,[savename,'_',objec_func,'_GEM']); close(fhga), pause(1);
end

% from best chromosome, compute final posterior probabilities and parameters
posteriors = zeros(n,k);
for mixcnt = 1:k
    ind = (best_chrom == mixcnt);
    mix_propors(mixcnt) = sum(ind)/n;
    if sum(ind) <= 1    % ensure we store cov = zeros(p) for singleton or missing clusters
        mix_covrmats(:,:,mixcnt) = 0; mix_meanvecs(:,:,mixcnt) = 0;
    else
        mix_meanvecs(:,:,mixcnt) = mean(Data(ind,:),1);
        mix_covrmats(:,:,mixcnt) = CovSmooth(Data(ind,:),smth,1,1,sum(ind));    % JAH 20080831 was smth,1,1, but no good reason to always smooth here
    end
    if isequal(mix_covrmats(:,:,mixcnt),I) || (MatrixProblem(mix_covrmats(:,:,mixcnt)) ~= 0)
        % just let the posterior probability be 0
    else
        posteriors(:,mixcnt) = mix_propors(mixcnt)*mvnpdf(Data,mix_meanvecs(:,:,mixcnt),mix_covrmats(:,:,mixcnt));
    end    
end     % mixtures loop
posteriors = posteriors./repmat(sum(posteriors,2),1,k);   % posterior probability of group membership

%{
JAH 20061203, checked for octave 3.4.3

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