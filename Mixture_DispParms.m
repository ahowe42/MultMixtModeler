% usage: Mixture_DispParms
% Use this script to display in a nice formatted "table" all the current
% parameters select in M3.

% display the parameters for the mixture modeling
disp(['Initialization Type: ',MMM.init_type])
disp(['Maximization Type: ',MMM.optim_type])
disp(sprintf('Evalute K in range: [%0.0f,%0.0f]',min(MMM.Krange),MMM.Kmax))
disp(sprintf('Population Size: %0.0f\nGeneration Size: %0.0f\nMininum # of Generations: %0.0f',MMM.popul_size,MMM.num_generns,MMM.nochange_terminate))
disp(sprintf('Crossover Rate: %0.2f\nMutation Rate: %0.2f',MMM.prob_xover,MMM.prob_mutate))
if MMM.elitism; disp('Elitism is: ON'); else; disp('Elitism is: OFF'); end;
disp(['Covariance Smoothing Code: ',MMM.regul_func])
disp('Information Criteria Functions:')
for infcnt = 1:10
    if not(isempty(MMM.InfCrit{infcnt})); disp(['     ',MMM.InfCrit{infcnt}]), end;
end
disp(['Maximization Function: ',MMM.optim_func])
if isequal(MMM.init_type,'GARM')
    disp(sprintf('Mahalanobis Scale Factor: %0.0f',MMM.regul_scale))
end
if isequal(MMM.optim_type,'EM')
    disp(sprintf('EM Convergence Criteria: %0.10f\t Maximum # of Iterations: %0.0f',MMM.convgcrit,MMM.maxiter))
end
switch MMM.data_type
    case -1     % simulated data
        disp(sprintf('%d observations of Simulated Data: %s',n,MMM.data_file))
    case 0      % data of known structure (with labels)
        disp(['Known Data: ',MMM.data_file])        
    case 1      % data of unknown structure
        disp(['Unknown Data: ',MMM.data_file])
end
disp(['Problem Specific - ',MMM.PStype])
switch MMM.PStype
    case {'Gaussian','PEKern'}
        % none
    case 'Kernel'
        disp(sprintf('Kernel Bandwidth Matrix Estimator Type: %d',MMM.htype))
    case 'EC'
        disp(['Elliptically-Contoured Distribution Subtype: ',MMM.htype])
end

% JAH 20070109, adapted for octave 3.4.3 20120319
% this code may be freely used, modified, and distributed (at no charge)
% as long as this footer remains unaltered
