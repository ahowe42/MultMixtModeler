 # M3 - Multivariate Mixture Modeler

This is updated Matlab code which I used to generate the mixture models and results in my doctoral thesis entitled **A New Generation of Mixture-Model Cluster Analysis with Information Complexity and the Genetic EM Algorithm**, as well as two publications:

-  J A Howe and H Bozdogan. Flexible Generalized Mixture Model Cluster Analysis with Elliptically-Contoured Distributions. *Journal of Pattern Recognition*, 1(1):5–22, May 2014

-  J A Howe and H Bozdogan.  Robust Mixture Model Cluster Analysis Using Adaptive Kernels. *Journal of Applied Statistics*, pages 320–336, November 2012.

Everything is run starting from the `Mixture_noGUI.m` script, so named because the version used to generate my results used a Matlab GUI.

## Model Parameters

There are 5 Matlab binary-formatted parameter files which can be used to setup and control the modeling simulations.  Any of them can be edited and saved to a new parameter file, so the modeling simulation can be run in different ways.

### `MixtureParams.mat`
This holds general modeling parameters as listed here:
- `convgcrit` - float; if the difference between successive iterations is less than this, an algorithm is deemed to have converged
- `elitism` - boolean; elitism causes, in the GA, the best solution for a generation to pass unchanged into the next generation
- `GAmax` - int; number of GA replications
- `htype` - NaN; placeholder
- `InfCrit` - cell of strings; placeholder
- `init_type` - string; type of initialization ('KM', 'GKM', 'GARM')
- `JustK` - boolean; instead of fitting mixture models for K = 2...Kmax, just fit the model with K=Kmax
- `Kmax` - int; largest mixture model to fit, depending on the value of JustK, will fit models for K = 2 ...Kmax
- `maxiter` - int; maximum number of iterations for K-means and EM algorithms
- `nochange_terminate` - int; number of GA iterations with no change in the objective function to allow early termination
- `num_generns` - int; maximum number of GA iterations
- `optim_func` - string; placeholder
- `optim_type` - string; type of optimization function to use ('EM', 'GEM')
- `popul_size` - int; number of individuals in a populaion in each generation of the GA
- `prob_mutate` - float; probability that any element in a GA solution will be mutated
- `prob_xover` - float; probability that a GA solution will undergo crossover
- `regul_func` - string; regularization function code for `CovSmooth`
- `regul_scale` - float; if `regul_func` = 'RDGREG', this is the regularization parameter
- `showplot` - boolean; optimization plots redrawn at each iteration

### `GaussianParams.mat`
This holds parameters specific to Gaussian mixture models (GMM)
- `InfCrit` - Cell with up to 10 string items specifying the names of the GMM information criteria function files
- `optim_func` - String specifying the name of the GMM optimization function (EM or GEM)
- `PStype` - String specifying the name of the type of mixture model being run ('Gaussian')

### `ECParams.mat`
This holds parameters specific to Elliptically-Contoured distribution mixture models (ECMM)
- `InfCrit` - Cell with up to 10 string items specifying the names of the ECMM information criteria function files
- `optim_func` - String specifying the name of the ECMM optimization function (EM or GEM)
- `ECtype` - String specifying the type of Elliptically-Contoured distribution to use (either 'KT' or PVII')
- `PStype` - String specifying the name of the type of mixture model being run ('EC')

### `KernelParams.mat`
This holds parameters specific to Kernel mixture models (KMM)
- `InfCrit` - Cell with up to 10 string items specifying the names of the KMM information criteria function files
- `optim_func` - String specifying the name of the KMM optimization function (EM or GEM)
- `htype` - Numeric code specifying the Kernel bandwidth estimator to use (see `MVKDE_Gauss.m`)
- `PStype` - String specifying the name of the type of mixture model being run ('Kernel')

### `PEKernParams.mat`
This holds parameters specific to Power Exponential Kernel mixture models (PEKMM)
- `InfCrit` - Cell with up to 10 string items specifying the names of the PEKMM information criteria function files
- `optim_func` - String specifying the name of the PEKMM optimization function (EM or GEM)
- `PStype` - String specifying the name of the type of mixture model being run ('PEKern')

## Advanced Analyses
My doctoral thesis included two advanced secondary techniques
- Subset Analysis 
- Influence Detection

The Matlab scripts that implemented these techniques, `Mixture_AllSubsAnal`, `Mixture_GASubsAnal`, and `Mixture_OutDetect` have not yet been converted to work without the previous GUI, and so aren't included in the repository.
