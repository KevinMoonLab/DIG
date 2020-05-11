% MAIN CODE MUL. Time series Manifold learning

function [Y, info] = DIG(data, varargin)
% DIG Runs DIG for visualizing noisy dynamical systems
%   Y = DIG(data) runs DIG on data (rows: Time series, columns: Time order data)
%   with default parameter settings and returns a 2 dimensional embedding.
%
%   If data is sparse PCA without mean centering will be done to maintain
%   low memory footprint. If data is dense then normal PCA (with mean
%   centering) is done.
%
%   Y = DIG(data, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to
%   specify optional parameter name/value pairs that control further details
%   of DIG.  Parameters are:
%
%   'ndim' - number of (output) embedding dimensions. Common values are 2
%   or 3. Defaults to 2.
%
%   'k' - number of nearest neighbors for bandwidth of adaptive alpha
%   decaying kernel or, when a=[], number of nearest neighbors of the knn
%   graph. For the unweighted kernel we recommend k to be a bit larger,
%   e.g. 10 or 15. Defaults to 5.
%
%   'a' - alpha of alpha decaying kernel. when a=[] knn (unweighted) kernel
%   is used. Defaults to 40.
%
%   't' - number of diffusion steps. Defaults to [] wich autmatically picks
%   the optimal t.
%
%   't_max' - maximum t for finding optimal t. if t = [] optimal t will be
%   computed by computing Von Neumann Entropy for each t <= t_max and then
%   picking the kneepoint. Defaults to 100.
%
%   'npca' - number of pca components for computing distances. Defaults to
%   100.
%
%   'mds_method' - method of multidimensional scaling. Choices are:
%
%       'mmds' - metric MDS (default)
%       'cmds' - classical MDS
%       'nmmds' - non-metric MDS
%
%   'distfun' - distance function. Default is 'euclidean'.
%
%   'distfun_mds' - distance function for MDS. Default is 'euclidean'.
%
%   'pot_method' - method of computing the PHATE potential dstance. Choices
%   are:
%
%       'log' - -log(P + eps). (default)
%
%       'sqrt' - sqrt(P). (not default but often produces superior
%       embeddings)
%
%       'Gamma' - 2/(1-\gamma)*P^((1-\gamma)/2)
%
%   'Gamma' - gamma value for gamma potential method. Value between -1 and
%   1. -1 is diffusion distance. 1 is log potential. 0 is sqrt. Smaller
%   gamma is a more locally sensitive embedding whereas larger gamma is a
%   more globally sensitive embedding. Defaults to 0.5.


set_parameters;

% Distances among time-windows of data
disp('Computing distances among time-windows of data')
switch series_distance
    case 'InfGeo'
      D = InfGeoDis(data, L);
      disp('Using information geometry distance, make sure the data is approximatly Mul. Gaussian with 0 mean')
    case 'EIG'
      [D,  info] = EIGDis(data, L, ndim, nbins, L2);
      k = 20; 
    case 'KL'
      disp('Using Kullback-Liebler divergence, make sure the data is approximatly Mul. Gaussian with 0 mean')
      D = KLDis(data,L);
end 

info.D = D;
% affinity matrix      
disp 'Computing Affinity Matrix'
Pmm = affinity_matrix(D,k,a, adaptativeK);
info.affinity_m = Pmm;

% Diffusion
disp 'Diffusing operator'
P_t = Pmm^t;

% Distances among rows in the of P_t     
disp 'Computing potential/gamma/inf_geometry distances'
PDX = P_distances(P_t, pot_method, Gamma);
info.Potencial_D = PDX;
disp 'MDS'
Y = multi_scale(PDX, ndim, series_distance);
 
end 