function [C, L, U] = clusterRDM_spectralClustering(W, maxclust, varargin)
%% 
% 
% 
%   C: An Nx1 vector of cluster number.
%   L: (normalized) Graph Laplacian, 
%   U: Laplacian eigen vectors
% 
% 
% Mostly relies on Ingo Buerk's code (2011,2012).
% (https://mathworks.com/matlabcentral/fileexchange/34412-fast-and-efficient-spectral-clustering)
% 
% 
% ayokoi (Jul 2018)


normalization = 'JordanWeiss'; % 'ShiMalik', 'non'
maxiter = 1000; % for kmeans option
replicates = 500; % for kmeans option
vararginoptions(varargin, {'normalization','maxiter','replicates'});

% init rng with fixed seed for replicability (kmeans inside spectralclustering relies on rng)
randSeed = 1234;
method = 'twister';
rng(randSeed, method); 
switch lower(normalization)
    case 'non'
        Type = 1;
    case 'shimalik';
        Type = 2;
    case 'jordanweiss'
        Type = 3;
    otherwise
        warning('Unsupported method. Use JordanWeiss, instead.');
        Type = 3;
end

% Call spectral clustering (inline function)
[C, L, U] = SpectralClustering(W, maxclust, Type, maxiter, replicates); 

% re-assgn cluster values if C is sparse matrix
% if nargout>3
%     C = full(Csp); C_=zeros(size(C,1),1);
%     for i=1:size(C,1)
%         C_(i,1) = find(C(i,:)==1);
%     end
%     C=C_;
% end

end

% Inline function by Ingo Buerk with slight modification
function [C, L, U] = SpectralClustering(W, k, Type, Maxiter, Replicates)
%SPECTRALCLUSTERING Executes spectral clustering algorithm
%   Executes the spectral clustering algorithm defined by
%   Type on the adjacency matrix W and returns the k cluster
%   indicator vectors as columns in C.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%   'Type' - Defines the type of spectral clustering algorithm
%            that should be used. Choices are:
%      1 - Unnormalized
%      2 - Normalized according to Shi and Malik (2000)
%      3 - Normalized according to Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
%     Statistics and Computing 17 (4), 2007
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% calculate degree matrix
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
L = D - W;

% compute normalized Laplacian if needed
switch Type
    case 2
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate inverse of D
        D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L;
    case 3
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate D^(-1/2)
        D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L * D;
end

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
dif   = eps;
[U, lambda, flag] = eigs(L, k, dif);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
% C = kmeans(U, k, 'start', 'cluster', ...
%                  'EmptyAction', 'singleton');
Opt = statset('UseParallel',0,'Maxiter', Maxiter);
C = kmeans(U, k, 'start', 'cluster', ...
    'EmptyAction', 'drop',...
    'Replicates', Replicates,...
    'Options', Opt);

% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
%C = sparse(1:size(D, 1), C, 1);

end