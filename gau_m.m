%% Version
% Last revision: May 2018 (Matlab R2017a)
% Author: Anastasios Dimas
%
%% Purpose
% The purpose of this code is to support the published paper:
% "Parameter Estimation for Hierarchical Channel Profiling", GlobalSIP, 2017 (pp. 234-238), IEEE,
% by Anastasios Dimas, Dionysios S. Kalogerias, Chryssalenia Koumpouzi, and Athina P. Petropulu.
%
% Any part of this code used in your work should cite the above publication.
%
% This code is provided "as is" to support the ideals of reproducible research. Any issues with this
% code should be reported by email to tasos.dimas@rutgers.edu. However, no guarantees are being made
% that the reported issues will be eventually fixed.
%
% The code is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
% available at https://creativecommons.org/licenses/by-nc-sa/4.0/
%
%%

function [v]=gau_m(params,y,dij,a,num_nodes)
%Gaussian density value
%
% Input parameters
%   - params: Vector of parameters
%   - y: sample measurments from all nodes at one time instant
%   - a: vector (see paper)
%   - dij: the Euclidean distance between nodes i and j (see paper)
%   - num_nodes: number of nodes deployed

% Output parameters
%   - v:
%
%%

% mean vector
m=a*params(1);
%covariance matrix
C_tq=(params(2)*exp(-dij/params(3)))+(params(4)*eye(num_nodes));

[~,p]=chol(C_tq);
if p>0 %check if covariance matrix is not positive definite
    C_tq=params(2)*exp(-dij)+params(4)*eye(num_nodes);
end

%multivariate normal dist.(calculation a little larger than function)
v = mvnpdf(y,m',C_tq);

end