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

function [ y ] = Synthetic_Data(num_nodes, dij,a, mu,theta_1,theta_2,sigmaSq,tau )
%This function generates the "original" power measurments of the nodes for a given state.
%Conditioned on this state, the measurments follow a Gaussian distribution (see paper)

% Input parameters
%   - num_nodes: the number of deployed nodes
%   - dij: the Euclidean distance between nodes i and j (see paper)
%   - a: vector of node distances from reference node (0,0)(see paper)
%   - mu:
%   - theta_1:
%   - theta_2:
%   - sigmaSq:
%   - tau: time duration of segment
%
% Output parameters
%   -y: generated power measurments

%%

%The covariance matrix
C_t=(theta_1*exp(-dij/theta_2))+(sigmaSq*eye(num_nodes));

%The multivariate Gaussian distribution
y = mvnrnd(a'*mu,C_t,tau);

end
