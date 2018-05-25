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

function [a,dij] = Create_Positions(axis_size,num_nodes)
%This function generates the positions at which the nodes are deployed.
%
% Input parameters
%   - axis_size: the length in meters of the considered square region
%   - num_nodes: the number of nodes we want to deploy (positive integer)
%
% Output parameters
%   -a: vector of node distances from reference node (0,0)(see paper)
%   -dij: the Euclidean distance between nodes i and j (see paper)
%
%%

%Integer spacing of nodes
pos_index=(1:axis_size)';

%All possible node coordinates in grid
grid_pos=zeros(axis_size*axis_size,2);
grid_pos(:,1)=kron(pos_index,ones(axis_size,1));
grid_pos(:,2)=kron(ones(axis_size,1),pos_index);

%Choose at random <num_nodes> positions from all possible node coordinates.
%This way the same coordinate can't be picked more than once.
node_pos_points=grid_pos(randperm(axis_size*axis_size,num_nodes),:);

%The distances between nodes i and j
dij=dist(node_pos_points');

%The distance of node i from reference node located at (0,0)
d=zeros(1,num_nodes);
for i=1:num_nodes
    d(i)=norm(node_pos_points(i,:)-[0,0]);
end

%Vector a
a=-10*log10(d(:));

end

