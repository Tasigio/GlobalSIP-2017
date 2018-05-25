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

function [ mu_LS,theta1_LS,theta2_LS,sigmasq_LS ] = ParameterEstimation(Y_segs,r,T,a,dij )
%This function calculated the generates the positions at which the nodes are deployed.

% Input parameters
%   - Y_segs: measurments of considered segment
%   - z: sum of all samples of considered segment
%   - T: time duration of considered segment
%   - a: vector of node distances from reference node (0,0) (see paper)
%   - dij: distances between nodes i and j
%
% Output parameters
%   - mu_LS: least squares estiamte of parameter mu
%   - theta1_LS: least squares estiamte of parameter theta1
%   - theta2_LS: least squares estiamte of parameter theta2
%   - sigmasq_LS: least squares estiamte of parameter sigmasq

%%

% Least squares estimate of parameter mu
mu_LS= ((a'*a)\a')*r';

%Least square estimate of parameter sigmasq
sigmasq_LS=5;

% Joint least squares estimate of theta1 and theta2
H=((Y_segs'-mu_LS*a)*(Y_segs'-mu_LS*a)')-sigmasq_LS*eye(12);
H=H/T;

%Take off diagonal elements and remove any negative ones
% Htr=nonzeros(triu(H,0)); %this could remove elements that are zero in upper triagular
Htr=H(triu(true(size(H)),0));
ind=find(Htr<0);
Htr(ind)=[];

%Take the respective dij's of the previous
Ds=dij(triu(true(size(dij)),0));
Ds(ind)=[];

D=[ones(size(Ds,1),1) -Ds];
r=((D'*D)\D')*log10(Htr); 

theta1_LS=exp(r(1,1));
theta2_LS=1/r(2,1);

end
