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

function [z]=viterbi_Par(K,P,params,Y_all,dij,a,num_nodes)

%Viterbi Algorithm for optimal sequence

%Input parameters
%   - K: State space dimension (candidate states are 1:K)
%   - P: transition matrix
%   - params: mean of each state
%   - Y_all: power measurments
%   - dij: the Euclidean distance between nodes i and j (see paper)
%   - a: vector of node distances from reference node (0,0)(see paper)
%   - num_nodes: the number of deployed nodes
%
% Output parameters
%   -z: Optimal State Sequence
% 
%%

T=length(Y_all);

T1=zeros(K,T);
T2=zeros(K,T);

%(observation likelihood of states at t=1?)
for i=1:K
    T1(i,1)=(K/T)*gau_m(params(:,i),Y_all(1,:),dij,a,num_nodes);
end

y_ti=zeros(K,1);
a1=zeros(K,1);
for i=2:T %for each value in the observation sequence
    
    %Pre calculate Gaussian density value (observation likelihood of state in current iteration
    for k=1:K
        y_ti(k)=gau_m(params(:,k),Y_all(i,:),dij,a,num_nodes);
    end
    for j=1:K %for every state
        T3=T1(:,i-1).*P(:,j).*y_ti;
        [a1(j,1),T2(j,i)]=max(T3);             
    end
   T1(:,i)=a1(:,1)/sum(a1);

end

z=zeros(T,1);
[~,z(T)]=max(T1(:,T));

for i=T:-1:2
    z(i-1)=T2(z(i),i);
end

end
