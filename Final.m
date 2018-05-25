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
close all;
clear;
clc;

%% General Parameters

%the length in meters of the considered square region
axis_size=10;
%maximum number of iterations
iter_max=7;
pvar=[300, 300 ,300, 300]; %sensitivity parameters
SigmaInv=diag(1./pvar);

%% Generate the Synthetic data

%Number of deployed nodes considered
num_nodes=12;
% Generate the positions at which the sensor nodes are deployed
[a,dij]=Create_Positions(axis_size,num_nodes);

%Number of original distinct segment considered
num_segments=10;
%The sample at which a new segment starts
changepoints_real=zeros(1,num_segments+1);
changepoints_real(1,1)=1;

%Array that holds the estimated number of segments in each iteration
num_segments_est=zeros(1,iter_max);
%Assume you start with segments of ~equal size
num_segments_est(1,1)=25;

%%

%For every segment, create the "original" synthetic data of power
%measurments and assign their "orginal" parameter values.
Y_all=[]; %Power measurments
MU=zeros(1,num_segments);
THETA1=zeros(1,num_segments);
THETA2=zeros(1,num_segments);
SIGMASQ=zeros(1,num_segments);

for q=1:num_segments
    %Original parameter values
    MU(1,q)=randi(10);
    THETA1(1,q)=1+randi(4);
    THETA2(1,q)=randi(3);
    SIGMASQ(1,q)=5;
    %segment will be divisible by initial estimate
    time_duration=(3000:num_segments_est(1,1):3500); 
    %assign random time duration of segment
    TAU(q,1)= time_duration(randi(length(time_duration)));
    
    [ y_tauq ] = Synthetic_Data(num_nodes,dij,a, MU(1,q),THETA1(1,q),THETA2(1,q),SIGMASQ(1,q),TAU(q,1) );
    
    Y_all=[Y_all; y_tauq];
    %time new segment starts
    changepoints_real(1,q+1)=changepoints_real(1,q)+TAU(q,1);
    
end

%Plot the generated data and segments
% plot(Y_all);
% line([changepoints_real' changepoints_real'], ylim,'Color','b');

%% Proposed approach

%Total duration of all segments
Sample_duration=length(Y_all(:,1));

%Divide duration of measurments into equal length segments
changepoints_est=(1:(Sample_duration/num_segments_est(1,1)):Sample_duration+1);

%Arrays for estimates of parameters
MU_EST=zeros(iter_max,num_segments_est(1,1));
THETA1_EST=zeros(iter_max,num_segments_est(1,1));
THETA2_EST=zeros(iter_max,num_segments_est(1,1));
SIGMASQ_EST=zeros(iter_max,num_segments_est(1,1));

E=zeros(1,iter_max); %Vector for escape probability

%Repeat iterative procedure for predefined number of rounds
for rounds=2:iter_max
    
    for q=1:num_segments_est(1,rounds-1) %for every segment
        tau=changepoints_est(q+1)-changepoints_est(q); %duration of segment
        Y_segs=Y_all(changepoints_est(q):changepoints_est(q+1)-1,:); %measurments of segment
        z=sum(Y_segs)/tau;
        [ mu_LS,theta1_LS,theta2_LS,sigmasq_LS ] = ParameterEstimation(Y_segs,z,tau,a,dij);
        
        MU_EST(rounds,q)= mu_LS;
        THETA1_EST(rounds,q)= theta1_LS;
        THETA2_EST(rounds,q)= theta2_LS;
        SIGMASQ_EST(rounds,q)= sigmasq_LS;
    end
    
    %Step 2 of iterative procedure
    %Escape probability
    E(1,rounds)=num_segments_est(1,rounds-1)/Sample_duration;
    
    %calculate g
    g=zeros(num_segments_est(1,rounds-1));
    P=zeros(num_segments_est(1,rounds-1));
    for i=1:num_segments_est(1,rounds-1) %for every segment
        s=zeros(num_segments_est(1,rounds-1),1);
        X=zeros(num_segments_est(1,rounds-1),4); %temporary vector of parameter differences
        for q=1:num_segments_est(1,rounds-1) %with every other segment
            X(q,:)=[MU_EST(rounds,i)-MU_EST(rounds,q),THETA1_EST(rounds,i)-THETA1_EST(rounds,q),THETA2_EST(rounds,i)-THETA2_EST(rounds,q),SIGMASQ_EST(rounds,i)-SIGMASQ_EST(rounds,q)];
            s(q)=X(q,:)*SigmaInv*X(q,:)';
        end
        
        %calculate denominator of g
        g_denom=sum(exp(-s));
        
        for q=1:num_segments_est(1,rounds-1)  %Calculate P
            g(i,q)=exp(-s(q))/g_denom;
            if i==q
                P(i,q)=1-E(1,rounds)+ (E(1,rounds)*g(i,q));
            else
                P(i,q)=E(1,rounds)*g(i,q);
            end
        end
    end
        
    params=[MU_EST(rounds,:);THETA1_EST(rounds,:);THETA2_EST(rounds,:);SIGMASQ_EST(rounds,:)];
  
    %Step 4 Viterbi----st=sequence of states
    fprintf('Viterbi round %d \n\n', rounds-1);
    st=viterbi_Par(num_segments_est(1,rounds-1),P,params,Y_all(:,:),dij,a,num_nodes);
    
    %New Segmentation
    %find how many different states there are
    u=length([1;find(diff(st)~=0)+1]);

    if rounds<iter_max
        changepoints_est= [1,(find(diff(st)~=0)+1)',Sample_duration+1];
        num_segments_est(rounds)=u;
    end    
    
end

plot_graphs(Y_all,iter_max,num_segments,num_segments_est,changepoints_real,changepoints_est,params,MU,THETA1,THETA2);
