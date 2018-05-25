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

function plot_graphs(Y_all,imax,num_segments,num_segments_est,changepoints_real,changepoints_est,params,MU,THETA1,THETA2)
%This function is used to create two plots:
%       - The estimated segments of the synthetic data
%       - The original (solid line) and estimated (dashed line) of the channel parameters.

%% Plot data and segments
figure('Name','Segments');
plot(Y_all);
% real segment boundaries
line([changepoints_real' changepoints_real'], ylim,'Color','b');
% estimated segment boudnaries
line([changepoints_est' changepoints_est'], ylim,'Color','k','LineStyle','--');

legend('actual segmentation','output segmentation --')
xlabel('Time sample','fontsize',14)
ylabel('Received Channel Strength (dB)','fontsize',14)

grid on

%% Plot parameter values and segments
figure('Name','Parameters');
  
for i=1:num_segments_est(imax-1)
    line([changepoints_real(i),changepoints_real(i+1)],[params(1,i),params(1,i)],'Color','r','LineStyle','--')
    line([changepoints_real(i),changepoints_real(i+1)],[params(2,i),params(2,i)],'Color','g','LineStyle','--')
    line([changepoints_real(i),changepoints_real(i+1)],[params(3,i),params(3,i)],'Color','m','LineStyle','--')
    line([changepoints_est(i) changepoints_est(i)], ylim,'Color','k','LineStyle','--');
end

for i=1:num_segments
    line([changepoints_real(i),changepoints_real(i+1)],[MU(1,i),MU(1,i)],'Color','r')
    line([changepoints_real(i),changepoints_real(i+1)],[THETA1(1,i),THETA1(1,i)],'Color','g')
    line([changepoints_real(i),changepoints_real(i+1)],[THETA2(1,i),THETA2(1,i)],'Color','m')
end

ylabel('Estimated Values of ChannelParameters','fontsize',14)
legend('MU','THETA1','THETA2')

title('Original and Estimated Values of Channel Parameters')
grid on

end
