%-----------------------------------------------------------------
%  SIMPLEREADTRIALS  - Analysis of simulation results
%  [D,M]=readtrials(file_trials, nb_trials, t_end)
%   
% read mean activity of all orientation columns in all trials
% compute mean, var, cov
% returns D: matrix of all trials 
%         M: vector of mean activity
%          
%
% to plot mean activity for trials 5 plot(D(5,:))
%------------------------------------------------------------------

function [D,M]=readtrials(file_trials, nb_trials, t_end)

Nx=1008;

% make it look nice
%set(0,'DefaultLineLineWidth',1.5)
%set(0, 'DefaultAxesFontName', 'Helvetica')
%set(0, 'DefaultAxesFontSize', 12)
%set(0, 'DefaultAxesFontWeight','bold')
%set(0,'DefaultLineMarkerSize',4)
%set(0,'DefaultLineMarkerFaceColor', [0 0 0])
Nori = 1./Nx:180./Nx:180; % orientation preference of all columns

%#################################
% read the results of all trials
%#################################

fid=fopen(file_trials);
fprintf('...  reading responses for all %d trials (t_end = %f) in file %s\n', nb_trials, t_end,file_trials);
[D,COUNT]=fscanf(fid, '%f', [nb_trials, Nx]);

%figure;
%subplot(3,1,1);
%contour(D);
%colorbar;
%ylabel('# trials');

%############################################
% compute and plot response mean and variance
%############################################

fprintf('Mean ...\n');
M=mean(D);
V=(1008/Nx)*var(D);

%figure;
%plot(Nori,M,'-', Nori,V,'-.');
%axis([1 180 0 (2*max(M))]);
%legend('mean','var');
%xlabel('orientation');
