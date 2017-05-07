%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for multisensory decision making
% Modified by HH 2017 @ UNIGE
% Adapted for the vestibular-visual multisensory heading discrimintation task
%
% lip_HH(para_override)
% para_override: {'name1',value1; 'name2',value2, ...}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lip_HH(para_override)

%clear
%clear global

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end
addpath(genpath(pwd));

rand('state',sum(100*clock));

if strcmp(version('-release'),'2014b')    % ION cluster;
    hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name
    if isempty(gcp('nocreate'))
        parpool(hostname(1:strfind(hostname,'.')-1),20);
    end
    ION_cluster = 1;
else % My own computer
    if strcmp(version('-release'),'2013a')
        if matlabpool('size') == 0
            matlabpool;
        end
    else
        if isempty(gcp('nocreate'))
            parpool;
        end
    end
    ION_cluster = 0;
end

% === Switches ===
if_debug = ~ION_cluster;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============ Sizes ============
% Input layers
N_targets = 2; % Target input
N_vis = 300; % Visual motion signal
N_vest = 300; % Vestibular motion signal

% Today we decide to add a perfect integrator layer here. HH20170317 in UNIGE
% This layer has a very long time constant and feeds it's input to LIP, whereas LIP,
% which has a short time constant, keeps receiving target input but does not integrate it.
% For most part of the code, I just rename the old "_lip" stuff to "_int"
N_int = 300;
N_lip = 300;  % Output layers (Real LIP layer)

% ============ Times ============
if ION_cluster
    dt = 5e-3; % Size of time bin in seconds
else 
    dt = 5e-3;
end
trial_dur_total = 1.7; % in s
stim_on_time = 0.2; % in s
motion_duration = 1.5; % in s

% ============ Decision bound ============
if_bounded = 1; % if_bounded = 1 means that the trial stops at the bound (reaction time version)
% f_bound = @(x) max(x);  % A bug: if lip_HH is a function, this anonymous function cannot be broadcasted into parfor loop
%  f_bound = @(x) abs(x(right_targ_ind)-x(left_targ_ind));

% Smoothed max
decis_thres = 34*[1 1 1+3/34]; % bound height, for different conditions
% Smoothed diff
% decis_thres = 13*[1 1 1+2/13]; % bound height, for different conditions

att_gain_stim_after_hit_bound = [0 0 0];

% =============== Conditions ===============
if ION_cluster
    unique_heading = [-8 -4 -2 -1 0 1 2 4 8];
    unique_stim_type = [1 2 3];
    N_rep = 100; % For each condition
else
    unique_heading = [-8 0 8];
    unique_stim_type = [1 2 3];
    N_rep = 5;
end

% =================== Stimuli ===================
% Temporal dynamics
num_of_sigma = 3.5; amp = 0.2;

% Parameters roughly follow the Yong Gu's MST data.
% Visual
coherence = 12;
r_spont_vis = 10;
b_pref_vis = 1.7;
b_null_vis = -0.2;
K_vis = 1.5;
K_cov_vis = 2;
var_vis = 1e-5;

% Vestibular
equivalent_conherence = coherence;
r_spont_vest = r_spont_vis;
b_pref_vest = b_pref_vis;
b_null_vest = b_null_vis;
K_vest = K_vis;
K_cov_vest = K_cov_vis;
var_vest = var_vis;

% Parameters for MT from Mazurek and Shadlen, except width.
% K_cov_mt and var_mt controls the correlations.
%   cov_mt(j,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(j))/360 *2*pi)-1));
%     w_mt = real(sqrtm(cov_mt));
%     aux_proba_in = proba_in + w_mt*randn(N_vis,1);
% r_spont_vis = 20;
% b_pref_vis = 0.4;
% b_null_vis = -0.2;
% K_vis = 4;
% K_cov_vis = 2;
% var_vis = 1e-5;

% Parameter of the activity evoked by the visual targets
max_rate_target = 30; %% in Hz
b_target = 0;
K_target = 4;
slope_norm_target=0.5;
dc_norm_target = 0.5;

% Gains
% Maybe should be set such that gain_vis:gain_vest = (integral of abs(a))/(integral of v)
gain_vel_vis = 7; % (m/s)^-1
gain_acc_vest = 2; %  gain_vel_vis * sum(vel)/sum(abs(acc)); % (m^2/s)^-1

% =================== Network configuration ===================
% -- Time constant for integration
time_const_int = 10000e-3; % in s
time_const_lip = 100e-3; % in s

% ---- Visual to INTEGRATOR ----
g_w_int_vis = 10;
dc_w_int_vis = 0;
k_int_vis = 5; % Larger, narrower
k_int_vis_along_vis = 0.1; % Larger, wider

% ---- Vestibular to INTEGRATOR ----
g_w_int_vest = g_w_int_vis;
dc_w_int_vest = dc_w_int_vis;
k_int_vest = k_int_vis;
k_int_vest_along_vest = k_int_vis_along_vis;

% ----- Targets to LIP ------
g_w_lip_targ= 8;
k_lip_targ= 5;
att_gain_targ = 1; % Drop in attention to visual target once motion stimulus appears.

% % ------- Recurrent connectivity in INTEGRATOR --------
% g_w_int_int = 35;
% K_int_int = 10;
% dc_w_int_int = -11;
%
% amp_I_int = 0;  % Mexico hat shape
% K_int_I = 2;
% bias_int = 0;
% Input-output function of INTEGRATOR
threshold_int = 0.0;

% ----- INTEGRATOR to the real LIP ----
g_w_lip_int = 20;
k_lip_int = 5;
dc_w_lip_int = -7;

amp_I_lip_int = 0;  % Mexico hat shape
k_lip_int_I = 2;

% ----- LIP recurrent connection ------
g_w_lip_lip = 5;
k_lip_lip = 5;
dc_w_lip_lip = -3;

amp_I_lip = 0;  % Mexico hat shape
k_lip_I = 2;
bias_lip = 0;

% Input-output function of LIP
threshold_lip = 0.0;

%%%%%%%%%%%%%%%% Override by input argument %%%%%%%%%%%%%%%
if nargin == 1
    len = size(para_override,1);
    for ll = 1:len
        if exist(para_override{ll,1},'var')
            eval([para_override{ll,1} '= para_override{ll,2};']);
            fprintf('Overriding %s = %s...\n',para_override{ll,1},num2str(para_override{ll,2}));
        else
            fprintf('Parameter ''%s'' not found...\n',para_override{ll,1});
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pure Parameters End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some preparation (after the pure parameter session in case some paras are overriden)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Network related ---- 
prefs_vis = linspace(-180,180,N_vis);
prefs_vest = linspace(-180,180,N_vis);
prefs_int = linspace(-180,180,N_int);
prefs_lip = linspace(-180,180,N_lip);
[~,right_targ_ind] = min(abs(prefs_lip - 90)); % Left and right for the heading task
[~,left_targ_ind] = min(abs(prefs_lip - -90));

gain_func_along_vis = @(x) (1./(1+exp(-(abs(90-abs(90-abs(x))))/k_int_vis_along_vis))-0.5)*2;
gain_func_along_vest = @(x) (1./(1+exp(-(abs(90-abs(90-abs(x))))/k_int_vest_along_vest))-0.5)*2;

% ---- Time related ---- 
trial_dur_total_in_bins = round(trial_dur_total/dt); %% in time bins (including pre_trial_dur)
stim_on_time_in_bins = stim_on_time/dt; % Time of motion start, in time bins.

if stim_on_time_in_bins>trial_dur_total_in_bins
    error('pre_trial_dur cannot be longer than trial_dur');
end

% --- Stimuli related ---
t_motion = 0:dt:trial_dur_total-stim_on_time; % in s
miu = motion_duration/2; sigma = motion_duration/2/num_of_sigma;
vel = exp(-(t_motion-miu).^2/(2*sigma^2));
vel = vel*amp/sum(vel*dt) ; % in m/s. Normalized to distance
acc = diff(vel)/dt; % in m/s^2

% To make sure t_motion, vel, and acc have the same length
t_motion(end) = [];
vel(end) = [];

if if_debug
    figure(111);clf
    set(gcf,'name','Motion profile','uni','norm','pos',[0.632       0.381       0.358       0.403]);
    h = plotyy(t_motion,vel,t_motion,acc);
    ylabel(h(1),'Velocity (m/s)');
    ylabel(h(2),'Acceleration (m^2/s)');
end

%%%%%%%%%%% Pack all paras %%%%%%%%%%%
para_list = who;
for i = 1:length(para_list)
    paras.(para_list{i}) = eval(para_list{i});
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization (I moved here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I reorganized the variable naming and storage conventions. HH20170302
% All data have the size of (neurons, time, trial, stimuli){train/test}
% Three types of data are used:
% 1. firing rate in Hz (rate_xxx),
% 2. firing probability for each bin (proba_xxx),
% 3. binary spike train for each bin (spikes_xxx)
% The network dynamics are computed like: Rate_in --> Prob_in --> Spike_in --> Spikes_out --> Rate_out

aux1_rate_int = zeros(N_int,trial_dur_total_in_bins+2);
aux2_rate_int = zeros(N_int,trial_dur_total_in_bins+2);

% spikes_target = cell(2,1); spikes_vis = cell(2,1); spikes_vest = cell(2,1);
% rate_int = cell(2,1); spikes_int = cell(2,1);
% rate_lip = cell(2,1); spikes_lip = cell(2,1);

% spikes_count_vis = cell(2,1); spikes_count_vest = cell(2,1); proba_count_int = cell(2,1); spikes_count_int = cell(2,1);

% for m = 1:if_test_set+1 % Train and Test

%     % Raw data
%     spikes_target = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%     spikes_vis = zeros(N_vis,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%     spikes_vest = zeros(N_vest,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%
%     rate_int = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%     spikes_int = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%     rate_lip = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%     spikes_lip = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));

% Spike count in sliding windows for Fisher information
%     spikes_count_vis = nan(N_vis, fisher_N_window, N_trial,length(unique_heading),length(unique_stim_type));
%     spikes_count_vest = nan(N_vest, fisher_N_window, N_trial,length(unique_heading),length(unique_stim_type));
%     proba_count_int = nan(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_stim_type));
%     spikes_count_int = nan(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_stim_type));
%     if shuffle_flag==1
%         spikes_count_int_shuff = zeros(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_stim_type));
%     end

% end

%  Network connections
w_int_vis = zeros(N_int,N_vis);
w_int_vest = zeros(N_int,N_vest);
% w_int_targ = zeros(N_int,N_int); % Not N_target, but N_int (supposing the same number as the lip neurons)
w_lip_targ = zeros(N_lip,N_lip); % Now integration layer does not receive target input, but LIP layer does. HH20170317
w_int_int = zeros(N_int,N_int);
w_lip_int = zeros(N_lip,N_int);
w_lip_lip = zeros(N_lip,N_lip);


for nn=1:N_int
    
    % -- VIS->Int, Vest->Int --
    w_int_vis(nn,:) = g_w_int_vis/N_vis *(exp(k_int_vis * (cos((prefs_int(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vis/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        + dc_w_int_vis/N_vis;   % Offset
    w_int_vest(nn,:) = g_w_int_vest/N_vest *(exp(k_int_vest * (cos((prefs_int(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vest/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        + dc_w_int_vest/N_vest;   % Offset
    
    % Added a K_int_vis_sin factor to tweak the slices of weight matrix along the vis/vest axis (sigmoid --> step)
    %     w_int_vis(nn,:) = g_w_int_vis/N_vis *(exp(k_int_vis * (cos((prefs_int(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
%         .* gain_func_along_vis(prefs_vis)... % Gaussian(theta_int - +/-90) * Sin(heading)
%         + dc_w_int_vis/N_vis;   % Offset
%     w_int_vest(nn,:) = g_w_int_vest/N_vest *(exp(k_int_vest * (cos((prefs_int(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
%         .* gain_func_along_vest(prefs_vest) ... % Gaussian(theta_int - +/-90) * Sin(heading)
%         + dc_w_int_vest/N_vest;   % Offset
    
    % -- Int->Int --
    %     w_int_int(nn,:) = g_w_int_int/N_int*...   %  Integrator recurrent
    %         ((exp(K_int_int*(cos((prefs_int-prefs_int(nn))/360*2*pi)-1)))-...
    %         amp_I_int*(exp(K_int_I*(cos((prefs_int-prefs_int(nn))/360*2*pi)-1))))...
    %         + dc_w_int_int/N_int;
end

for nn=1:N_lip
    
    % -- Targ->LIP --
    w_lip_targ(nn,:) = g_w_lip_targ/N_int *(exp(k_lip_targ*(cos((prefs_lip-prefs_lip(nn))/360 *2*pi)-1)));  %  Target input
    
    % -- Int->LIP --
    w_lip_int(nn,:) = g_w_lip_int/N_int*...
        ((exp(k_lip_int*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1)))-...
        amp_I_lip_int*(exp(k_lip_int_I*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1))))...
        + dc_w_lip_int/N_int;
    
    % -- LIP->LIP --
    w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...
        ((exp(k_lip_lip*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
        amp_I_lip*(exp(k_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
        + dc_w_lip_lip/N_lip;
end

if if_debug
    figure(90); clf;
    set(gcf,'uni','norm','pos',[0.006       0.099       0.779       0.766],'name','Weights & Raster plot');
    
    subplot(4,3,1);
    imagesc(prefs_lip,prefs_lip,w_lip_lip'); colorbar; axis  xy; title('LIP->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{lip}');
    set(gca,'xtick',-180:90:180);
    hold on; plot(prefs_lip,w_lip_lip(:,end/2)/range(w_lip_lip(:,end/2))*100,'linew',3,'color','c');
    plot(xlim,[0 0],'--c');
    
    subplot(4,3,4);
    imagesc(prefs_lip,prefs_int,w_lip_int'); colorbar; axis  xy; title('Int->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180);
    hold on; plot(prefs_lip,w_lip_int(:,end/2)/range(w_lip_int(:,end/2))*100,'linew',3,'color','c');
    plot(xlim,[0 0],'--c');
    
    subplot(4,3,7);
    % surf(prefs_int,prefs_int,w_int_int');
    imagesc(prefs_int,prefs_int,w_int_int');    colorbar; axis xy; title('Int->Int');
    xlabel('\theta_{int}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    hold on; plot(prefs_int,w_int_int(:,end/2)/range(w_int_int(:,end/2))*100,'linew',3,'color','c');
    plot(xlim,[0 0],'--c');
    
    %     subplot(4,3,7);
    %     imagesc(prefs_int,prefs_vest,w_int_vest'); colorbar; axis  xy; title('Vest->Int');
    %     xlabel('\theta_{int}'); ylabel('\theta_{vest}'); ylim([-20 20]);
    %     set(gca,'xtick',-180:90:180);
    
    subplot(4,3,10);
    imagesc(prefs_int,prefs_vis,w_int_vis'); colorbar; axis  xy; title('Vis/Vest->Int');
    xlabel('\theta_{int}'); ylim([-20 20]);
    set(gca,'xtick',-180:90:180);
    hold on; plot(prefs_int,w_int_vis(:,end/2)/max(w_int_vis(:,end/2))*20,'linew',3,'color','c');
    plot(xlim,[0 0],'--c');
    
    temp_ax = axes('Pos',[0.05 0.124 0.053 0.134]);
    plot(prefs_vis,gain_func_along_vis(prefs_vis)); hold on;
    plot([-20 -20],ylim,'r-'); plot([20 20],ylim,'r-');
    view(270,90);     set(gca,'xtick',-180:90:180);
    xlabel('\theta_{vis/vest}');
    
    colormap hot;
end

%  Sensory correlation matrix
[xx,yy] = meshgrid(prefs_vis,prefs_vis);
cov_vis = var_vis*exp(K_cov_vis*(cos((xx-yy)/360 *2*pi)-1));

[xx,yy] = meshgrid(prefs_vest,prefs_vest);
cov_vest = var_vest*exp(K_cov_vest * (cos((xx-yy)/360 *2*pi)-1));

w_cov_vis = real(sqrtm(cov_vis));
w_cov_vest = real(sqrtm(cov_vest));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_time = tic;

% === To reduce the overheads, I put all conditions*heading*reps into a large matrix. HH20170417 ===
[x,y]=meshgrid(1:length(unique_stim_type),1:length(unique_heading));
z=[x(:) y(:)];
ss_hh_cheatsheet = reshape(repmat(z,1,N_rep)',2,[])'; % All combinations of ss and hh
n_parfor_loops = size(ss_hh_cheatsheet,1);

proba_vis_for_each_heading = zeros(N_vis,length(unique_heading));
proba_vest_for_each_heading = zeros(N_vest,length(unique_heading));

% -- Data to save in parloop --
% % Raw data
% spikes_target = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
% spikes_vis = zeros(N_vis,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
% spikes_vest = zeros(N_vest,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
%
% rate_int = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
% spikes_int = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
% rate_lip = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));
% spikes_lip = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_stim_type));

spikes_target = zeros(N_lip,trial_dur_total_in_bins,n_parfor_loops);
spikes_vis = zeros(N_vis,trial_dur_total_in_bins,n_parfor_loops);
spikes_vest = zeros(N_vest,trial_dur_total_in_bins,n_parfor_loops);

rate_int = zeros(N_int,trial_dur_total_in_bins,n_parfor_loops);
spikes_int = zeros(N_int,trial_dur_total_in_bins,n_parfor_loops);
rate_lip = zeros(N_lip,trial_dur_total_in_bins,n_parfor_loops);
spikes_lip = zeros(N_lip,trial_dur_total_in_bins,n_parfor_loops);

RT = (trial_dur_total-stim_on_time) * ones(n_parfor_loops,1);

% --- Generate stuffs that are ss/hh- dependent --
for hh = 1:length(unique_heading)  % Motion directions
    
    % ==== Not trial- or time- dependent stuffs ===
    %         fprintf('cond = %g, heading = %g\n  ',cond_this, stim_this);
    
    % -- Sensory responses --
    % proba_vis: proba of firing of visual neurons in response to motion
    max_rate_vis = r_spont_vis + b_pref_vis * coherence;
    b_vis = r_spont_vis + b_null_vis * coherence;
    proba_vis_for_each_heading(:,hh) = ((max_rate_vis-b_vis)*exp(K_vis*(cos((prefs_vis'-unique_heading(hh))/360*2*pi)-1))+b_vis)*dt;
    
    %     max_rate_vest = r_spont_vest + b_pref_vis * coherence;
    %     b_vis = r_spont_vis + b_null_vis * coherence;
    
    % Vestibular response does not depend on coherence.
    % Here I just set the vestibular activity similar to visual response under 'equivalent_coh' coh
    max_rate_vest = r_spont_vest + b_pref_vest * equivalent_conherence;
    b_vest = r_spont_vest + b_null_vest * equivalent_conherence;
    proba_vest_for_each_heading(:,hh) = ((max_rate_vest-b_vest)*exp(K_vest*(cos((prefs_vest'-unique_heading(hh))/360*2*pi)-1))+b_vest)*dt;
    
    %         if if_debug
    %             figure(91);clf;
    %             plot(prefs_vis,proba_vis/dt,'r');
    %             plot(prefs_vest,proba_vest/dt,'b');
    %             title('Sensory population tuning');
    %         end
end

for ss = 1:length(unique_stim_type)
    % === Some slicing stuffs necessary for parfor ===
    att_gain_this = att_gain_stim_after_hit_bound(unique_stim_type(ss));
    decis_thres_this = decis_thres(unique_stim_type(ss));
end

% --- Generate other stuffs that are not ss/hh- dependent ---
% proba_targ: proba of firing of target neurons in response to visual targets
% pos_targ = stim_this + [0:360/N_targets:359.9];
pos_targ = [-90 90]; % Always these two for fine discrimination task. HH2017
proba_target_tuning = zeros(N_lip,1);
for nn=1:N_targets
    proba_target_tuning = proba_target_tuning + ((max_rate_target-b_target)*...
        exp(K_target*(cos((prefs_lip'-pos_targ(nn))/360*2*pi)-1))+b_target)*dt;
end
proba_target_tuning = proba_target_tuning/(slope_norm_target*(N_targets/2)+dc_norm_target);   %  Divisive normalization

aux_proba_target = proba_target_tuning*ones(1,trial_dur_total_in_bins); % Expand along the time axis

aux_proba_target(:,1:round(stim_on_time/dt)) = 1.7*aux_proba_target(:,1:round(stim_on_time/dt)); % Saliency of target

%% == Now I only have to call parfor for just one time. ==
parfor_progress(n_parfor_loops);

parfor tt = 1:n_parfor_loops % For each trial
    
    % -- Get data from cheat sheet --
    cond_this = ss_hh_cheatsheet(tt,:);
    ss_this = cond_this(1);
    hh_this = cond_this(2);
    
    % -- Target input spike train --
    spikes_target_this = rand(N_lip,trial_dur_total_in_bins)<(aux_proba_target);
    
    % -- Visual input spike train --
    aux_proba_vis = proba_vis_for_each_heading(:,hh_this)*[zeros(1,stim_on_time_in_bins) vel*gain_vel_vis]...
        + w_cov_vis*randn(N_vis,trial_dur_total_in_bins) ...
        .*repmat([zeros(1,stim_on_time_in_bins) ones(size(vel))],N_vis,1);
    
    % -- Vestibular ---
    % With the temporal gain of abs(acc)
    aux_proba_vest = proba_vest_for_each_heading(:,hh_this)*[zeros(1,stim_on_time_in_bins) abs(acc)*gain_acc_vest]...
        + w_cov_vest*randn(N_vest,trial_dur_total_in_bins) ...
        .*repmat([zeros(1,stim_on_time_in_bins) ones(size(vel))],N_vis,1);
    
    % -- Stimulus condition selection --
    if unique_stim_type(ss_this) == 1
        aux_proba_vis = 0*aux_proba_vis; % Shut down visual activity
    elseif unique_stim_type(ss_this) ==2
        aux_proba_vest = 0*aux_proba_vest; % Shut down vestibular activity
    end
    
    spikes_vis_this = rand(N_vis,trial_dur_total_in_bins)<(aux_proba_vis);
    spikes_vest_this = rand(N_vest,trial_dur_total_in_bins)<(aux_proba_vest);
    
    % -- Reset the attention for motion stimuli --
    att_gain_stim = 1; % Drop in attention to motion stimuli due to hitting whatever decision bound.
    
    % == Real network dynamics after stimlus onset ==
    %                 k = stim_on_time_in_bins; % Time begins at stim_on_time
    k = 1;  % Time begins at the very beginning
    
    % -- Variables for this parfor loop --
    rate_int_this = zeros(N_int,trial_dur_total_in_bins);
    rate_lip_this = zeros(N_lip,trial_dur_total_in_bins);
    spikes_int_this = zeros(N_int,trial_dur_total_in_bins);
    spikes_lip_this = zeros(N_lip,trial_dur_total_in_bins);
    decision_ac = zeros(N_int,trial_dur_total_in_bins+2);
    
    while k<=trial_dur_total_in_bins-1
        
        % -- Update INTEGRATOR layer --
        %                     rate_int(:,k+1,tt,hh,ss) = bias_int + (1-dt/time_const_int)*rate_int(:,k,tt,hh,ss)...   %  Self dynamics.  in Hz!
        %                         + 1/time_const_int * (...
        %                               w_int_int * spikes_int(:,k,tt,hh,ss)... %  INTEGRATOR recurrent
        %                             + att_gain_stim * w_int_vis * spikes_vis(:,k,tt,hh,ss)...     %  Visual input
        %                             + att_gain_stim * w_int_vest * spikes_vest(:,k,tt,hh,ss)...     % Vestibular input
        %                         ... % + att_gain_targ * w_int_targ * spikes_target(:,k,tt,hh,ss)...  % No longer here. HH20170317
        %                         );
        
        % Just let the INTEGRATOR to be ideal. (straight sum)
        rate_int_this(:,k+1) = rate_int_this(:,k)...   %  Self dynamics.  in Hz!
            + att_gain_stim * w_int_vis * spikes_vis_this(:,k)...     %  Visual input
            + att_gain_stim * w_int_vest * spikes_vest_this(:,k);     % Vestibular input
        
        % -- Update LIP layer --
        rate_lip_this(:,k+1) = bias_lip + (1-dt/time_const_lip)*rate_lip_this(:,k)...   %  Self dynamics.  in Hz!
            + 1/time_const_lip * (...
            w_lip_lip * spikes_lip_this(:,k)...  %  LIP recurrent
            + att_gain_targ * w_lip_targ * spikes_target_this(:,k)...  % Target input moved here. HH20170317
            + w_lip_int * spikes_int_this(:,k)... %  INTEGRATOR->LIP
            );
        
        % -- Turn rate to binary spike for the next step --
        spikes_int_this(:,k+1)  = (rand(N_int,1) < dt*(rate_int_this(:,k+1)-threshold_int));
        spikes_lip_this(:,k+1)  = (rand(N_lip,1) < dt*(rate_lip_this(:,k+1)-threshold_lip));
        
        % -- Variable used for stopping the integration --
        
        %                     if mm == 1 % Only apply to train set
        decision_ac(:,k+1) = rate_lip_this(:,k+1);
        %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
        %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out(:,k));
        
        % -- Termination --
        if if_bounded && k*dt > stim_on_time + 0.5 && att_gain_stim == 1 ...
                ... && (max(decision_ac(:,k+1)) > decis_thres_this)   % Only Max
                && max(mean(mean(decision_ac(left_targ_ind-5:left_targ_ind+5,max(1,k-20):k+1))),...    % Smoothed Max
                mean(mean(decision_ac(right_targ_ind-5:right_targ_ind+5,max(1,k-20):k+1)))) > decis_thres_this
                ...&& abs(mean(mean(decision_ac(left_targ_ind-5:left_targ_ind+5,max(1,k-20):k+1))) - ...    % Smoothed diff
                ...mean(mean(decision_ac(right_targ_ind-5:right_targ_ind+5,max(1,k-20):k+1)))) > decis_thres_this
            
            
            % Set the attention for motion stimuli to zero
            att_gain_stim = att_gain_this;
            RT(tt) = (k-stim_on_time_in_bins)*dt;
            % last_proba(:,count) = rate_int(:,k,tt,hh,ss);
        end
        %                     end
        
        k=k+1;
    end  % of network dynamics
    
    % == Collect data at the end of this parfor ==
    spikes_target(:,:,tt) = spikes_target_this;
    spikes_vis(:,:,tt)= spikes_vis_this;
    spikes_vest(:,:,tt)= spikes_vest_this;
    spikes_int(:,:,tt)  = spikes_int_this;
    spikes_lip(:,:,tt)  = spikes_lip_this;
    rate_int(:,:,tt) = rate_int_this;
    rate_lip(:,:,tt) = rate_lip_this;
    
    parfor_progress;
    
end % of parfor trial
parfor_progress(0);
fprintf('\n');
par_time = toc(par_time)

% Reorganize data saved in parloop
spikes_target = reshape(spikes_target,N_lip,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
spikes_vis= reshape(spikes_vis,N_vest,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
spikes_vest= reshape(spikes_vest,N_vest,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
spikes_int  = reshape(spikes_int,N_int,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
spikes_lip  = reshape(spikes_lip,N_lip,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
rate_int = reshape(rate_int,N_int,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
rate_lip = reshape(rate_lip,N_lip,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
RT = reshape(RT,N_rep,length(unique_heading),length(unique_stim_type));

% % Pack and save data
% disp('Saving data...');
% if if_debug
%     save('./result/last_result.mat','paras','spikes_target','spikes_vis','spikes_vest','spikes_int','spikes_lip','rate_int','rate_lip','RT','-v7.3');
% else
%     save('./result/last_result.mat','paras','rate_lip','rate_int','RT','-v7.3');
% end

if if_debug
    assignin('base','rate_lip', rate_lip);
    assignin('base','paras', paras);
    assignin('base','rate_int', rate_int);
    assignin('base','spikes_target', spikes_target);
    assignin('base','spikes_vis', spikes_vis);
    assignin('base','spikes_vest', spikes_vest);
    assignin('base','spikes_int', spikes_int);
    assignin('base','spikes_lip', spikes_lip);
    assignin('base','proba_target_tuning', proba_target_tuning);
    assignin('base','RT', RT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Moved to AnalysisResult.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plot results...');
AnalysisResult;


