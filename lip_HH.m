%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for multisensory decision making
% Modified by HH 2017 @ UNIGE
% Adapted for the vestibular-visual multisensory heading discrimintation task
%
% lip_HH(para_override)
% para_override: {'name1',value1; 'name2',value2, ...}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function userResult = lip_HH(para_override,output_result)
%clear
%clear global

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end
%addpath(genpath(pwd));

% rand('state',sum(100*clock));
rng('shuffle')

hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name

if ~isempty(strfind(hostname,'node')) || ~isempty(strfind(hostname,'clc')) % contains(hostname,{'node','clc'})  % ION cluster;
    if isempty(gcp('nocreate'))
        % parpool(hostname(1:strfind(hostname,'.')-1),20);
        parpool('local',20);
    end
    ION_cluster = 1;
    addpath(genpath(pwd));

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

pctRunOnAll warning off stats:obsolete:ReplaceThisWithMethodOfObjectReturnedBy;

if ~exist('./result/','dir')
    mkdir('./result/');
end

% Override
% ION_cluster = 1;

% === Switches ===
if_debug = ~ION_cluster;

h_grouped = [];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============ Sizes ============
% Input layers
N_targets = 2; % Target input

% ============ Times ============
if ION_cluster
    N_vis = 200; % Visual motion signal
    N_vest = 200; % Vestibular motion signal
    
    % Today we decide to add a perfect integrator layer here. HH20170317 in UNIGE
    % This layer has a very long time constant and feeds it's input to LIP, whereas LIP,
    % which has a short time constant, keeps receiving target input but does not integrate it.
    % For most part of the code, I just rename the old "_lip" stuff to "_int"
    N_int = 300;
    N_lip = 300;  % Output layers (Real LIP layer)
    
    % dt = 1e-3; % Size of time bin in seconds
    dt = 5e-3; % Size of time bin in seconds
else
    N_vis = 100; % Visual motion signal
    N_vest = 100; % Vestibular motion signal
    
    N_int = 100;
    N_lip = 100;  % Output layers (Real LIP layer)
    dt = 2e-3;

%{
    N_vis = 200; % Visual motion signal
    N_vest = 200; % Vestibular motion signal
    
    % Today we decide to add a perfect integrator layer here. HH20170317 in UNIGE
    % This layer has a very long time constant and feeds it's input to LIP, whereas LIP,
    % which has a short time constant, keeps receiving target input but does not integrate it.
    % For most part of the code, I just rename the old "_lip" stuff to "_int"
    N_int = 300;
    N_lip = 300;  % Output layers (Real LIP layer)
    
    % dt = 1e-3; % Size of time bin in seconds
    dt = 5e-3; % Size of time bin in seconds

%}
end
trial_dur_total = 1.7; % in s
stim_on_time = 0.2; % in s
motion_duration = 1.5; % in s
% Add 150 ms delay for both and another 100 ms delay for visual to better match the data
delay_for_all = 0.15; % in s
delay_another_for_visual = 0.1; % in s

% ============ Decision bound ============
if_bound_RT = 1; % if_bounded = 1 means that the trial stops at the bound (reaction time version); 
                 % otherwise, use fixed_RT below
read_out_at_the_RT = 0; % Readout decision at RT instead of at the end of each trial
% d = @(x) max(x);  % A bug: if lip_HH is a function, this anonymous function cannot be broadcasted into parfor loop
%  f_bound = @(x) abs(x(right_targ_ind)-x(left_targ_ind));

% Smoothed max
decis_bound = 37*[1 1 1] + 5*[0 0 1]; % bound height, for different conditions
% decis_thres = [inf inf inf]; % bound height, for different conditions

%  decis_bound = 40*[1 1 1+8/29]; % bound height, for different conditions

% Smoothed diff
% decis_bound = 13*[1 1 1+2/13]; % bound height, for different conditions

% ============ if if_bound_RT = 0, use the fixed_RTs ======
fixed_RT = [1.5 1.5 1.5]; % Seconds after stim on

att_gain_stim_after_hit_bound = [0 0 0];

% =============== Conditions ===============
if ION_cluster
    unique_heading = [-8 -4 -2 -1 0 1 2 4 8];
    conflict_heading = 0; % Vis - Vest
    unique_stim_type = [1 2 3];
    N_rep = 200; % For each condition
else
    unique_heading = [-8 -4 -2 -1 0 1 2 4 8];
%     unique_heading = [-8 0 8];
    conflict_heading = 0; % Vis - Vest
    unique_stim_type = [1 2 3];
    N_rep = 50;
%     N_rep = 100; % For each condition

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

A_eta = 1e-5; % Described in Beck 2008 Method
% A_eta = 0.15; % The one I predicted by my estimation that will lead to "average corr less than 90 = 0.2". For dt = 1 ms !! HH20180503

var_vis = A_eta; % Although the Methods of Beck 2008 did not mention the "sqrtm", I assume there was.


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
gain_vel_vis = 7; %   Corresponding to max vel of 0.93 * 7 = 2.45
gain_acc_vest = 2.6; %  gain_vel_vis * sum(vel)/sum(abs(acc)); % (m^2/s)^-1  % Corresponding to max acc of 2.42

% =================== Network configuration ===================
% -- Time constant for integration
time_const_int = 10000e-3; % in s
time_const_lip = 100e-3; % in s

% ---- Visual to INTEGRATOR ----
g_w_int_vest = 20; % Let's vary the gain separately`
g_w_int_vis = 20; 
dc_w_int_vis = 0;
k_int_vis = 4; % Larger, narrower
k_int_vis_along_vis = 0.1; % Larger, wider

gamma = 0; % Time_varying_weight_factor = reliability ^ gamma; when gamma = 0, weights are constant over time

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
g_w_lip_int = 15;
k_lip_int = 10;
dc_w_lip_int = -3.6;

amp_I_lip_int = 0;  % Mexico hat shape
k_lip_int_I = 2;

% ----- LIP recurrent connection ------
g_w_lip_lip = 5;
k_lip_lip = 10;
dc_w_lip_lip = -3;

amp_I_lip = 0;  % Mexico hat shape
k_lip_I = 2;
bias_lip = 0;

% Input-output function of LIP
threshold_lip = 0.0;

% Others
save_folder = '';

% ==== Heterogeneity ====
heter_enable = 1; % Master switch of heterogeneity

% --- "POST": Add heterogeneity in post synaptic parameters ---
% Variability in:    g       k     dc                   
heter_post =  heter_enable * [0.2     0.5     0 ;  % vest -> int
                     0.2     0.5     0 ;  % vis -> int
                     0       0.2   0.1 ;  % int -> lip
                     0       0.2   0.2 ]; % lip -> lip
                 
heter_post (:) = 0;

% --- "Normal": Add normal distributed noise in the weight matrix ---
heter_normal = heter_enable * 0 * [1 1 1 1];  % vest -> int, vis -> int, int -> lip, lip -> lip

% --- "Dropout": Increase sparseness in the weight matrix ---
heter_dropout = heter_enable * 00 * [1 1 1 1]; % vest -> int, vis -> int, int -> lip, lip -> lip

% --- "LogNormal": log normal distribution for each group of weights (diagonal) ---
heter_lognormal  = heter_enable * 1 * [1 1 1 1]; % vest -> int, vis -> int, int -> lip, lip -> lip

vis_vest_weight_noise_cor = heter_enable * -0.8 * 0;  % Controls the correlation between noise of weight in vest -> int and vis -> int

%%%%%%%%%%%%%%%% Override by input argument %%%%%%%%%%%%%%%
para_override_txt = '';                 
weights_override = 0;
 
if nargin >= 1
    len = size(para_override,1);
    for ll = 1:len
        if exist(para_override{ll,1},'var')
            eval([para_override{ll,1} '= para_override{ll,2};']);
            
            fprintf('Overriding %s = %s...\n',para_override{ll,1},num2str(para_override{ll,2}));
            if strcmp(para_override{ll,1},'save_folder')
                if ~exist(['./result/' save_folder],'dir')
                    mkdir(['./result/' save_folder]);
                end
            elseif strcmp(para_override{ll,1},'weights_override')
                weights_override = 1;
                load(para_override{ll,2});
            else
                para_override_txt = [para_override_txt sprintf('_%s_%s',para_override{ll,1},num2str(para_override{ll,2}))];
            end
        else
            fprintf('Parameter ''%s'' not found...\n',para_override{ll,1});
        end
    end
end

% [~,gitlog] = system('git describe --always');
% gitlog = sprintf('%s',gitlog(1:7));
cc = clock;
cc = sprintf('%g%02g%02g%g%g',cc(1:5));
% save_folder = [save_folder cc '_' gitlog '_'];
save_folder = [save_folder cc '_' ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Pure Parameters End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some preparation (after the pure parameter session in case some paras are overriden)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Some weights ---
% Vestibular
equivalent_conherence = coherence;
r_spont_vest = r_spont_vis;
b_pref_vest = b_pref_vis;
b_null_vest = b_null_vis;
K_vest = K_vis;
K_cov_vest = K_cov_vis;
var_vest = var_vis;

% ---- Vestibular to INTEGRATOR ----
dc_w_int_vest = dc_w_int_vis;
k_int_vest = k_int_vis;
k_int_vest_along_vest = k_int_vis_along_vis;


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
miu = motion_duration/2 + delay_for_all; 
sigma = motion_duration/2/num_of_sigma;

use_real_profile = 1;

if use_real_profile
    % -- Use the real velocity profile instead. HH20170808 --
    load('Gaussian_vel_real_sigma3p5.mat');
    real_t = Gaussian_vel_real_sigma3p5(:,1)/1000;
    real_v = Gaussian_vel_real_sigma3p5(:,2);
    real_v_interp = interp1(real_t,real_v,t_motion);
    vel = real_v_interp;
    
    vel = vel*amp/sum(vel*dt) ; % in m/s. Normalized to distance
    acc = diff(vel)/dt; % in m/s^2
    
    % Shift vel with visual delay
    shift_bins = round(delay_another_for_visual/dt);
    vel_vis = [zeros(1,shift_bins) vel(1:end-shift_bins)];
    
    
else  % -- Use the ideal velocity profile ---
    
    vel = exp(-(t_motion-miu).^2/(2*sigma^2));
    vel = vel*amp/sum(vel*dt) ; % in m/s. Normalized to distance
    acc = diff(vel)/dt; % in m/s^2
    
    % Shift vel with visual delay
    vel_vis = exp(-(t_motion-miu-delay_another_for_visual).^2/(2*sigma^2)); % -- Use the ideal velocity profile ---
    vel_vis = vel_vis*amp/sum(vel_vis*dt) ; % in m/s. Normalized to distance
end


% To make sure t_motion, vel, and acc have the same length
t_motion(end) = [];
vel_vis(end) = [];
vel(end) = [];

% Time-varying synaptic weights (sensory -> INT) in proportion to reliability
%%
% gamma = 5; % I put it earlier
w_int_vest_timevarying_factor = [zeros(1,stim_on_time_in_bins) (abs(acc)).^gamma/max((abs(acc)).^gamma)];
w_int_vis_timevarying_factor = [zeros(1,stim_on_time_in_bins) vel_vis.^gamma/max(vel_vis.^gamma)];

% Plot time-varying factor
if if_debug
    ts = ((0:trial_dur_total_in_bins-1)-stim_on_time_in_bins)*dt; % in s
    figure(1439);  clf; hold on;
    plot(ts,w_int_vest_timevarying_factor,'linew',2);
    plot(ts,w_int_vis_timevarying_factor,'linew',2);
    plot(xlim,[1 1],'k--');
    SetFigure(15); ylabel('Time-varying factor'); title(['gamma = ' num2str(gamma)]);
    xlim([0 1.5]); ylim([0 1.1])
end

%%
if if_debug
    figure(111); clf
    set(gcf,'name','Motion profile','uni','norm','pos',[0.632       0.381       0.358       0.403]);
    h = plotyy(t_motion,acc,t_motion,vel_vis);
    ylabel(h(2),'Velocity (m/s)');
    ylabel(h(1),'Acceleration (m^2/s)');
   	% keyboard;
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
if weights_override  % Use the specific weights. HH20170705
    w_int_vest = weights_saved.w_int_vest;
    w_int_vis = weights_saved.w_int_vis;
    w_int_int = weights_saved.w_int_int;
    w_lip_int = weights_saved.w_lip_int;
    w_lip_lip = weights_saved.w_lip_lip;
    w_lip_targ = weights_saved.w_lip_targ;
    disp('Weights overridden...');
else
    
    w_int_vis = zeros(N_int,N_vis);
    w_int_vest = zeros(N_int,N_vest);
    % w_int_targ = zeros(N_int,N_int); % Not N_target, but N_int (supposing the same number as the lip neurons)
    w_lip_targ = zeros(N_lip,N_lip); % Now integration layer does not receive target input, but LIP layer does. HH20170317
    w_int_int = zeros(N_int,N_int);
    w_lip_int = zeros(N_lip,N_int);
    w_lip_lip = zeros(N_lip,N_lip);
    
    % ==== Heterogeneity Method 1: Add heterogeneity in post synaptic parameters ("POST") ====
    
    g_w_int_vest_heter = (1+randn(1,N_int)*heter_post(1,1)) * g_w_int_vest; % (mean = std)
    k_int_vest_heter = (1+randn(1,N_int)*heter_post(1,2)) * k_int_vest; k_int_vest_heter(k_int_vest_heter<=1) = 1;
    dc_w_int_vest_heter = (1+randn(1,N_int)*heter_post(1,3)) * dc_w_int_vest - randn(1,N_int)*0.1; dc_w_int_vest_heter(dc_w_int_vest_heter>0) = 0;
    
    g_w_int_vis_heter = (1+randn(1,N_int)*heter_post(2,1)) * g_w_int_vis; % (mean = std)
    k_int_vis_heter = (1+randn(1,N_int)*heter_post(2,2)) * k_int_vis;  k_int_vis_heter(k_int_vis_heter<=1) = 1;
    dc_w_int_vis_heter = (1+randn(1,N_int)*heter_post(2,3)) * dc_w_int_vis - randn(1,N_int)*0.1;  dc_w_int_vis_heter(dc_w_int_vis_heter>0) = 0;
    
    g_w_lip_int_heter = (1+randn(1,N_lip)*heter_post(3,1)) * g_w_lip_int; % (mean = std)
    k_lip_int_heter = (1+randn(1,N_lip)*heter_post(3,2)) * k_lip_int; k_lip_int_heter(k_lip_int_heter<=1) = 1;
    dc_w_lip_int_heter = (1+randn(1,N_lip)*heter_post(3,3)) * dc_w_lip_int;
    
    g_w_lip_lip_heter = (1+randn(1,N_lip)*heter_post(4,1)) * g_w_lip_lip; % (mean = std)
    k_lip_lip_heter = (1+randn(1,N_lip)*heter_post(4,2)) * k_lip_lip; k_lip_lip_heter(k_lip_lip_heter<=1) = 1;
    dc_w_lip_lip_heter = (1+randn(1,N_lip)*heter_post(4,3)) * dc_w_lip_lip;
    
    for nn=1:N_int
        
        % -- VIS->Int, Vest->Int --
        %     w_int_vis(nn,:) = g_w_int_vis/N_vis *(exp(k_int_vis * (cos((prefs_int(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
        %         .* abs(sin(prefs_vis/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        %         + dc_w_int_vis/N_vis;   % Offset
        %
        %     w_int_vest(nn,:) = g_w_int_vest/N_vest *(exp(k_int_vest * (cos((prefs_int(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
        %         .* abs(sin(prefs_vest/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        %         + dc_w_int_vest/N_vest;   % Offset
        
        w_int_vis(nn,:) = g_w_int_vis_heter(nn)/N_vis *(exp(k_int_vis_heter(nn) * (cos((prefs_int(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
            .* abs(sin(prefs_vis/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
            + dc_w_int_vis_heter(nn)/N_vis;   % Offset
        
        w_int_vest(nn,:) = g_w_int_vest_heter(nn)/N_vest *(exp(k_int_vest_heter(nn) * (cos((prefs_int(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
            .* abs(sin(prefs_vest/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
            + dc_w_int_vest_heter(nn)/N_vest;   % Offset
        
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
    
    
    for nn=1:N_lip % For each POST-synaptic cell
        
        % -- Targ->LIP --
        w_lip_targ(nn,:) = g_w_lip_targ/N_int *(exp(k_lip_targ*(cos((prefs_lip-prefs_lip(nn))/360 *2*pi)-1)));  %  Target input
        
        % -- Int->LIP --
        %     w_lip_int(nn,:) = g_w_lip_int/N_int*...
        %         ((exp(k_lip_int*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1)))-...
        %         amp_I_lip_int*(exp(k_lip_int_I*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1))))...
        %         + dc_w_lip_int/N_int;
        w_lip_int(nn,:) = g_w_lip_int_heter(nn)/N_int*...
            ((exp(k_lip_int_heter(nn)*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1)))-...
            amp_I_lip_int*(exp(k_lip_int_I*(cos((prefs_int-prefs_lip(nn))/360*2*pi)-1))))...
            + dc_w_lip_int_heter(nn)/N_int;
        
        % -- LIP->LIP --
        %     w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...
        %         ((exp(k_lip_lip*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
        %         amp_I_lip*(exp(k_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
        %         + dc_w_lip_lip/N_lip;
        w_lip_lip(nn,:) = g_w_lip_lip_heter(nn)/N_lip*...
            ((exp(k_lip_lip_heter(nn)*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
            amp_I_lip*(exp(k_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
            + dc_w_lip_lip_heter(nn)/N_lip;
    end
    
    if heter_enable
        
        ws = {'int_vest','int_vis','lip_int','lip_lip'};
        
        w_old_int_vest = w_int_vest;
        w_old_int_vis = w_int_vis;
        
        for ww = 1:length(ws)
            
            % ==== Heterogeneity Method 3 : Add Gaussian noise in weights ("Normal") ====
            % w_int_vest = w_int_vest + randn(size(w_int_vest))*std(w_int_vest(:))*heter_normal(1);
            eval(sprintf('w_%s =  w_%s + randn(size(w_%s))*std(w_%s(:))*heter_normal(%g);',ws{ww},ws{ww},ws{ww},ws{ww},ww));
            
            % ==== Heterogeneity Method 4: Log ormal (proportional to the mean value) ====
            % w_int_vest = sign(w_int_vest).* abs(w_int_vest).^(1+randn(size(w_int_vest)).* abs(w_int_vest) * heter_lognormal(1));
            % eval(sprintf('w_%s = sign(w_%s).* abs(w_%s).^(1+randn(size(w_%s)).* abs(w_%s) * heter_lognormal(%g));',...
            %     ws{ww},ws{ww},ws{ww},ws{ww},ws{ww},ww));
            
            % Calculate designed mean and variance
            eval(sprintf(' m = abs(w_%s); v = (abs(w_%s) * heter_lognormal(%g)).^2;',ws{ww},ws{ww},ww));
            
            % Calculate parameters for log normal
            miu = log(m.^2./sqrt(v+m.^2));
            sigma = sqrt(log(v./m.^2+1));
            
            % Generate lognormal and recover the sign
            eval(sprintf('w_%s = sign(w_%s) .* exp( miu + randn(size(w_%s)) .* sigma );',ws{ww},ws{ww},ws{ww}));
            
            % ==== Heterogeneity Method 5: Dropout while keep the mean unchanged ====
            % w_int_vest(rand(size(w_int_vest))<(1)) = 0;
            eval(sprintf('w_%s(rand(size(w_%s)) < heter_dropout(%g)) = 0;',ws{ww},ws{ww},ww));
            % Scale to keep the mean unchanged
            eval(sprintf('w_%s = w_%s / (1-heter_dropout(%g));',ws{ww},ws{ww},ww));
            
        end
        
        % %{
        % Just to verify the distribution of diagonals of w_lip_lip
        % userResult.w_lip_lip = w_lip_lip;
        
        figure(1318); clf; hold on;
        delta_theta = 0;
        
        off_diag = round(delta_theta/(360/size(w_lip_lip,1)));
        x = diag(w_lip_lip,off_diag);
        
        if range(x)>0
            [N,edges] = hist(x,min(x):0.0003:max(x));
            plot([edges edges(end)+edges(2)-edges(1)],[smooth(N,1);0]/sum(N)/(edges(2)-edges(1)),'k-','linew',2);
        end
        
        ylabel('Prob density');
        xlabel(['LIP --> LIP, \Delta\theta = ' num2str(delta_theta)]);
        
        max_y = max(ylim);
        plot(mean(x)*ones(1,2),[0 max_y*1.05],'k--')
        
        %      h_grouped = [h_grouped gca];
        
        %  -- Noise correlation in int_vest and int_vis (quick and dirty) 20170613--
        %%
        figure(128); clf
        vest_noise = w_int_vest - w_old_int_vest;
        vis_noise = w_int_vis - w_old_int_vis;
        plot(vest_noise,vis_noise,'k.');
        
        vest_correlated = (1-vis_vest_weight_noise_cor) * vest_noise + vis_vest_weight_noise_cor * vis_noise;
        vis_correlated =  vis_vest_weight_noise_cor * vest_noise + (1-vis_vest_weight_noise_cor) * vis_noise;
        
        hold on; plot(vest_correlated,vis_correlated,'r.');
        [r,p] = corr(vest_correlated(:),vis_correlated(:));
        title(sprintf('r = %g, p = %g',r,p));
        %%
        w_int_vest = w_old_int_vest + vest_correlated;
        w_int_vis = w_old_int_vis + vis_correlated;
        
        %% Make the average of the weight matrix symmetric for heter
        if heter_enable
            temp = mat2cell(w_int_vest,size(w_int_vest,1)/2*[1 1],size(w_int_vest,2)/2*[1 1]);
            aver_temp = cellfun(@mean,cellfun(@mean,temp,'UniformOutput',false));
            temp{1} = temp{1}/aver_temp(1) * mean(aver_temp([1 4]));
            temp{4} = temp{4}/aver_temp(4) * mean(aver_temp([1 4]));
            temp{2} = temp{2}/aver_temp(2) * mean(aver_temp([2 3]));
            temp{3} = temp{3}/aver_temp(3) * mean(aver_temp([2 3]));
            w_int_vest = cell2mat(temp);
            
            temp = mat2cell(w_int_vis,size(w_int_vis,1)/2*[1 1],size(w_int_vis,2)/2*[1 1]);
            aver_temp = cellfun(@mean,cellfun(@mean,temp,'UniformOutput',false));
            temp{1} = temp{1}/aver_temp(1) * mean(aver_temp([1 4]));
            temp{4} = temp{4}/aver_temp(4) * mean(aver_temp([1 4]));
            temp{2} = temp{2}/aver_temp(2) * mean(aver_temp([2 3]));
            temp{3} = temp{3}/aver_temp(3) * mean(aver_temp([2 3]));
            w_int_vis = cell2mat(temp);
        end
        
    end % If heter
    
    % Return immediately for debugging
    % if nargin == 2 % Save result
    %     for rr = 1:length(output_result)
    %         eval(['result.(output_result{rr}) =' output_result{rr}]);
    %     end
    % else
    %     result =[];
    % end
    %
    % return;
    %}

end

if if_debug
    figure(9111); clf;
    set(gcf,'uni','norm','pos',[0.0053571     0.25429      0.9869     0.27048],'name','Weights & Raster plot');
    
    subplot(1,4,4);
%     axes('position',[0.08 0.725 0.23 0.241]);
    h1 = imagesc(prefs_lip,prefs_lip,w_lip_lip'); colorbar; axis  xy; title('LIP->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{lip}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    set(h1,'ButtonDownFcn',{@plot_weight,prefs_lip,prefs_lip,w_lip_lip'});
    axis square
    
    subplot(1,4,3);
%     axes('position',[0.08 0.384 0.23 0.241]);
    h2 = imagesc(prefs_int,prefs_lip,w_lip_int'); colorbar; axis xy; title('Int->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    set(h2,'ButtonDownFcn',{@plot_weight,prefs_lip,prefs_int,w_lip_int'});
    axis square
    %     subplot(4,3,7);
    %     % surf(prefs_int,prefs_int,w_int_int');
    %     imagesc(prefs_int,prefs_int,w_int_int');    colorbar; axis xy; title('Int->Int');
    %     xlabel('\theta_{int}'); ylabel('\theta_{int}');
    %     set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    %     hold on; plot(prefs_int,w_int_int(:,end/2)/range(w_int_int(:,end/2))*100,'linew',3,'color','c');
    %     plot(xlim,[0 0],'--c');
    
    %     subplot(4,3,7);
    %     imagesc(prefs_int,prefs_vest,w_int_vest'); colorbar; axis  xy; title('Vest->Int');
    %     xlabel('\theta_{int}'); ylabel('\theta_{vest}'); ylim([-20 20]);
    %     set(gca,'xtick',-180:90:180);
    
    subplot(1,4,2);
%     axes('position',[0.08 0.055 0.23 0.241]);
    h3 = imagesc(prefs_vis,prefs_int,w_int_vis'); colorbar; axis xy; title('Vis->Int');
    xlabel('\theta_{int}'); ylim([-180 180]);
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    set(h3,'ButtonDownFcn',{@plot_weight,prefs_int,prefs_vis,w_int_vis'});
    ylabel('\theta_{vis}');
    axis square
%     temp_ax = axes('Pos',[0.05 0.124 0.053 0.134]);
%     plot(prefs_vis,gain_func_along_vis(prefs_vis)); hold on;
%     plot([-30 -30],ylim,'r-'); plot([30 30],ylim,'r-');
%     view(270,90);     set(gca,'xtick',-180:90:180);
%     xlabel('\theta_{vis/vest}');

    subplot(1,4,1);
%     axes('position',[0.08 0.055 0.23 0.241]);
    h4 = imagesc(prefs_vest,prefs_int,w_int_vest'); colorbar; axis xy; title('Vest->Int');
    xlabel('\theta_{int}'); ylim([-180 180]);
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    set(h4,'ButtonDownFcn',{@plot_weight,prefs_int,prefs_vest,w_int_vest'});
    ylabel('\theta_{vest}');
    axis square
  
    colormap hot;
    % keyboard;
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

rate_vest = zeros(N_vest,trial_dur_total_in_bins,n_parfor_loops);
rate_vis = zeros(N_vis,trial_dur_total_in_bins,n_parfor_loops);
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
    rate_vis_for_each_heading(:,hh) = ((max_rate_vis-b_vis)*exp(K_vis*(cos((prefs_vis'-(unique_heading(hh)+conflict_heading/2))/360*2*pi)-1))+b_vis);
    
    %     max_rate_vest = r_spont_vest + b_pref_vis * coherence;
    %     b_vis = r_spont_vis + b_null_vis * coherence;
    
    % Vestibular response does not depend on coherence.
    % Here I just set the vestibular activity similar to visual response under 'equivalent_coh' coh
    max_rate_vest = r_spont_vest + b_pref_vest * equivalent_conherence;
    b_vest = r_spont_vest + b_null_vest * equivalent_conherence;
    rate_vest_for_each_heading(:,hh) = ((max_rate_vest-b_vest)*exp(K_vest*(cos((prefs_vest'-(unique_heading(hh)-conflict_heading/2))/360*2*pi)-1))+b_vest);
    
    %         if if_debug
    %             figure(91);clf;
    %             plot(prefs_vis,proba_vis/dt,'r');
    %             plot(prefs_vest,proba_vest/dt,'b');
    %             title('Sensory population tuning');
    %         end
end

for ss = 1:length(unique_stim_type)
    % === Some slicing stuffs necessary for parfor ===
    att_gain_this(ss) = att_gain_stim_after_hit_bound(unique_stim_type(ss));
    decis_bound_this(ss) = decis_bound(unique_stim_type(ss));
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

    % I realize that I should set random seed for each parallel workers! (set at the beginning of the file does not count)
    % Otherwise, the workers are using actually the same series of random numbers across trials, which will lead to "bumps" in PSTH
    % HH20180622
    % rand('state',sum(100*clock));
    
    % -- Get data from cheat sheet --
    cond_this = ss_hh_cheatsheet(tt,:);
    ss_this = cond_this(1);
    hh_this = cond_this(2);
    
    % -- Target input spike train --
    spikes_target_this = rand(N_lip,trial_dur_total_in_bins)<(aux_proba_target);
    
    % -- Visual input spike train --
    rate_vis_this = rate_vis_for_each_heading(:,hh_this)*[zeros(1,stim_on_time_in_bins) vel_vis*gain_vel_vis]...
        + w_cov_vis / dt * randn(N_vis,trial_dur_total_in_bins) ... % Note that the noise term depends on dt!! HH20171128
        .* 1 ;% repmat([zeros(1,stim_on_time_in_bins) ones(size(vel_vis))],N_vis,1);  % No reason to turn off noise before stimuli. HH20171128
    
    % -- Vestibular ---
    % With the temporal gain of abs(acc)
    rate_vest_this = rate_vest_for_each_heading(:,hh_this)*[zeros(1,stim_on_time_in_bins) abs(acc)*gain_acc_vest]...
        + w_cov_vest / dt * randn(N_vest,trial_dur_total_in_bins) ...   % Note that the noise term depends on dt!!
        .* 1; % repmat([zeros(1,stim_on_time_in_bins) ones(size(vel_vis))],N_vis,1);
    
    % -- Stimulus condition selection --
    % Or alternatively, I should only turn off the sensory input but let the correlated noise unchanged?
    if unique_stim_type(ss_this) == 1
        rate_vis_this = 0 * rate_vis_this; % Shut down visual activity
    elseif unique_stim_type(ss_this) ==2
        rate_vest_this = 0 * rate_vest_this; % Shut down vestibular activity
    end
    
    spikes_vis_this = rand(N_vis,trial_dur_total_in_bins)<(rate_vis_this) * dt; % Put dt here. HH20171128
    spikes_vest_this = rand(N_vest,trial_dur_total_in_bins)<(rate_vest_this) * dt;
    
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
    decision_ac = zeros(N_lip,trial_dur_total_in_bins+2);
    
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
%         rate_int_this(:,k+1) = rate_int_this(:,k)...   %  Self dynamics.  in Hz!
%             + att_gain_stim * w_int_vis * spikes_vis_this(:,k)...     %  Visual input 
%             + att_gain_stim * w_int_vest * spikes_vest_this(:,k);     % Vestibular input

        rate_int_this(:,k+1) = rate_int_this(:,k)...   %  Self dynamics.  in Hz!
            + att_gain_stim * w_int_vis * w_int_vis_timevarying_factor(k) * spikes_vis_this(:,k)...     % Visual input with reliablity-dependent weights HH20171026
            + att_gain_stim * w_int_vest * w_int_vest_timevarying_factor(k) * spikes_vest_this(:,k);     % Vestibular input with reliablity-dependent weights
        
        
        % -- Update LIP layer --
        rate_lip_this(:,k+1) = bias_lip + (1-dt/time_const_lip)*rate_lip_this(:,k)...   %  Self dynamics.  in Hz!
            + 1/time_const_lip * (...
            w_lip_lip * spikes_lip_this(:,k)...  %  LIP recurrent
            + att_gain_targ * w_lip_targ * spikes_target_this(:,k)...  % Target input moved here. HH20170317
            + w_lip_int * spikes_int_this(:,k)... %  INTEGRATOR->LIP
            );
        
%         %%%%%%%%%%%%%%%%%% Try normalization across modality %%%%%%%%%%%%%%%%%%%%%
%         if unique_stim_type(ss_this) == 3
%             rate_lip_this(:,k+1) = rate_lip_this(:,k+1)/1.1;
%         end
%         %%%%%%%%%%%%%%%%%% Try normalization across modality %%%%%%%%%%%%%%%%%%%%%
        
        % -- Turn rate to binary spike for the next step --
        spikes_int_this(:,k+1)  = (rand(N_int,1) < dt*(rate_int_this(:,k+1)-threshold_int));
        spikes_lip_this(:,k+1)  = (rand(N_lip,1) < dt*(rate_lip_this(:,k+1)-threshold_lip));
        
        % -- Variable used for stopping the integration --
        
        %                     if mm == 1 % Only apply to train set
        decision_ac(:,k+1) = rate_lip_this(:,k+1);
        %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
        %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out(:,k));
        
        % -- Termination --
        
        if if_bound_RT        % - Rate termination --
            if att_gain_stim == 1 && k*dt > stim_on_time + 0.5  ...
                    ... && (max(decision_ac(:,k+1)) > decis_bound_this)   % Only Max
                    && max(mean(mean(decision_ac(left_targ_ind-5:left_targ_ind+5,max(1,k-20):k+1))),...    % Smoothed Max
                    mean(mean(decision_ac(right_targ_ind-5:right_targ_ind+5,max(1,k-20):k+1)))) > decis_bound_this(ss_this)
                ...&& abs(mean(mean(decision_ac(left_targ_ind-5:left_targ_ind+5,max(1,k-20):k+1))) - ...    % Smoothed diff
                    ...mean(mean(decision_ac(right_targ_ind-5:right_targ_ind+5,max(1,k-20):k+1)))) > decis_bound_this(ss_this)
                    %}
                % Set the attention for motion stimuli to zero
                att_gain_stim = att_gain_this(ss_this);
                RT(tt) = (k-stim_on_time_in_bins)*dt;
                % last_proba(:,count) = rate_int(:,k,tt,hh,ss);
            end
        else    % - Time termination (Just for para_scanning other parameters!) -
            if att_gain_stim == 1 && ...
                    ((unique_stim_type(ss_this)==1 && k*dt > stim_on_time + fixed_RT(1))...
                    ||(unique_stim_type(ss_this)==2 && k*dt > stim_on_time + fixed_RT(2))...
                    ||(unique_stim_type(ss_this)==3 && k*dt > stim_on_time + fixed_RT(3)))
                % Set the attention for motion stimuli to zero
                att_gain_stim = att_gain_this(ss_this);
                RT(tt) = (k-stim_on_time_in_bins)*dt;
                % last_proba(:,count) = rate_int(:,k,tt,hh,ss);
            end
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
    
    rate_vest(:,:,tt) = rate_vest_this/dt;
    rate_vis(:,:,tt) = rate_vis_this/dt;
    rate_int(:,:,tt) = rate_int_this;
    rate_lip(:,:,tt) = rate_lip_this;
    
    try
        parfor_progress;
    end
    
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
rate_vest = reshape(rate_vest,N_vest,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
rate_vis = reshape(rate_vis,N_vis,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
rate_int = reshape(rate_int,N_int,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
rate_lip = reshape(rate_lip,N_lip,trial_dur_total_in_bins,N_rep,length(unique_heading),length(unique_stim_type));
RT = reshape(RT,N_rep,length(unique_heading),length(unique_stim_type));

% Generate the correct answer: LEFT = -1, RIGHT = 1; random for 0 headings
correct_ans_all = (rand(n_parfor_loops,1) < (sign(unique_heading(ss_hh_cheatsheet(:,2))') + 1)/2) * 2 - 1;
correct_ans_all = reshape(correct_ans_all,N_rep,length(unique_heading),length(unique_stim_type));

% % Pack and save data
% disp('Saving data...');
% if if_debug
%     save('./result/last_result.mat','paras','spikes_target','spikes_vis','spikes_vest','spikes_int','spikes_lip','rate_int','rate_lip','RT','-v7.3');
% else
%     save('./result/last_result.mat','paras','rate_lip','rate_int','RT','-v7.3');
% end

if if_debug || 1
    to_assignin = {'ION_cluster','rate_lip','rate_vest','rate_vis','paras','rate_int','spikes_target','spikes_vis','spikes_vest','spikes_int','spikes_lip',...
                   'proba_target_tuning','RT','correct_ans_all','w_lip_lip','w_lip_int','w_int_vis','w_int_vest'};
    for tta = 1:length(to_assignin)
        assignin('base',to_assignin{tta},eval(to_assignin{tta}));
    end
end

% Save all the weight matrices to replicate the results afterwards
weights_saved.w_int_vest = w_int_vest;
weights_saved.w_int_vis = w_int_vis;
weights_saved.w_int_int = w_int_int;
weights_saved.w_lip_int = w_lip_int;
weights_saved.w_lip_lip = w_lip_lip;
weights_saved.w_lip_targ = w_lip_targ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Moved to AnalysisResult.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plot results...');
AnalyzeResult;

%% Save result
if nargin == 2 % Save result
    for rr = 1:length(output_result)
        eval(['userResult.(output_result{rr}) =' output_result{rr}]);
    end
else
    userResult = [];
end

save(sprintf('./result/%sresults.mat',save_folder),'userResult',...
                'weights_saved','paras','diff_PSTH_correct_mean_headings','ps_psycho');  % Mandatory

if ~ION_cluster
    tmp1 = whos;
    for i = 1:length(tmp1)
        assignin('base',tmp1(i).name,eval(tmp1(i).name));
    end
    % keyboard;
end

end

function plot_weight(~,~,x,y,w)
persistent h;
pos = get(gca,'CurrentPo');
xx = pos(1);
yy = pos(3);
[~,x_ind] = min(abs(xx-x));
[~,y_ind] = min(abs(yy-y));

hold on;  if ishandle(h); delete(h); end
h(1) = plot(x,w(y_ind,:)/max(abs(w(:)))*max(ylim)/2,'linew',3,'color','c');
h(2) = plot(xlim,[0 0],':c');
h(3) = plot(xlim,[yy yy],'--c');

h(4) = plot(w(:,x_ind)/max(abs(w(:)))*max(xlim)/2,y,'linew',3,'color','c');
h(5) = plot([0 0],ylim,':c');
h(6) = plot([xx xx],ylim,'--c');

end


