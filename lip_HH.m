%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
% Modified by HH 2017 @ UNIGE
% Adapted for the vestibular-visual multisensory heading discrimintation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear
clear global

load h

tic;

global Y_targ X_train;
rand('state',sum(100*clock))

% %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === Switches ===
if_debug = 1;
if_test_set = 0;

% === Sizes ===
% Input layers
N_vis = 100; % Visual motion signal
prefs_vis = linspace(-180,180,N_vis);

N_vest = 100; % Vestibular motion signal
prefs_vest = linspace(-180,180,N_vis);

N_targets = 2; % Target input

% Today we decide to add a perfect integrator layer here. HH20170317 in UNIGE
% This layer has a very long time constant and feeds it's input to LIP, whereas LIP,
% which has a short time constant, keeps receiving target input but does not integrate it.
% For most part of the code, I just rename the old "_lip" stuff to "_int"
N_int = 100;
prefs_int = linspace(-180,180,N_int);

% Output layers (Real LIP layer)
N_lip = 100;
prefs_lip = linspace(-180,180,N_lip);

[~,right_targ_ind] = min(abs(prefs_lip - 90)); % Left and right for the heading task
[~,left_targ_ind] = min(abs(prefs_lip - -90));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision bound
if_bounded = 0; % if_bounded = 1 means that the trial stops at the bound (reaction time version)
f_bound = @(x) abs(x(right_targ_ind)-x(left_targ_ind));
decis_thres = [3 6 3]; % bound height, for different conditions
att_gain_stim_after_hit_bound = [0 0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === Times ===
dt = 2e-3; % Size of time bin in seconds
trial_dur_total = 1.7; % in s
stim_on_time = 0.2; % in s
motion_duration = 1.5; % in s

trial_dur_total_in_bins = round(trial_dur_total/dt); %% in time bins (including pre_trial_dur)
stim_on_time_in_bins = stim_on_time/dt; % Time of motion start, in time bins.

ts = ((0:trial_dur_total_in_bins-1)-stim_on_time_in_bins)*dt; % in s

if stim_on_time_in_bins>trial_dur_total_in_bins
    error('pre_trial_dur cannot be longer than trial_dur');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Stimuli ===
% unique_heading = [0 1 2 4 8];
% unique_condition = [1 2 3];
unique_heading = [ 8];
unique_condition = [3];
N_trial = 10; % For each condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


colors = [0 0 1; 1 0 0; 0 0.8 0.4];


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
equivalent_conherence = 12;
r_spont_vest = 10;
b_pref_vest = 1.7;
b_null_vest = -0.2;
K_vest = 1.5;
K_cov_vest = 2;
var_vest = 1e-5;

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
max_rate_target = 42; %% in Hz
b_target = 0;
K_target = 4;
slope_norm_target=0.5;
dc_norm_target = 0.5;

% == Temporal dynamics ==
t_motion = 0:dt:trial_dur_total-stim_on_time; % in s
num_of_sigma = 3.5; amp = 0.2;
miu = motion_duration/2; sigma = motion_duration/2/num_of_sigma;
vel = exp(-(t_motion-miu).^2/(2*sigma^2));
vel = vel*amp/sum(vel*dt) ; % in m/s. Normalized to distance
acc = diff(vel)/dt; % in m/s^2

% To make sure t_motion, vel, and acc have the same length
t_motion(end) = [];
vel(end) = [];

% Gains
gain_vel_vis = 5; % (m/s)^-1
gain_acc_vest = 2; % (m^2/s)^-1

if if_debug
    figure(111);clf
    set(gcf,'name','Motion profile','uni','norm','pos',[0.632       0.381       0.358       0.403]);
    h = plotyy(t_motion,vel,t_motion,acc);
    ylabel(h(1),'Velocity (m/s)');
    ylabel(h(2),'Acceleration (m^2/s)');
end


% === Network configuration ===

% -- Time constant for integration
time_const_int = 10000e-3; % in s
time_const_lip = 100e-3; % in s

% -- Visual to INTEGRATOR
g_w_int_vis = 25;
K_int_vis = 5;
dc_w_int_vis = -0;

% -- Vestibular to INTEGRATOR
g_w_int_vest = 25;
K_int_vest = 5;
dc_w_int_vest = -0;

% --- Targets to LIP ----
g_w_lip_targ= 30;
K_lip_targ= 5;
att_gain_targ = 1; % Drop in attention to visual target once motion stimulus appears.

% --- Recurrent connectivity in INTEGRATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_int_int = 35;
K_int_int = 10;
dc_w_int_int = -11;

amp_I_int = 0;  % Mexico hat shape
K_int_I = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input-output function of INTEGRATOR
bias_int = 0;
threshold_int = 0.0;

% --- INTEGRATOR to the real LIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_lip_int = 25;
K_lip_int = 10;
dc_w_lip_int = 0;

amp_I_lip_int = 0;  % Mexico hat shape
K_lip_int_I = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- LIP recurrent connection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_lip_lip = 25;
K_lip_lip = 8;
dc_w_lip_lip = -11;

amp_I_lip = 0;  % Mexico hat shape
K_lip_I = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input-output function of LIP
bias_lip = 0;
threshold_lip = 0.0;


% === For statistics ===

shuffle_flag = 0; %shuffles trials and computes shuffled Fisher information
info_out_flag = 1; %controls whether Fisher information in the output layer is being computed
info_in_flag = 1; %controls whether Fisher information in the input layer is being computed

% duration of the window used to compute Fisher information
fisher_wind_size = 50;
fisher_wind_offset = 0;
fisher_N_window = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization (I moved here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I reorganized the variable naming and storage conventions. HH20170302
% All data have the size of (neurons, time, trial, stimuli){train/test}
% Three types of data are used:
% 1. firing rate in Hz (rate_xxx),
% 2. firing probability for each bin (proba_xxx),
% 3. binary spike train for each bin (spikes_xxx)
% The network dynamics are computed like: Rate_in --> Prob_in --> Spike_in --> Spike_out --> Rate_out

aux1_rate_int = zeros(N_int,trial_dur_total_in_bins+2);
aux2_rate_int = zeros(N_int,trial_dur_total_in_bins+2);

decision_ac = zeros(N_int,trial_dur_total_in_bins+2);

spikes_target = cell(2,1); spikes_vis = cell(2,1); spikes_vest = cell(2,1); 
rate_int = cell(2,1); spikes_int = cell(2,1);
rate_lip = cell(2,1); spikes_lip = cell(2,1);

spikes_count_vis = cell(2,1); spikes_count_vest = cell(2,1); proba_count_int = cell(2,1); spikes_count_int = cell(2,1);

for m = 1:if_test_set+1 % Train and Test
    
    % Raw data
    spikes_target{m} = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    spikes_vis{m} = zeros(N_vis,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    spikes_vest{m} = zeros(N_vest,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    
    rate_int{m} = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    spikes_int{m} = zeros(N_int,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    rate_lip{m} = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    spikes_lip{m} = zeros(N_lip,trial_dur_total_in_bins,N_trial,length(unique_heading),length(unique_condition));
    
    % Spike count in sliding windows for Fisher information
    spikes_count_vis{m} = nan(N_vis, fisher_N_window, N_trial,length(unique_heading),length(unique_condition));
    spikes_count_vest{m} = nan(N_vest, fisher_N_window, N_trial,length(unique_heading),length(unique_condition));
    proba_count_int{m} = nan(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_condition));
    spikes_count_int{m} = nan(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_condition));
    if shuffle_flag==1
        spikes_count_int_shuff{m} = zeros(N_int, fisher_N_window, N_trial, length(unique_heading),length(unique_condition));
    end
    
end

% Output variables
RT = (trial_dur_total-stim_on_time) * ones(N_trial,length(unique_heading),length(unique_condition));

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
    w_int_vis(nn,:) = g_w_int_vis/N_vis *(exp(K_int_vis * (cos((prefs_int(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vis/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        + dc_w_int_vis/N_vis;   % Offset
    w_int_vest(nn,:) = g_w_int_vest/N_vest *(exp(K_int_vest * (cos((prefs_int(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vest/180*pi))... % Gaussian(theta_int - +/-90) * Sin(heading)
        + dc_w_int_vest/N_vest;   % Offset
    
    % -- Targ->LIP --
    w_lip_targ(nn,:) = g_w_lip_targ/N_int *(exp(K_lip_targ*(cos((prefs_int-prefs_int(nn))/360 *2*pi)-1)));  %  Target input
    
    % -- Int->Int --
    w_int_int(nn,:) = g_w_int_int/N_int*...   %  Integrator recurrent
        ((exp(K_int_int*(cos((prefs_int-prefs_int(nn))/360*2*pi)-1)))-...
        amp_I_int*(exp(K_int_I*(cos((prefs_int-prefs_int(nn))/360*2*pi)-1))))...
        + dc_w_int_int/N_int;
    
end

for nn=1:N_lip
    % -- Int->LIP --
    w_lip_int(nn,:) = g_w_lip_int/N_int*...
        ((exp(K_lip_int*(cos((prefs_lip-prefs_int(nn))/360*2*pi)-1)))-...
        amp_I_lip_int*(exp(K_lip_int_I*(cos((prefs_lip-prefs_int(nn))/360*2*pi)-1))))...
        + dc_w_lip_int/N_int;

    % -- LIP->LIP --
    w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...
        ((exp(K_lip_lip*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
        amp_I_lip*(exp(K_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
        + dc_w_lip_lip/N_lip;
end    

if if_debug
    figure(90); clf;
    set(gcf,'uni','norm','pos',[0.006       0.099       0.779       0.766],'name','Weights & Raster plot');
    
    subplot(4,3,1);
    imagesc(prefs_lip,prefs_lip,w_lip_lip'); colorbar; axis  xy; title('LIP->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{lip}');
    set(gca,'xtick',-180:90:180);

    subplot(4,3,4);
    imagesc(prefs_lip,prefs_int,w_lip_int'); colorbar; axis  xy; title('Int->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180);
    
    subplot(4,3,7);
    % surf(prefs_int,prefs_int,w_int_int');
    imagesc(prefs_int,prefs_int,w_int_int');    colorbar; axis xy; title('Int->Int');
    xlabel('\theta_{int}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    
    %     subplot(4,3,7);
    %     imagesc(prefs_int,prefs_vest,w_int_vest'); colorbar; axis  xy; title('Vest->Int');
    %     xlabel('\theta_{int}'); ylabel('\theta_{vest}'); ylim([-20 20]);
    %     set(gca,'xtick',-180:90:180);
    
    subplot(4,3,10);
    imagesc(prefs_int,prefs_vis,w_int_vis'); colorbar; axis  xy; title('Vis/Vest->Int');
    xlabel('\theta_{int}'); ylabel('\theta_{vis/vest}'); ylim([-20 20]);
    set(gca,'xtick',-180:90:180);
    
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
count = 1;

for cc = 1:length(unique_condition)
    for hh = 1:length(unique_heading)  % Motion directions
        
        % ==== Not trial- or time- dependent stuffs ===
        stim_this = unique_heading(hh);
        cond_this = unique_condition(cc);
        
        fprintf('cond = %g, heading = %g\n  ',cond_this, stim_this);
        
        % -- Sensory responses --
        % proba_vis: proba of firing of visual neurons in response to motion
        max_rate_vis = r_spont_vis + b_pref_vis * coherence;
        b_vis = r_spont_vis + b_null_vis * coherence;
        proba_vis = ((max_rate_vis-b_vis)*exp(K_vis*(cos((prefs_vis'-stim_this)/360*2*pi)-1))+b_vis)*dt;
        
        %     max_rate_vest = r_spont_vest + b_pref_vis * coherence;
        %     b_vis = r_spont_vis + b_null_vis * coherence;
        
        % Vestibular response does not depend on coherence.
        % Here I just set the vestibular activity similar to visual response under 'equivalent_coh' coh
        max_rate_vest = r_spont_vest + b_pref_vest * equivalent_conherence;
        b_vest = r_spont_vest + b_null_vest * equivalent_conherence;
        proba_vest = ((max_rate_vest-b_vest)*exp(K_vest*(cos((prefs_vest'-stim_this)/360*2*pi)-1))+b_vest)*dt;
        
        %         % -- Stimulus condition selection --
        %         if cond_this == 1
        %             proba_vis = 0*proba_vis; % Shut down visual activity (but there are still noisy spikes if select here)
        %         elseif cond_this ==2
        %             proba_vest = 0*proba_vest; % Shut down vestibular activity
        %         end
        
        if if_debug
            figure(91);clf;
            plot(prefs_vis,proba_vis/dt,'r');
            plot(prefs_vest,proba_vest/dt,'b');
            title('Sensory population tuning');
        end
        
        % proba_targ: proba of firing of target neurons in response to visual targets
        
        % pos_targ = stim_this + [0:360/N_targets:359.9];
        pos_targ = [-90 90]; % Always these two for fine discrimination task. HH2017
        
        proba_target_tuning = zeros(N_int,1);
        for nn=1:N_targets
            proba_target_tuning = proba_target_tuning + ((max_rate_target-b_target)*...
                exp(K_target*(cos((prefs_int'-pos_targ(nn))/360*2*pi)-1))+b_target)*dt;
        end
        proba_target_tuning = proba_target_tuning/(slope_norm_target*(N_targets/2)+dc_norm_target);   %  Divisive normalization
        
        
        % === Begin trials ===
        for mm = 1:if_test_set+1 % Train and test sets
            for tt = 1:N_trial % For each trial
                
                if round(tt/50)==tt/50
                    fprintf('%d  ',tt);
                end
                
                % == Generate input spike trains ==
                
                %                 % -- Target input spike train --
                %                 % Note that trial_dur_total now includes the pre-trial time.
                %                 % There is no need (actually erroneous?) to separate pre-trial and trial.
                %                 % Because in this case, the INTEGRATOR dynamics are missing during pre-trial,
                %                 % and the MT's noise term was also missing.
                %
                %                 % Ohh, I got it. We don't want to accumulate target information during pre-trial,
                %                 % so we deliberately ignore INTEGRATOR dynamics during that period (although a little bit weird)!!
                
                aux_proba_target = proba_target_tuning*ones(1,trial_dur_total_in_bins); % Expand along the time axis
                spikes_target{mm}(:,:,tt,hh,cc) = rand(N_int,trial_dur_total_in_bins)<(aux_proba_target);
                
                % -- Visual input spike train --
                % The stim_on_time controls when motion information comes in.
                % aux_proba_vis = proba_vis*[zeros(1,stim_on_time) ones(1,trial_dur_total-stim_on_time)]...
                %    + w_cov_vis*randn(N_vis,trial_dur_total);
                
                % Temporal dynamics added. HH2017
                aux_proba_vis = proba_vis*[zeros(1,stim_on_time_in_bins) vel*gain_vel_vis]...
                    + w_cov_vis*randn(N_vis,trial_dur_total_in_bins);
                
                % -- Vestibular ---
                % With the temporal gain of abs(acc)
                aux_proba_vest = proba_vest*[zeros(1,stim_on_time_in_bins) abs(acc)*gain_acc_vest]...
                    + w_cov_vest*randn(N_vest,trial_dur_total_in_bins);
                
                % -- Stimulus condition selection --
                if cond_this == 1
                    aux_proba_vis = 0*aux_proba_vis; % Shut down visual activity
                elseif cond_this ==2
                    aux_proba_vest = 0*aux_proba_vest; % Shut down vestibular activity
                end
                
                spikes_vis{mm}(:,1:trial_dur_total_in_bins,tt,hh,cc)= rand(N_vis,trial_dur_total_in_bins)<(aux_proba_vis);
                spikes_vest{mm}(:,1:trial_dur_total_in_bins,tt,hh,cc)= rand(N_vest,trial_dur_total_in_bins)<(aux_proba_vest);
                
                % -- Reset the attention for motion stimuli --
                att_gain_stim = 1; % Drop in attention to motion stimuli due to hitting whatever decision bound.
                
                %                 % == Override INTEGRATOR dynamics during pre-trial by setting the spikes manually (weird thing) ==
                % !! Now we have the dedicated INTEGRATION layer so we can use all the dynamics from the very beginning. HH20170317
                
                %                 rate_int{mm}(:,1:stim_on_time_in_bins,tt,hh,cc) = aux_proba_target(:,1:stim_on_time_in_bins) / dt;
                %                 spikes_int{mm}(:,1:stim_on_time_in_bins,tt,hh,cc) = spikes_target{mm}(:,1:stim_on_time_in_bins,tt,hh,cc);
                %
                %                 aux1_rate_int(:,stim_on_time_in_bins) = rate_int{mm}(:,stim_on_time_in_bins,tt,hh,cc);
                %                 aux2_rate_int(:,stim_on_time_in_bins) = 0;
                
                % == Real network dynamics after stimlus onset ==
                %                 k = stim_on_time_in_bins; % Time begins at stim_on_time
                k = 1;  % Time begins at the very beginning
                
                while k<=trial_dur_total_in_bins-1
                    
                    %                     % -- I don't understand yet why there are two variables aux1 and aux2 --
                    %                     aux1_rate_int(:,k+1) = (1-dt/time_const_vis)*aux1_rate_int(:,k)...   %  Self dynamics.  in Hz!
                    %                         + 1/time_const_vis*(...
                    %                         att_gain_stim * w_int_vis * spikes_vis{mm}(:,k,tt,hh,cc)...     %  Visual input
                    %                         + att_gain_stim * w_int_vest * spikes_vest{mm}(:,k,tt,hh,cc)...     % Vestibular input
                    %                         ... % + att_gain_targ * w_int_targ * spikes_target{mm}(:,k,tt,hh,cc)...  % No longer here. HH20170317
                    %                         );             %  Target visual input (att_gain2 has been set to 0)
                    %
                    %                     aux2_rate_int(:,k+1) = (1-dt/time_const_int)*aux2_rate_int(:,k)+...  %  Self dynamics.  in Hz!
                    %                         +1/time_const_int*(w_int_int*spikes_int{mm}(:,k,tt,hh,cc));  %  INTEGRATOR recurrent. in Hz!
                    %
                    %                     rate_int{mm}(:,k+1,tt,hh,cc) = aux1_rate_int(:,k+1) + aux2_rate_int(:,k+1) + bias_int;  %  Bias term: b_out
                    %
                    
                    % -- Update INTEGRATOR layer --
                    rate_int{mm}(:,k+1,tt,hh,cc) = bias_int + (1-dt/time_const_int)*rate_int{mm}(:,k,tt,hh,cc)...   %  Self dynamics.  in Hz!
                        + 1/time_const_int * (...
                              w_int_int * spikes_int{mm}(:,k,tt,hh,cc)... %  INTEGRATOR recurrent
                            + att_gain_stim * w_int_vis * spikes_vis{mm}(:,k,tt,hh,cc)...     %  Visual input
                            + att_gain_stim * w_int_vest * spikes_vest{mm}(:,k,tt,hh,cc)...     % Vestibular input
                        ... % + att_gain_targ * w_int_targ * spikes_target{mm}(:,k,tt,hh,cc)...  % No longer here. HH20170317
                        );            
                    
                    % -- Update LIP layer --
                    rate_lip{mm}(:,k+1,tt,hh,cc) = bias_lip + (1-dt/time_const_lip)*rate_lip{mm}(:,k,tt,hh,cc)...   %  Self dynamics.  in Hz!
                        + 1/time_const_lip * (...
                              w_lip_lip * spikes_lip{mm}(:,k,tt,hh,cc)...  %  LIP recurrent
                            + att_gain_targ * w_lip_targ * spikes_target{mm}(:,k,tt,hh,cc)...  % Moved here. HH20170317
                            + w_lip_int * spikes_int{mm}(:,k,tt,hh,cc)... %  INTEGRATOR->LIP
                        );            
                    
                    
                    % -- Turn rate to binary spike for the next step --
                    spikes_int{mm}(:,k+1,tt,hh,cc)  = (rand(N_int,1) < dt*(rate_int{mm}(:,k+1,tt,hh,cc)-threshold_int));
                    spikes_lip{mm}(:,k+1,tt,hh,cc)  = (rand(N_lip,1) < dt*(rate_lip{mm}(:,k+1,tt,hh,cc)-threshold_lip));
                    
                    % -- Variable used for stopping the integration --
                    
                    %{
                    if mm == 1 % Only apply to train set
                        decision_ac(:,k+1) = rate_int{mm}(:,k+1,tt,hh,cc);
                        %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
                        %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
                        
                        % -- Termination --
                        if ((range(if_bounded*decision_ac(:,k+1)) > decis_thres) && (RT(tt,hh,cc)==0)) || (k==trial_dur_total_in_bins-1)
                            RT(tt,hh,cc) = k;
                            if k>fisher_wind_size-1   %  Enough data for Fisher info of at least one window
                                last_spikes__last_Fisher_win(:,count) = sum(spikes_int{1}(:,k-fisher_wind_size+1:k)')';
                            else
                                last_spikes__last_Fisher_win(:,count) = sum(spikes_int{1}(:,1:k)')';
                            end
                            last_proba(:,count) = rate_int{mm}(:,k,tt,hh,cc);
                            k=trial_dur_total_in_bins;  %  Terminate trial immediately
                        end
                    end
                    %}
                    
                    if mm == 1 % Only apply to train set
                        decision_ac(:,k+1) = rate_int{mm}(:,k+1,tt,hh,cc);
                        %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
                        %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
                        
                        % -- Termination --
                        if if_bounded && att_gain_stim == 1 && (f_bound(decision_ac(:,k+1)) > decis_thres(unique_condition(cc)))
                            % Set the attention for motion stimuli to zero
                            att_gain_stim = att_gain_stim_after_hit_bound(unique_condition(cc));
                            RT(tt,hh,cc) = (k-stim_on_time_in_bins)*dt;
                            % last_proba(:,count) = rate_int{mm}(:,k,tt,hh,cc);
                        end
                    end
                    
                    k=k+1;
                end  % of network dynamics
                
                %{
            %%% collect spike and firing probabilities for Tin and Tout neurons
            
            proba_out_av1 = proba_out_av1 + rate_int{1};
            proba_out50(count,:) = rate_int{1}(50,:);
            proba_out75(count,:) = rate_int{1}(75,:);
            proba_out100(count,:) = rate_int{1}(100,:);
            spikes_out50(count,:) = spikes_int{1}(50,:);
            spikes_out75(count,:) = spikes_int{1}(75,:);
            spikes_out100(count,:) = spikes_int{1}(100,:);
            spikes_in50(count,:) = spikes_vis{1}(50,:);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% training and testing data. Spike counts and average probabilities.
            %% size nu_unit x 2 N_trial (N_trial for each of the 2 positions)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spikes_in_trial11(:,count) = sum(spikes_vis{1}')';
            spikes_in_trial12(:,count) = sum(spikes_vis{2}')';
            
            spikes_out_trial01(:,count) = sum(spikes_int{1}')';
            spikes_out_trial02(:,count) = sum(spikes_int{2}')';
            
            proba_out_trial01(:,count) = sum(rate_int{1}(:,1:trial_dur_total)')'/fisher_wind_size;
            proba_out_trial02(:,count) = sum(rate_int{2}(:,1:trial_dur_total)')'/fisher_wind_size;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% collect counts in small time windows
            %% m=1 training data, m=2 testing data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for n=1:fisher_N_window
                proba_count_int{m,n}(:,count)= ...
                    sum(rate_int{m}(:,fisher_wind_offset+fisher_wind_size*(n-1)+1:fisher_wind_size*n)')'/fisher_wind_size;
                spikes_count_int{m,n}(:,count)= ...
                    sum(spikes_int{m}(:,fisher_wind_offset+fisher_wind_size*(n-1)+1:fisher_wind_size*n)')';
                spikes_count_vis{m,n}(:,count)= ...
                    sum(spikes_vis{m}(:,1:fisher_wind_offset+fisher_wind_size*n)')';
            end
                %}
                count = count+1;
            end % of trial
        end % of train/test
        
        fprintf('\n');
        
    end % of heading
end % of stimulus condition


%% == Plotting example network activities ==
to_plot_heading = length(unique_heading);
to_plot_cond = length(unique_condition);

rate_real_vis_aver = nanmean(spikes_vis{1}(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
rate_real_vest_aver = nanmean(spikes_vest{1}(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
rate_expected_int_aver = nanmean(rate_int{1}(:,:,:,to_plot_heading,to_plot_cond),3);
rate_expected_lip_aver = nanmean(rate_lip{1}(:,:,:,to_plot_heading,to_plot_cond),3);
rate_real_int_aver = nanmean(spikes_int{1}(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
rate_real_lip_aver = nanmean(spikes_lip{1}(:,:,:,to_plot_heading,to_plot_cond),3)/dt;


%{
%%  ====== Animation ======
figure(1001); clf

for ttt = 1:length(ts)
    % LIP layer
    subplot(3,1,1);
    plot(prefs_lip,rate_expected_lip_aver(:,ttt),'r');  hold on;
    plot(prefs_lip,rate_real_lip_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',min(prefs_lip:90:max(prefs_lip)));
    xlim([min(prefs_lip) max(prefs_lip)]);  ylim([-20 50]);
    title(ttt*dt*1000);

    % INTEGRATOR layer
    subplot(3,1,2);
    plot(prefs_int,rate_expected_int_aver(:,ttt),'r');  hold on;
    plot(prefs_int,rate_real_int_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',min(prefs_int):90:max(prefs_int));
    xlim([min(prefs_int) max(prefs_int)]);  ylim([-20 50]);
    
    % Visual layer
    subplot(3,1,3);
    plot(prefs_vis,rate_real_vis_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',min(prefs_vis):90:max(prefs_vis));
    xlim([min(prefs_vis) max(prefs_vis)]);  ylim([-20 50]);
    
    drawnow;
end


%}

%% ====== Raster plot ======
if if_debug
    figure(90);
    
    subplot(4,3,[2 3]);
    imagesc(ts,prefs_lip,rate_real_lip_aver*dt); axis xy;
    ylabel('LIP');
    title(sprintf('Firing prob. for each bin, averaged over %g trials, cond = %g, heading = %g',...
        N_trial,unique_condition(to_plot_cond),unique_heading(to_plot_heading)));  colorbar

    subplot(4,3,[5 6]);
    imagesc(ts,prefs_int,rate_real_int_aver*dt); axis xy; colorbar
    ylabel('INTEGRATOR');
    
    subplot(4,3,[8 9]);
    imagesc(ts,prefs_vest,rate_real_vest_aver*dt); axis xy; colorbar
    ylabel('Vest');
    
    subplot(4,3,[11 12]);
    imagesc(ts,prefs_vis,rate_real_vis_aver*dt); axis xy; colorbar
    ylabel('Vis');
    
end
%}

% ====== Example Pref and Null traces ======
%%{
set(figure(1002),'name','Example PSTHs'); clf;
set(gcf,'uni','norm','pos',[0.003       0.039       0.555       0.878]);

for cc = 1:length(unique_condition)
    % --- LIP ---
    subplot(3,2,2*unique_condition(cc)-1);
    
    for trialtoplot = 1:round(N_trial/10):N_trial
        plot(ts,rate_lip{1}(right_targ_ind,:,trialtoplot,to_plot_heading,cc),'color',colors(unique_condition(cc),:),'linewid',2); hold on;
        plot(ts,rate_lip{1}(left_targ_ind,:,trialtoplot,to_plot_heading,cc),'--k','linewid',1);
    end
    plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
    axis tight;
    
    %     if if_bounded
    %         plot(xlim,[decis_thres decis_thres],'k--');
    %         ylim([min(ylim),decis_thres*1.1]);
    %     end
    
    title(sprintf('rate\\_lip, heading = %g, coh = %g',unique_heading(to_plot_heading),coherence));
    
    % --- Int ---
    subplot(3,2,2*unique_condition(cc));
    
    for trialtoplot = 1:round(N_trial/10):N_trial
        plot(ts,rate_int{1}(right_targ_ind,:,trialtoplot,to_plot_heading,cc),'color',colors(unique_condition(cc),:),'linewid',2); hold on;
        plot(ts,rate_int{1}(left_targ_ind,:,trialtoplot,to_plot_heading,cc),'--k','linewid',1);
    end
    plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
    axis tight;
    
    %     if if_bounded
    %         plot(xlim,[decis_thres decis_thres],'k--');
    %         ylim([min(ylim),decis_thres*1.1]);
    %     end
    
    title(sprintf('rate\\_int, heading = %g, coh = %g',unique_heading(to_plot_heading),coherence));
    
    
    %     subplot(3,2,2);
    %     plot(ts,nanmean(rate_lip{1}(right_targ_ind,:,:,to_plot_heading,cc),3)...
    %         -nanmean(rate_lip{1}(left_targ_ind,:,:,to_plot_heading,cc),3),'color',colors(unique_condition(cc),:),'linew',2);
    %     hold on;
    %     title(sprintf('averaged of all %g trials (correct + wrong), pref-null',N_trial));
    %     axis tight;
end

%}

% ====== Behavior performance ======
%%{
figure(99); clf; hold on; axis([-8 8 0 1]);
set(gcf,'name','Behavior','uni','norm','pos',[0.32       0.071       0.357        0.83]);

% Psychometric
subplot(2,1,1);
for cc = 1:length(unique_condition)
    % Make decisions
    rate_int_at_decision = squeeze(rate_lip{1}(:,end,:,:,cc)); % Using lip
    
    [~, pos_max_rate_int_at_decision] = max(rate_int_at_decision,[],1);
    choices{cc} = prefs_int(shiftdim(pos_max_rate_int_at_decision)) >= 0; % 1 = rightward, 0 = leftward
    
    % I just flip the psychometric curve to the negative headings
    psychometric = [unique_heading' sum(choices{cc}(:),1)'/N_trial];
    fake_psy = flipud(psychometric(unique_heading>0,:));
    fake_psy(:,1) = -fake_psy(:,1);
    fake_psy(:,2) = 1-fake_psy(:,2);
    psychometric = [fake_psy ; psychometric];
    
    plot(psychometric(:,1),psychometric(:,2),'o','color',colors(unique_condition(cc),:),'markerfacecolor',colors(unique_condition(cc),:),'markersize',13); % Psychometric
    hold on;
    xxx = -max(unique_heading):0.1:max(unique_heading);
    
    [bias, threshold] = cum_gaussfit_max1(psychometric);
    psycho(cc,:) = [bias,threshold];
    plot(xxx,normcdf(xxx,bias,threshold),'-','color',colors(unique_condition(cc),:),'linewid',4);
    
    set(text(min(xlim),0.4+0.1*cc,sprintf('threshold = %g\n',threshold)),'color',colors(unique_condition(cc),:));
end

if length(unique_condition) == 3
    pred_thres = (psycho(1,2)^(-2)+psycho(2,2)^(-2))^(-1/2);
    set(text(min(xlim),0.4+0.1*(cc+1),sprintf('pred = %g\n',pred_thres)),'color','k');
end

% Chronometric
subplot(2,1,2);
for cc = 1:length(unique_condition)
    % Make decisions
    RT_this_cc_mean = mean(RT(:,:,cc),1);
    RT_this_cc_se = std(RT(:,:,cc),[],1);
    
    errorbar(unique_heading+(cc-2)*0.2,RT_this_cc_mean,RT_this_cc_se,'o-','color',colors(unique_condition(cc),:),...
        'markerfacecolor',colors(unique_condition(cc),:),'markersize',10,'linew',2);
    hold on;
end
plot(vel/max(vel)*max(xlim)/5,t_motion,'k--');
title(sprintf('RT, all trials (correct + wrong), Bound = %s',num2str(decis_thres)));

%}

%% ===== Averaged PSTH, different angles =====
%%{

PSTH_condition_absheading_choice = nan(length(ts),length(unique_condition),length(unique_heading),2);

set(figure(435),'name','PSTH'); clf;
set(gcf,'uni','norm','pos',[0.326       0.041       0.668       0.681]);
% abs_unique_heading = unique(abs(unique_heading));
sub_x = fix(sqrt(length(unique_heading)));
sub_y = ceil(length(unique_heading)/sub_x);
y_max = -inf; y_min = inf;

for hh = 1:length(unique_heading)
    
    % == Raw PSTH (stim type), heading, correct only
    figure(435); subplot(sub_x,sub_y,hh); hold on;
    
    for cc = 1:length(unique_condition)
        
        % I assume "PREF" = "LEFT" here. In this case, with headings [0 1 2 4 8], I
        % plot the pref=90 neuron and the pref=-90 neuron in correct trials to simulate
        % pref and null activity in correct only trials for each abs(heading) in my experiments.
        select_correct = (choices{cc}(:,hh) == (unique_heading(hh)>=0));
        
        % Although note here the PREF and NULL activity become correlated.
        PSTH_condition_absheading_choice(:,cc,hh,1) = squeeze(nanmean(rate_lip{1}(right_targ_ind,:,select_correct,hh,cc),3));
        PSTH_condition_absheading_choice(:,cc,hh,2) = squeeze(nanmean(rate_lip{1}(left_targ_ind,:,select_correct,hh,cc),3));
        
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,1),'color',colors(unique_condition(cc),:),'linew',2.5);
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,2),'--','color',colors(unique_condition(cc),:),'linew',2.5);
        title(sprintf('abs(heading) = %g, %g trials, correct only',abs(unique_heading(hh)),N_trial));
    end
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    axis tight;
    
    %     if if_bounded
    %         plot(xlim,[decis_thres decis_thres],'k--');
    %         ylim([min(ylim),decis_thres*1.1]);
    %     end
    
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    
end

set(findobj(435,'type','axes'),'ylim',[y_min y_max]);

% =====  Delta-PSTH (stim type), heading, correct only ====
set(figure(436),'name','Delta PSTH'); clf;
set(gcf,'uni','norm','pos',[0.326       0.041       0.668       0.681]);
y_max = -inf; y_min = inf;

for hh = 1:length(unique_heading)
    
    figure(436);  subplot(sub_x,sub_y,hh); hold on;
    
    for cc = 1:length(unique_condition)
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,1)-PSTH_condition_absheading_choice(:,cc,hh,2),'color',colors(unique_condition(cc),:),'linew',2.5);
        title(sprintf('abs(heading) = %g, %g trials, correct only',abs(unique_heading(hh)),N_trial));
    end
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    axis tight;
end
set(findobj(436,'type','axes'),'ylim',[y_min y_max]);

%}

toc


