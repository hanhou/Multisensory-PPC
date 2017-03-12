%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
% Modified by HH 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear
clear global

load h

global Y_targ X_train;
rand('state',sum(100*clock))

% %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Sizes ===
% Input layers
N_vis = 200; % Visual motion signal
prefs_vis = 0:360/N_vis:359.9999;

N_targets = 2; % Target input

% Output layers
N_lip = 100;
prefs_lip = 0:360/N_lip:359.9999;

% === Time ===
dt = 1e-3; % Size of time bin in seconds
trial_dur_total = 300; %% in time bins (including pre_trial_dur)
stim_on_time = 50; % Time of motion start, in time bins.
N_trial = 100; % For each condition

if stim_on_time>trial_dur_total
    error('pre_trial_dur cannot be longer than trial_dur');
end

decis_flag = 1; %% decis_flag=1 means that the trial stops at the bound (reaction time version)
decis_thres = 58; %% bound height

% === Stimuli ===
stim_positions = [90];
coherence = 30;

% Parameters for MT from Mazurek and Shadlen, except width.
% K_cov_mt and var_mt controls the correlations.
%  << cov_mt(j,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(j))/360 *2*pi)-1));
%     w_mt = real(sqrtm(cov_mt));
%     aux_proba_in = proba_in + w_mt*randn(N_vis,1); >>
r_spont_vis = 20;
b_pref_vis = 0.4;
b_null_vis = -0.2;
K_vis = 4;
K_cov_vis = 2;
var_vis = 1e-5;

% Parameter of the activity evoked by the visual targets
max_rate_target = 42; %% in Hz
b_target = 0;
K_target = 4;
slope_norm_target=0.5;
dc_norm_target = 0.5;


% === Network configuration ===

% -- Time constant for integration
time_const_lip = 1000e-3; % in s
time_const_vis  = time_const_lip;

% -- Visual to LIP
g_w_lip_vis =65;
K_lip_vis = 5;

% --- Targets to LIP
g_w_lip_targ= 40;
K_lip_targ= 5;
att_gain2_targ = 0; % Drop in attention to visual target once motion stimulus appears.

% --- Recurrent connectivity in LIP
g_w_lip_lip = 21;
K_lip_lip = 10;
dc_w_lip_lip = -0.17;

amp_I = 0.0;  % << What's this? HH >>
K_lip_I = 2;

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

aux1_rate_lip = zeros(N_lip,trial_dur_total+2);
aux2_rate_lip = zeros(N_lip,trial_dur_total+2);

spikes_in_trial11 = zeros(N_vis,N_trial*2);
spikes_in_trial12 = zeros(N_vis,N_trial*2);

decision_ac = zeros(N_lip,trial_dur_total+2);


spikes_target = cell(2,1); spikes_vis = cell(2,1); rate_lip = cell(2,1); spikes_lip = cell(2,1);
spikes_count_vis = cell(2,1); proba_count_lip = cell(2,1); spikes_count_lip = cell(2,1);

for m = 1:2 % Train and Test
    
    % Raw data
    spikes_target{m} = nan(N_lip,trial_dur_total,N_trial,length(stim_positions));
    spikes_vis{m} = nan(N_vis,trial_dur_total,N_trial,length(stim_positions));
    rate_lip{m} = nan(N_lip,trial_dur_total,N_trial,length(stim_positions));
    spikes_lip{m} = nan(N_lip,trial_dur_total,N_trial,length(stim_positions));
    
    % Spike count in sliding windows for Fisher information
    spikes_count_vis{m} = zeros(N_vis, fisher_N_window, N_trial,length(stim_positions));
    proba_count_lip{m} = zeros(N_lip, fisher_N_window, N_trial, length(stim_positions));
    spikes_count_lip{m} = zeros(N_lip, fisher_N_window, N_trial, length(stim_positions));
    if shuffle_flag==1
        spikes_count_lip_shuff{m} = zeros(N_lip, fisher_N_window, N_trial, length(stim_positions));
    end
    
end

% Output variables
proba_out01__pretrial  = zeros(N_lip,stim_on_time);
proba_out02__pretrial  = zeros(N_lip,stim_on_time);
spikes_out01__pretrial_train = zeros(N_lip,stim_on_time);
spikes_out02__pretrial_test = zeros(N_lip,stim_on_time);

proba_out_av0__pretrial  = zeros(N_lip,stim_on_time);
proba_out_av1  = zeros(N_lip,trial_dur_total);

spikes_out_trial01 = zeros(N_lip,N_trial*2);
proba_out_trial01  = zeros(N_lip,N_trial*2);
spikes_out_trial02 = zeros(N_lip,N_trial*2);
proba_out_trial02  = zeros(N_lip,N_trial*2);

proba_out50_pre=zeros(N_trial*2,stim_on_time);
proba_out1_pre=zeros(N_trial*2,stim_on_time);
proba_out50=zeros(N_trial*2,trial_dur_total);
proba_out1=zeros(N_trial*2,trial_dur_total);

RT = zeros(N_trial,length(stim_positions));

spikes_out50 = zeros(N_trial*2,trial_dur_total);
spikes_out75 = zeros(N_trial*2,trial_dur_total);
spikes_out100 = zeros(N_trial*2,trial_dur_total);
last_spikes = zeros(N_lip,2*N_trial);
last_proba = zeros(N_lip,2*N_trial);


% << Network connections >>
w_lip_vis = zeros(N_lip,N_vis);
w_lip_targ = zeros(N_lip,N_lip); % Not N_target, but N_lip (supposing the same number as the lip neurons)
w_lip_lip = zeros(N_lip,N_lip);

for nn=1:N_lip
    w_lip_vis(nn,:) = g_w_lip_vis/N_vis *(exp(K_lip_vis*(cos((prefs_vis-prefs_lip(nn))/360 *2*pi)-1)));  % << MT input >>
    w_lip_targ(nn,:) = g_w_lip_targ/N_vis *(exp(K_lip_targ*(cos((prefs_lip-prefs_lip(nn))/360 *2*pi)-1)));  % << Target input >>
    w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...   % << LIP recurrent >>
        ((exp(K_lip_lip*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
        amp_I*(exp(K_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
        + dc_w_lip_lip;
end

% << MT correlation matrix >>
cov_vis = zeros(N_vis,N_vis);
for nn=1:N_vis
    cov_vis(nn,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(nn))/360 *2*pi)-1));
end
w_cov_vis = real(sqrtm(cov_vis));


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;

for ss = 1:length(stim_positions)  % Motion directions
    
    % ==== Not trial- or time- dependent stuffs ===
    stim_this = stim_positions(ss);
    fprintf('%dth  pos: %4.2f\n  ',ss, stim_this);
    
    % proba_vis: proba of firing of visual neurons in response to motion
    max_rate_vis = r_spont_vis + b_pref_vis * coherence;
    b_vis = r_spont_vis + b_null_vis * coherence;
    proba_vis_tuning = ((max_rate_vis-b_vis)*exp(K_vis*(cos((prefs_vis'-stim_this)/360*2*pi)-1))+b_vis)*dt;
    
    % proba_targ: proba of firing of target neurons in response to visual targets
    pos_targ = stim_this + [0:360/N_targets:359.9];
    proba_target_tuning = zeros(N_lip,1);
    for nn=1:N_targets
        proba_target_tuning = proba_target_tuning + ((max_rate_target-b_target)*...
            exp(K_target*(cos((prefs_lip'-pos_targ(nn))/360*2*pi)-1))+b_target)*dt;
    end
    proba_target_tuning = proba_target_tuning/(slope_norm_target*(N_targets/2)+dc_norm_target);   % << Divisive normalization >>
    
    
    % === Begin trials ===
    for mm = 1:2 % Train and test sets
        for tt = 1:N_trial % For each trial
            
            if round(tt/50)==tt/50
                fprintf('%d  ',tt);
            end
            
            % == Generate input spike trains ==
            
            % -- Target input spike train --
            % Note that trial_dur_total now includes the pre-trial time.
            % There is no need (actually erroneous?) to separate pre-trial and trial.
            % Because in this case, the LIP dynamics are missing during pre-trial,
            % and the MT's noise term was also missing.
            
            % Ohh, I got it. We don't want to accumulate target information during pre-trial,
            % so we deliberately ignore LIP dynamics during that period (although a little bit weird)!!
            
            aux_proba_target = proba_target_tuning*ones(1,trial_dur_total); % Expand along the time axis
            spikes_target{mm}(:,:,tt,ss) = rand(N_lip,trial_dur_total)<(aux_proba_target);
            
            % -- Visual input spike train --
            % The stim_on_time controls when motion information comes in.
            aux_proba_vis = proba_vis_tuning*[zeros(1,stim_on_time) ones(1,trial_dur_total-stim_on_time)]...
                + w_cov_vis*randn(N_vis,trial_dur_total);
            spikes_vis{mm}(:,1:trial_dur_total,tt,ss)= rand(N_vis,trial_dur_total)<(aux_proba_vis);
            
            % -- Vestibular: too be added ---
            
            % == Override LIP dynamics during pre-trial by setting the spikes manually (weird thing) ==
            rate_lip{mm}(:,1:stim_on_time,tt,ss) = aux_proba_target(:,1:stim_on_time) / dt;
            spikes_lip{mm}(:,1:stim_on_time,tt,ss) = spikes_target{mm}(:,1:stim_on_time,tt,ss);
            
            aux1_rate_lip(:,stim_on_time) = rate_lip{mm}(:,stim_on_time,tt,ss);
            aux2_rate_lip(:,stim_on_time) = 0;
            
            % == Real network dynamics after stimlus onset ==
            k = stim_on_time; % Time begins at stim_on_time
            
            while k<=trial_dur_total-1
                
                % -- I don't understand yet why there are two variables aux1 and aux2 --
                aux1_rate_lip(:,k+1) = (1-dt/time_const_vis)*aux1_rate_lip(:,k)...   % << Self dynamics >> in Hz!
                    +1/time_const_vis*(w_lip_vis*spikes_vis{mm}(:,k,tt,ss)...     % << MT motion input >>
                    +w_lip_targ*att_gain2_targ*spikes_target{mm}(:,k,tt,ss));             % << Target visual input (att_gain2 has been set to 0) >>
                
                aux2_rate_lip(:,k+1) = (1-dt/time_const_lip)*aux2_rate_lip(:,k)+...
                    +1/time_const_lip*(w_lip_lip*spikes_lip{mm}(:,k,tt,ss));  % << LIP recurrent >> in Hz!
                
                rate_lip{mm}(:,k+1,tt,ss) = aux1_rate_lip(:,k+1) + aux2_rate_lip(:,k+1) + bias_lip;  % << Bias term: b_out >>
                
                % -- Note that proba_out is in Hz!! --
                spikes_lip{mm}(:,k+1,tt,ss)  = (rand(N_lip,1) < dt*(rate_lip{mm}(:,k+1,tt,ss)-threshold_lip));
                
                % -- Variable used for stopping the integration --
                
                if mm == 1 % Only apply to train set
                    decision_ac(:,k+1) = rate_lip{mm}(:,k+1,tt,ss);
                    %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
                    %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
                    
                    % -- Termination --
                    if ((max(decis_flag*decision_ac(:,k+1)) > decis_thres) && (RT(tt,ss)==0)) || (k==trial_dur_total-1)
                        RT(tt,ss) = k;
                        if k>fisher_wind_size-1   % << Enough data for Fisher info of at least one window >>
                            last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,k-fisher_wind_size+1:k)')';
                        else
                            last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,1:k)')';
                        end
                        last_proba(:,count) = rate_lip{mm}(:,k,tt,ss);
                        k=trial_dur_total;  % << Terminate trial immediately>>
                    end
                end
                k=k+1;
            end  % of network dynamics
            
        %{
            %%% collect spike and firing probabilities for Tin and Tout neurons
            
            proba_out_av1 = proba_out_av1 + rate_lip{1};
            proba_out50(count,:) = rate_lip{1}(50,:);
            proba_out75(count,:) = rate_lip{1}(75,:);
            proba_out100(count,:) = rate_lip{1}(100,:);
            spikes_out50(count,:) = spikes_lip{1}(50,:);
            spikes_out75(count,:) = spikes_lip{1}(75,:);
            spikes_out100(count,:) = spikes_lip{1}(100,:);
            spikes_in50(count,:) = spikes_vis{1}(50,:);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% training and testing data. Spike counts and average probabilities.
            %% size nu_unit x 2 N_trial (N_trial for each of the 2 positions)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spikes_in_trial11(:,count) = sum(spikes_vis{1}')';
            spikes_in_trial12(:,count) = sum(spikes_vis{2}')';
            
            spikes_out_trial01(:,count) = sum(spikes_lip{1}')';
            spikes_out_trial02(:,count) = sum(spikes_lip{2}')';
            
            proba_out_trial01(:,count) = sum(rate_lip{1}(:,1:trial_dur_total)')'/fisher_wind_size;
            proba_out_trial02(:,count) = sum(rate_lip{2}(:,1:trial_dur_total)')'/fisher_wind_size;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% collect counts in small time windows
            %% m=1 training data, m=2 testing data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for n=1:fisher_N_window
                proba_count_lip{m,n}(:,count)= ...
                    sum(rate_lip{m}(:,fisher_wind_offset+fisher_wind_size*(n-1)+1:fisher_wind_size*n)')'/fisher_wind_size;
                spikes_count_lip{m,n}(:,count)= ...
                    sum(spikes_lip{m}(:,fisher_wind_offset+fisher_wind_size*(n-1)+1:fisher_wind_size*n)')';
                spikes_count_vis{m,n}(:,count)= ...
                    sum(spikes_vis{m}(:,1:fisher_wind_offset+fisher_wind_size*n)')';
            end
        %}
            count = count+1;
        end % of trial
    end % of train/test
    
    fprintf('\n');
end % of stimulus position



%% == Plotting example network activities ==
%%{
figure(1001); clf

rate_real_vis_aver = nanmean(spikes_vis{1}(:,:,:,1),3)/dt;

rate_expected_lip_aver = nanmean(rate_lip{1}(:,:,:,1),3);
rate_real_lip_aver = nanmean(spikes_lip{1}(:,:,:,1),3)/dt;

for ttt = 1:trial_dur_total
    % LIP layer
    subplot(2,1,1);
    plot(prefs_lip,rate_expected_lip_aver(:,ttt),'r');  hold on;
    plot(prefs_lip,rate_real_lip_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',0:90:360);
    xlim([0 360]); ylim([-10 100]);
    title(ttt);
    
    % Visual layer
    subplot(2,1,2);
    plot(prefs_vis,rate_real_vis_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',0:90:360);
    xlim([0 360]); ylim([-10 100]);
    
    drawnow;
end


%}

%%{
figure(1002); clf;
to_plot_stim = 1;
[~,pref_ind] = min(abs(prefs_lip - stim_positions(to_plot_stim)));
[~,null_ind] = min(abs(prefs_lip - (stim_positions(to_plot_stim)+180)));

for cc = 1:10:N_trial
    plot(rate_lip{1}(pref_ind,:,cc,to_plot_stim),'r');
    hold on; plot(rate_lip{1}(null_ind,:,cc,to_plot_stim),'b--');
end
plot(xlim,[decis_thres decis_thres],'k--');
%}





