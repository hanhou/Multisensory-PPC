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
N_vis = 100; % Visual motion signal
prefs_vis = 0:180/N_vis:179.9999;

N_targets = 2; % Target input

% Output layers
N_lip = 100;
prefs_lip = 0:180/N_lip:179.9999;

% === Time ===
dt = 1e-3; % Size of time bin in seconds
trial_dur_total = 300; %% in time bins (including pre_trial_dur)
stim_on_time = 50; % Time of motion start, in time bins.
N_trial = 100; % For each condition

if stim_on_time>trial_dur_total
    error('pre_trial_dur cannot be longer than trial_dur');
end

decis_flag = 0; %% decis_flag=1 means that the trial stops at the bound (reaction time version)
decis_thres = 58; %% bound height

% === Stimuli ===
stim_positions = [90 92];
coherence = 6.4;

% Parameters for MT from Mazurek and Shadlen, except width.
% K_cov_mt and var_mt controls the correlations.
%  << cov_mt(j,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(j))/180 *2*pi)-1));
%     w_mt = real(sqrtm(cov_mt));
%     aux_proba_in = proba_in + w_mt*randn(N_vis,1); >>
r_spont_vis = 20;
b_pref_vis = 0.4;
b_null_vis = -0.2;
K_vis = 4;
K_cov_vis = 2;
var_vis = 1e-5;

% Parameter of the activity evoked by the visual targets
max_rate_targ = 42; %% in Hz
b_targ = 0;
K_targ = 4;
slope_norm_targ=0.5;
dc_norm_targ = 0.5;


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
wind_size = 50;
wind_offset = 0;
N_window = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization (I moved here)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I reorganized the variable naming and storage conventions. HH20170302
% All data have the size of (neurons, time, trial, stimuli){train/test}
% Three types of data are used:
% 1. firing rate in Hz (rate_xxx),
% 2. firing probability for each bin (proba_xxx),
% 3. spike or not for each bin (spikes_xxx)
% The network dynamics are computed like: Rate_in --> Prob_in --> Spike_in --> Spike_out --> Rate_out

aux1_proba1 = zeros(N_lip,trial_dur_total+2);
aux2_proba1 = zeros(N_lip,trial_dur_total+2);

spikes_in_trial11 = zeros(N_vis,N_trial*2);
spikes_in_trial12 = zeros(N_vis,N_trial*2);

decision_ac = zeros(N_lip,trial_dur_total+2);

aux1_proba2 = zeros(N_lip,trial_dur_total+2);
aux2_proba2 = zeros(N_lip,trial_dur_total+2);


spikes_target = cell(2,1); spikes_vis = cell(2,1); proba_lip = cell(2,1); spikes_lip = cell(2,1);
spikes_count_vis = cell(2,1); proba_count_lip = cell(2,1); spikes_count_lip = cell(2,1);

for m = 1:2 % Train and Test
    
    % Raw data
    spikes_target{m} = zeros(N_lip,trial_dur_total,N_trial,length(stim_positions));
    
    
    spikes_vis{m} = zeros(N_vis,trial_dur_total,N_trial,length(stim_positions));
    proba_lip{m} = zeros(N_lip,trial_dur_total,N_trial,length(stim_positions));
    spikes_lip{m} = zeros(N_lip,trial_dur_total,N_trial,length(stim_positions));
    
    % Spike count in sliding windows for Fisher information
    spikes_count_vis{m} = zeros(N_vis, N_window, N_trial,length(stim_positions));
    proba_count_lip{m} = zeros(N_lip, N_window, N_trial, length(stim_positions));
    spikes_count_lip{m} = zeros(N_lip, N_window, N_trial, length(stim_positions));
    if shuffle_flag==1
        spikes_count_lip_shuff{m} = zeros(N_lip, N_window, N_trial, length(stim_positions));
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

RT = zeros(N_trial*2,1);

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
    w_lip_vis(nn,:) = g_w_lip_vis/N_vis *(exp(K_lip_vis*(cos((prefs_vis-prefs_lip(nn))/180 *2*pi)-1)));  % << MT input >>
    w_lip_targ(nn,:) = g_w_lip_targ/N_vis *(exp(K_lip_targ*(cos((prefs_vis-prefs_lip(nn))/180 *2*pi)-1)));  % << Target input >>
    w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...   % << LIP recurrent >>
        ((exp(K_lip_lip*(cos((prefs_lip-prefs_lip(nn))/180*2*pi)-1)))-...
        amp_I*(exp(K_lip_I*(cos((prefs_lip-prefs_lip(nn))/180*2*pi)-1))))...
        + dc_w_lip_lip;
end

% << MT correlation matrix >>
cov_vis = zeros(N_vis,N_vis);
for nn=1:N_vis
    cov_vis(nn,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(nn))/180 *2*pi)-1));
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
    proba_vis = ((max_rate_vis-b_vis)*exp(K_vis*(cos((prefs_vis'-stim_this)/180*2*pi)-1))+b_vis)*dt;
    
    % proba_targ: proba of firing of target neurons in response to visual targets
    pos_targ = stim_this + [0:180/N_targets:179];
    proba_target = zeros(N_lip,1);
    for nn=1:N_targets
        proba_target = proba_target + ((max_rate_targ-b_targ)*...
            exp(K_targ*(cos((prefs_vis'-pos_targ(nn))/180*2*pi)-1))+b_targ)*dt;
    end
    proba_target = proba_target/(slope_norm_targ*(N_targets/2)+dc_norm_targ);   % << Divisive normalization >>
    
    
    % === Begin trials ===
    for mm = 1:2 % Train and test sets
        
        % --- Generate input spike trains ---
        for tt = 1:N_trial % For each trial
            
            if round(tt/50)==tt/50
                fprintf('%d  ',tt);
            end
            
            % << Generate Poisson spike train using definition. >>
            % Note that trial_dur_total now includes the pre-trial time.
            spikes_target{mm}(:,:,tt,ss) = rand(N_vis,trial_dur_total)<(proba_target*ones(1,trial_dur_total));
            
            %%% proba_out01: probability of firing of output neuron during
            %%% pretrial in response to visual target alone. Set to proba_in2.
            %%% proba_out02: same for testing set
            proba_out01__pretrial = proba_target*ones(1,stim_on_time)/dt; % << Why /delta_t here? Turn to firing rate? Yes. >>
            proba_out02__pretrial = proba_target*ones(1,stim_on_time)/dt;
            
            %%% spikes_out01: pretrial output spike train for training set
            %%% spikes_out11: pretrial output spike train for testing set
            %%% during pretrial, output neurons simply respond like visual units
            % << Because we actually want to initialize LIP response, not MT response >>
            spikes_out01__pretrial_train = spikes_target{mm}(:,1:stim_on_time);
            spikes_out02__pretrial_test = spikes_in22__target_test(:,1:stim_on_time);
            
            proba_out_av0__pretrial = proba_out_av0__pretrial + proba_out01__pretrial;  % << Averaged firing rate across N_trials >>
            proba_out50_pre(count,:) = proba_out01__pretrial(50,:);
            proba_out1_pre(count,:) = proba_out01__pretrial(1,:);
            
            %%%trial
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % proba_out{1}: probability of firing nu_neurons x nu_time_steps. Training
            % set
            % spikes_out{1}: spike train of LIP neurons [nu_neurons x nu_time_steps] Training
            % set. Computed from proba_out{1}.
            % proba_out{2}: probability of firing nu_neurons x nu_time_steps. Testing
            % set.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            aux1_proba1 = zeros(N_lip,trial_dur_total+2);
            aux1_proba1(:,1)= proba_out01__pretrial(:,stim_on_time);  % << Initial value of proba_out >>
            aux2_proba1 = zeros(N_lip,trial_dur_total+2);
            
            aux1_proba2 = zeros(N_lip,trial_dur_total+2);
            aux1_proba2(:,1)= proba_out02__pretrial(:,stim_on_time);
            aux2_proba2 = zeros(N_lip,trial_dur_total+2);
            
            proba_lip{1} = zeros(N_lip,trial_dur_total);
            proba_lip{2} = zeros(N_lip,trial_dur_total);
            
            proba_lip{1}(:,1)= proba_out01__pretrial(:,stim_on_time);  % << Initial value of proba_out. in Hz!! >>
            proba_lip{2}(:,1)= proba_out02__pretrial(:,stim_on_time);
            
            spikes_lip{1}(:,:)= zeros(N_lip,trial_dur_total);
            spikes_lip{2}(:,:)= zeros(N_lip,trial_dur_total);
            spikes_lip{1}(:,1)= spikes_out01__pretrial_train(:,stim_on_time);
            spikes_lip{2}(:,1)= spikes_out02__pretrial_test(:,stim_on_time);
            
            % spike_in: input spike trains due to motion. Two sets for
            % training and testing sets.
            
            % <<
            %    Add correlated noise to MT response.
            %    Here is the bug Alex mentioned in his email:
            %       "This code may have one bug: when I draw the randow
            %       seed for the firing rate of the MT neurons to introduce correlated
            %       noise, I used to do so only once at the beginning of each trial.
            %       Instead, the random seed should be drawn again at every time step on
            %       every trial, otherwise, you induce temporal correlations in the MT
            %       responses."
            %    However, it seems that not only the random seed fails to be updated, but even the noise is fixed! [randn(N_vis,1)]
            % >>
            %             aux_proba_in = proba_in__MTmotion + w_mt*randn(N_vis,1);
            
            %             spikes_in{m}(:,1:trial_dur)= rand(N_vis,trial_dur)<(aux_proba_in*ones(1,trial_dur));
            
            aux_proba_in = proba_vis*ones(1,trial_dur_total) + w_cov_vis*randn(N_vis,trial_dur_total);
            spikes_vis{m}(:,1:trial_dur_total)= rand(N_vis,trial_dur_total)<(aux_proba_in);
            
            %          aux_proba_in = proba_in_temp + w_mt*randn(N_vis,1);
            %          spikes_in{m}(:,101:200)= rand(N_vis,100)<(aux_proba_in*ones(1,100));
            
            spikes_target(:,1:trial_dur_total) = rand(N_vis,trial_dur_total)<(proba_target*ones(1,trial_dur_total));
            spikes_in22__target_test(:,1:trial_dur_total) = rand(N_vis,trial_dur_total)<(proba_target*ones(1,trial_dur_total));
            
            count_time50__this_trial=1;
            count_time100__this_trial=1;
            k=1;
            
            while k<=trial_dur_total-1
                
                % << == Network dynamics == >>
                % << I don't understand yet why there are two variables aux1 and aux2 >>
                aux1_proba1(:,k+1) = (1-dt/time_const_vis)*aux1_proba1(:,k)...   % << Self dynamics >>
                    +1/time_const_vis*(w_lip_vis*spikes_vis{1}(:,k)...     % << MT motion input >>
                    +w_lip_targ*att_gain2_targ*spikes_target(:,k));             % << Target visual input (att_gain2 has been set to 0) >>
                
                aux2_proba1(:,k+1) = (1-dt/time_const_lip)*aux2_proba1(:,k)+...
                    +1/time_const_lip*(w_lip_lip*spikes_lip{1}(:,k));  % << LIP recurrent >>
                
                proba_lip{1}(:,k+1) = aux1_proba1(:,k+1) + aux2_proba1(:,k+1) + bias_lip;  % << Bias term: b_out >>
                
                aux1_proba2(:,k+1) = (1-dt/time_const_vis)*aux1_proba2(:,k)...
                    +1/time_const_vis*(w_lip_vis*spikes_vis{2}(:,k)...
                    +w_lip_targ*att_gain2_targ*spikes_in22__target_test(:,k));
                aux2_proba2(:,k+1) = (1-dt/time_const_lip)*aux2_proba2(:,k)...
                    +1/time_const_lip*(w_lip_lip*spikes_lip{2}(:,k));
                proba_lip{2}(:,k+1) = aux1_proba2(:,k+1) + aux2_proba2(:,k+1) + bias_lip;
                
                % << Poisson with ReLU. Note that proba_out is in Hz!! >>
                spikes_lip{1}(:,k+1)  = (rand(N_lip,1) < dt*(proba_lip{1}(:,k+1)-threshold_lip));
                spikes_lip{2}(:,k+1) = (rand(N_lip,1) < dt*(proba_lip{2}(:,k+1)-threshold_lip));
                
                %% variable used for stopping the integration
                %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
                %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
                decision_ac(:,k+1) = proba_lip{1}(:,k+1);
                
                % collect spike times into a cell
                % << Only for the 50th and the 100th neuron >>
                if spikes_lip{1}(50,k+1) == 1
                    spike_time50{count}(count_time50__this_trial) = k+1;
                    count_time50__this_trial = count_time50__this_trial+1;
                end
                
                if spikes_lip{1}(100,k+1) == 1
                    spike_time100{count}(count_time100__this_trial) = k+1;
                    count_time100__this_trial = count_time100__this_trial+1;
                end
                
                %termination
                if ((max(decis_flag*decision_ac(:,k+1))> decis_thres) && (RT(count)==0)) || (k==trial_dur_total-1)
                    RT(count,1) = k;
                    if k>wind_size-1   % << Enough data for Fisher info of at least one window >>
                        last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,k-wind_size+1:k)')';
                    else
                        last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,1:k)')';
                    end
                    last_proba(:,count) = proba_lip{1}(:,k);
                    k=trial_dur_total;  % << Terminate trial immediately>>
                end
                k=k+1;
            end
            
            %%% collect spike and firing probabilities for Tin and Tout neurons
            
            proba_out_av1 = proba_out_av1 + proba_lip{1};
            proba_out50(count,:) = proba_lip{1}(50,:);
            proba_out75(count,:) = proba_lip{1}(75,:);
            proba_out100(count,:) = proba_lip{1}(100,:);
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
            
            proba_out_trial01(:,count) = sum(proba_lip{1}(:,1:trial_dur_total)')'/wind_size;
            proba_out_trial02(:,count) = sum(proba_lip{2}(:,1:trial_dur_total)')'/wind_size;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% collect counts in small time windows
            %% m=1 training data, m=2 testing data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for n=1:N_window
                proba_count_lip{m,n}(:,count)= ...
                    sum(proba_lip{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')'/wind_size;
                spikes_count_lip{m,n}(:,count)= ...
                    sum(spikes_lip{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')';
                spikes_count_vis{m,n}(:,count)= ...
                    sum(spikes_vis{m}(:,1:wind_offset+wind_size*n)')';
            end
            
            count = count+1;
        end % of trial
    end % of train/test
    
    fprintf('\n');
end


