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
N_trial = 100;

dt = 1e-3; % Size of time bin in seconds
pre_trial_dur = 50; % in time bins. must be 50 or more
trial_dur = 200; %% in time bins

if pre_trial_dur>trial_dur
    error('pre_trial_dur cannot be longer than trial_dur');
end

decis_flag = 0; %% decis_flag=1 means that the trial stops at the bound (reaction time version)
decis_thres = 58; %% bound height

% === Stimuli ===
stim_pos = [90 92]; 
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

spikes_in_trial11 = zeros(N_vis,N_trial*2);
spikes_in_trial12 = zeros(N_vis,N_trial*2);

spikes_in21 = zeros(N_vis,trial_dur);
spikes_in22 = zeros(N_vis,trial_dur);

aux1_proba1 = zeros(N_lip,trial_dur+2);
aux2_proba1 = zeros(N_lip,trial_dur+2);


decision_ac = zeros(N_lip,trial_dur+2);

aux1_proba2 = zeros(N_lip,trial_dur+2);
aux2_proba2 = zeros(N_lip,trial_dur+2);


spikes_in = cell(2,1);
proba_out = cell(2,1);
spikes_out = cell(2,1);
for m=1:2
   spikes_in{m} = zeros(N_lip,trial_dur);
   proba_out{m} = zeros(N_lip,trial_dur);
   spikes_out{m} = zeros(N_lip,trial_dur);
end

proba_count = cell(2,N_window);
spikes_count = cell(2,N_window);
for m=1:2
for n=1:N_window    
   proba_count{m,n} = zeros(N_lip,N_trial*2);
   spikes_count{m,n} = zeros(N_lip,N_trial*2);
   spikes_count_in{m,n} = zeros(N_vis,N_trial*2);
end
end

if shuffle_flag==1
   spikes_count_shuff = cell(2,N_window);
   for m=1:2
   for n=1:N_window    
      spikes_count_shuff{m,n} = zeros(N_lip,N_trial*2);
   end
   end
end

proba_out01  = zeros(N_lip,pre_trial_dur);
proba_out02  = zeros(N_lip,pre_trial_dur);
spikes_out01 = zeros(N_lip,pre_trial_dur);
spikes_out02 = zeros(N_lip,pre_trial_dur);

proba_out_av0__pretrial  = zeros(N_lip,pre_trial_dur);
proba_out_av1  = zeros(N_lip,trial_dur);

spikes_out_trial01 = zeros(N_lip,N_trial*2);
proba_out_trial01  = zeros(N_lip,N_trial*2);
spikes_out_trial02 = zeros(N_lip,N_trial*2);
proba_out_trial02  = zeros(N_lip,N_trial*2);

proba_out50_pre=zeros(N_trial*2,pre_trial_dur);
proba_out1_pre=zeros(N_trial*2,pre_trial_dur);
proba_out50=zeros(N_trial*2,trial_dur);
proba_out1=zeros(N_trial*2,trial_dur);

w_oi = zeros(N_lip,N_vis);
w_oi2 = zeros(N_lip,N_vis);
w_oo = zeros(N_lip,N_lip);

% << Network connections >>
for j=1:N_lip
    w_oi(j,:) = g_w_lip_vis/N_vis *(exp(K_lip_vis*(cos((prefs_vis-prefs_lip(j))/180 *2*pi)-1)));  % << MT input >>
    w_oi2(j,:) = g_w_lip_targ/N_vis *(exp(K_lip_targ*(cos((prefs_vis-prefs_lip(j))/180 *2*pi)-1)));  % << Target input >>
    w_oo(j,:) = g_w_lip_lip/N_lip*...   % << LIP recurrent >>
                ((exp(K_lip_lip*(cos((prefs_lip-prefs_lip(j))/180*2*pi)-1)))-...
                amp_I*(exp(K_lip_I*(cos((prefs_lip-prefs_lip(j))/180*2*pi)-1))))...
                + dc_w_lip_lip;
end

% << MT correlation matrix >>
for j=1:N_vis
    cov_mt(j,:) = var_vis*exp(K_cov_vis*(cos((prefs_vis-prefs_vis(j))/180 *2*pi)-1));
end
w_mt = real(sqrtm(cov_mt));

RT = zeros(N_trial*2,1);


spikes_out50 = zeros(N_trial*2,trial_dur);
spikes_out75 = zeros(N_trial*2,trial_dur);
spikes_out100 = zeros(N_trial*2,trial_dur);
last_spikes = zeros(N_lip,2*N_trial);
last_proba = zeros(N_lip,2*N_trial);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%
count =1;
%% there are two main loops, for two different positions of the input
for ss = 1:length(stim_pos)  % Motion directions 
    stim_this = stim_pos(ss);
    fprintf('%dth  pos: %4.2f\n  ',ss, stim_this);
    
    %    %%% proba_in: proba of firing of input neurons in response to motion
    max_rate_in = r_spont_vis + b_pref_vis*coherence;
    b_in = r_spont_vis + b_null_vis*coherence;
    proba_in__MTmotion = ((max_rate_in-b_in)*exp(K_vis*(cos((prefs_vis'-stim_this)/180*2*pi)-1))+b_in)*dt;
    
    max_rate_in = r_spont_vis + b_pref_vis*coherence*2;
    b_in = r_spont_vis + b_null_vis*coherence*2;
    proba_in_temp = ((max_rate_in-b_in)*exp(K_vis*(cos((prefs_vis'-stim_this)/180*2*pi)-1))+b_in)*dt;
    
    %      max_rate_in = b_pref*coherence;
    %      proba_in = (max_rate_in/2*exp(K_in*(cos((prefs_vis'-pos)/180*2*pi)-1)))*delta_t;
    
    %     proba_in_temp = (max_rate_in/2*exp(K_in*(cos((prefs_vis'-pos)/180*2*pi)-1)))*delta_t; ;
    %     proba_in_temp2 = (max_rate_in*exp(K_in*(cos((prefs_vis'-pos)/180*2*pi)-1)))*delta_t; ;
    
    %    max_rate_in = r_spont_MT + b_pref*coherence;
    %    b_in = r_spont_MT + b_null*coherence;
    %    proba_in_temp = ((max_rate_in-b_in)*exp(K_in*(cos((prefs_vis'-pos)/180*2*pi)-1))+b_in)*delta_t;
    
    %%% proba_in2: proba of firing of input neurons  in response to visual
    %%% targets
    
    % << Initialize firing rate of LIP neurons during "pretrial" by driving LIP with proba_in2. >>
    pos_targ = stim_this+[0:180/N_targets:179];
    proba_in2__target = zeros(N_lip,1);
    for m=1:N_targets
        proba_in2__target = proba_in2__target + ((max_rate_targ-b_targ)*...
            exp(K_targ*(cos((prefs_vis'-pos_targ(m))/180*2*pi)-1))+b_targ)*dt;
    end
    proba_in2__target = proba_in2__target/(slope_norm_targ*(N_targets/2)+dc_norm_targ);   % << Divisive normalization >>
    
    
    for j=1:N_trial
        if round(j/50)==j/50
            fprintf('%d  ',j);
        end
        
        %%% pretrial
        
        %% input spike trains during pretrial.
        %% 21 and 22: input from visual targets for training and testing set
        %% trial_dur is used here out of lazinness. It should be
        %% pre_trial_dur but we recycle the variable that will be used again
        %% during the trial.
        
        % << Generate Poisson spike train for pre_trial_dur using definition. >>
        spikes_in21__target_train(:,1:trial_dur) = rand(N_vis,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        spikes_in22__target_test(:,1:trial_dur) = rand(N_vis,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        
        %%% proba_out01: probability of firing of output neuron during
        %%% pretrial in response to visual target alone. Set to proba_in2.
        %%% proba_out02: same for testing set
        proba_out01__pretrial = proba_in2__target*ones(1,pre_trial_dur)/dt; % << Why /delta_t here? Turn to firing rate? Yes. >>
        proba_out02__pretrial = proba_in2__target*ones(1,pre_trial_dur)/dt;
        
        %%% spikes_out01: pretrial output spike train for training set
        %%% spikes_out11: pretrial output spike train for testing set
        %%% during pretrial, output neurons simply respond like visual units
        % << Because we actually want to initialize LIP response, not MT response >>
        spikes_out01__pretrial_train = spikes_in21__target_train(:,1:pre_trial_dur);
        spikes_out02__pretrial_test = spikes_in22__target_test(:,1:pre_trial_dur);
        
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
        aux1_proba1 = zeros(N_lip,trial_dur+2);
        aux1_proba1(:,1)= proba_out01__pretrial(:,pre_trial_dur);  % << Initial value of proba_out >>
        aux2_proba1 = zeros(N_lip,trial_dur+2);
        
        aux1_proba2 = zeros(N_lip,trial_dur+2);
        aux1_proba2(:,1)= proba_out02__pretrial(:,pre_trial_dur);
        aux2_proba2 = zeros(N_lip,trial_dur+2);
        
        proba_out{1} = zeros(N_lip,trial_dur);
        proba_out{2} = zeros(N_lip,trial_dur);
        
        proba_out{1}(:,1)= proba_out01__pretrial(:,pre_trial_dur);  % << Initial value of proba_out. in Hz!! >>
        proba_out{2}(:,1)= proba_out02__pretrial(:,pre_trial_dur);
        
        spikes_out{1}(:,:)= zeros(N_lip,trial_dur);
        spikes_out{2}(:,:)= zeros(N_lip,trial_dur);
        spikes_out{1}(:,1)= spikes_out01__pretrial_train(:,pre_trial_dur);
        spikes_out{2}(:,1)= spikes_out02__pretrial_test(:,pre_trial_dur);
        
        %% spike_in: input spike trains due to motion. Two sets for
        %% training and testing sets.
        for m=1:2
            
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
            
            aux_proba_in = proba_in__MTmotion*ones(1,trial_dur) + w_mt*randn(N_vis,trial_dur);            
            spikes_in{m}(:,1:trial_dur)= rand(N_vis,trial_dur)<(aux_proba_in);

            %          aux_proba_in = proba_in_temp + w_mt*randn(N_vis,1);
            %          spikes_in{m}(:,101:200)= rand(N_vis,100)<(aux_proba_in*ones(1,100));
        end
        
        spikes_in21__target_train(:,1:trial_dur) = rand(N_vis,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        spikes_in22__target_test(:,1:trial_dur) = rand(N_vis,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        
        count_time50__this_trial=1;
        count_time100__this_trial=1;
        k=1;
        
        while k<=trial_dur-1
            
            % << == Network dynamics == >>
            % << I don't understand yet why there are two variables aux1 and aux2 >>
            aux1_proba1(:,k+1) = (1-dt/time_const_vis)*aux1_proba1(:,k)...   % << Self dynamics >>
                +1/time_const_vis*(w_oi*spikes_in{1}(:,k)...     % << MT motion input >>
                +w_oi2*att_gain2_targ*spikes_in21__target_train(:,k));             % << Target visual input (att_gain2 has been set to 0) >>
            
            aux2_proba1(:,k+1) = (1-dt/time_const_lip)*aux2_proba1(:,k)+...
                +1/time_const_lip*(w_oo*spikes_out{1}(:,k));  % << LIP recurrent >>
            
            proba_out{1}(:,k+1) = aux1_proba1(:,k+1) + aux2_proba1(:,k+1) + bias_lip;  % << Bias term: b_out >>
            
            aux1_proba2(:,k+1) = (1-dt/time_const_vis)*aux1_proba2(:,k)...
                +1/time_const_vis*(w_oi*spikes_in{2}(:,k)...
                +w_oi2*att_gain2_targ*spikes_in22__target_test(:,k));
            aux2_proba2(:,k+1) = (1-dt/time_const_lip)*aux2_proba2(:,k)...
                +1/time_const_lip*(w_oo*spikes_out{2}(:,k));
            proba_out{2}(:,k+1) = aux1_proba2(:,k+1) + aux2_proba2(:,k+1) + bias_lip;
            
            % << Poisson with ReLU. Note that proba_out is in Hz!! >>
            spikes_out{1}(:,k+1)  = (rand(N_lip,1) < dt*(proba_out{1}(:,k+1)-threshold_lip));
            spikes_out{2}(:,k+1) = (rand(N_lip,1) < dt*(proba_out{2}(:,k+1)-threshold_lip));
            
            %% variable used for stopping the integration
            %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
            %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
            decision_ac(:,k+1) = proba_out{1}(:,k+1);
            
            % collect spike times into a cell
            % << Only for the 50th and the 100th neuron >>
            if spikes_out{1}(50,k+1) == 1
                spike_time50{count}(count_time50__this_trial) = k+1;
                count_time50__this_trial = count_time50__this_trial+1;
            end
            
            if spikes_out{1}(100,k+1) == 1
                spike_time100{count}(count_time100__this_trial) = k+1;
                count_time100__this_trial = count_time100__this_trial+1;
            end
            
            %termination
            if ((max(decis_flag*decision_ac(:,k+1))> decis_thres) && (RT(count)==0)) || (k==trial_dur-1)
                RT(count,1) = k;
                if k>wind_size-1   % << Enough data for Fisher info of at least one window >>
                    last_spikes__last_Fisher_win(:,count) = sum(spikes_out{1}(:,k-wind_size+1:k)')';  
                else
                    last_spikes__last_Fisher_win(:,count) = sum(spikes_out{1}(:,1:k)')'; 
                end
                last_proba(:,count) = proba_out{1}(:,k);
                k=trial_dur;  % << Terminate trial immediately>>
            end
            k=k+1;
        end
        
        %%% collect spike and firing probabilities for Tin and Tout neurons
        
        proba_out_av1 = proba_out_av1 + proba_out{1};
        proba_out50(count,:) = proba_out{1}(50,:);
        proba_out75(count,:) = proba_out{1}(75,:);
        proba_out100(count,:) = proba_out{1}(100,:);
        spikes_out50(count,:) = spikes_out{1}(50,:);
        spikes_out75(count,:) = spikes_out{1}(75,:);
        spikes_out100(count,:) = spikes_out{1}(100,:);
        spikes_in50(count,:) = spikes_in{1}(50,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% training and testing data. Spike counts and average probabilities.
        %% size nu_unit x 2 N_trial (N_trial for each of the 2 positions)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spikes_in_trial11(:,count) = sum(spikes_in{1}')';
        spikes_in_trial12(:,count) = sum(spikes_in{2}')';
        
        spikes_out_trial01(:,count) = sum(spikes_out{1}')';
        spikes_out_trial02(:,count) = sum(spikes_out{2}')';
        
        proba_out_trial01(:,count) = sum(proba_out{1}(:,1:trial_dur)')'/wind_size;
        proba_out_trial02(:,count) = sum(proba_out{2}(:,1:trial_dur)')'/wind_size;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% collect counts in small time windows
        %% m=1 training data, m=2 testing data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for m=1:2
            for n=1:N_window
                proba_count{m,n}(:,count)= ...
                    sum(proba_out{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')'/wind_size;
                spikes_count{m,n}(:,count)= ...
                    sum(spikes_out{m}(:,wind_offset+wind_size*(n-1)+1:wind_size*n)')';
                spikes_count_in{m,n}(:,count)= ...
                    sum(spikes_in{m}(:,1:wind_offset+wind_size*n)')';
            end
        end
        
        count = count+1;
    end
    fprintf('\n');
end


