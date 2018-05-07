%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear
clear global

load h

global Y_targ X_train;
rand('state',sum(100*clock))

% %{

N_in = 100;
pos_in=[0:180/N_in:179.9999];

N_out = 100;
step_out = 180/N_out;
pos_out=[0:step_out:179.9999];
nu_targets = 2;

decis_thres = 58; %% bound height
decis_flag = 0; %% decis_flag=1 means that the trial stops at the bound

peak = 90;
pre_trial_dur= 50; %% must be 50 or more
trial_dur = 200; %%all times are in ms

if pre_trial_dur>trial_dur
    error('pre_trial_dur cannot be longer than trial_dur');
end

delta_t = 1e-3; %size of time bin in seconds
step_test =2; % difference in s between the two test stimuli

nu_trial = 100;
coherence = 6.4;

% Parameters for MT from Mazurek and Shadlen, except width.
% K_cov_mt and var_mt controls the correlations.
%  << cov_mt(j,:) = var_mt*exp(K_cov_mt*(cos((pos_in-pos_in(j))/180 *2*pi)-1));
%     w_mt = real(sqrtm(cov_mt));
%     aux_proba_in = proba_in + w_mt*randn(N_in,1); >>
r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;
K_cov_mt = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% var_mt = 1e-5; % Original. Use bugged code, cor90 = 0.06
% var_mt = 7e-5; % Bugged code, cor90 = 0.2; Corrected code, cor90 = 0.0027
var_mt = 0.15; % Predicted by my estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter of the activity evoked by the visual targets
max_rate_in2 = 42; %% in Hz
b_in2 = 0;
K_in2 = 4;

% Weight parameters
b_out = 0;
K_oin = 5;
K_oin2= 5;
K_oo = 10;
K_oI = 2;
threshold = 0.0;

%time constant pre trial
pre_time_const_out = 1000e-3;
pre_time_const_in  = pre_time_const_out;

%time constant for integration
time_const_out = 1000e-3;
time_const_in  = time_const_out;

% MT weights
g_w_oi =65;

% Visual evoked weights
g_w_oi2= 40;
%%% drop in attention to visual target once motion stimulus appears.
att_gain2 = 0;

% recurrent connectivity in LIP
g_w_oo = 21;
amp_i = 0.0;
dc_w_oo = -0.17;

ext_noise=0;

shuffle_flag = 0; %shuffles trials and computes shuffled Fisher information
info_out_flag = 1; %controls whether Fisher information in the output layer is being computed
info_in_flag = 1; %controls whether Fisher information in the input layer is being computed


%%% parameters of the divisive normalization for visually evoked activity
slope_norm=0.5;
dc_norm = 0.5;

%%% duration of the window used to compute Fisher information
wind_size = 50;
wind_offset = 0;
nu_window = 4;

%}

%{
% main parameters for the paper

N_in = 100;
pos_in=[0:180/N_in:179.9999];
N_out = 100;
step_out = 180/N_out;
pos_out=[0:step_out:179.9999];
nu_targets = 2;
decis_thres = 58;
decis_flag = 1;

peak = 90;
pre_trial_dur= 50; %% must be 50 or more
trial_dur = 1000; %%all times are in ms

if pre_trial_dur>trial_dur
   error('pre_trial_dur cannot be longer than trial_dur');
end

delta_t = 1e-3; %size of time bin in second
% difference in s between the two test stimuli
step_test =3;

nu_trial = 400;
coherence = 51.2;

% Parameters for MT from Mazurek and Shadlen, except width. K_cov_mt
% and var_mt controls the correlations.
r_spont_MT = 20;
b_pref = 0.4;
b_null = -0.2;
K_in = 4;
K_cov_mt = 2;
var_mt = 1e-5;

% parameter of the visually evoked activity
max_rate_in2 = 42;
b_in2 = 0;
K_in2 = 4;
% no normalisation depending on target number (as 2 targets anyway)
% slope_norm*(nu_targets/2)+dc_norm = 1 for nu_targets = 2
slope_norm = 0.5;
dc_norm = 0.5;

% Weight parameters
b_out = 0;
K_oin = 5;
K_oin2= 5;
K_oo = 10;
K_oI = 2;
threshold = 0.0;

%time constant pre trial
pre_time_const_out = 1000e-3;
pre_time_const_in  = pre_time_const_out;

%time constant for integration
time_const_out = 1000e-3;
time_const_in  = time_const_out;
 
% MT weight
g_w_oi = 65;

% Visual evoked weights
g_w_oi2= 40;
%%% drop in attention to visual target once motion stimulus appears.
att_gain2 = 0;


%g_w_oo = 20;%%% parameter for Cosyne
g_w_oo = 21;6
amp_i = 0.0;
dc_w_oo = -0.17;

ext_noise=0;

shuffle_flag = 0;
info_out_flag = 1;
info_in_flag = 1;

%%% duration of the window used to compute Fisher
wind_size = 50;
wind_offset = 0;
nu_window = 4;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialization



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%
count =1;
%% there are two main loops, for two different positions of the input
for mm = 1:2  % Motion directions 
    pos = peak+step_test*(mm-step_test/2);
    fprintf('%d  pos: %4.2f\n  ',mm, pos);
    
    %    %%% proba_in: proba of firing of input neurons in response to motion
    max_rate_in = r_spont_MT + b_pref*coherence;
    b_in = r_spont_MT + b_null*coherence;
    proba_in_MTmotion = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t;
    
    max_rate_in = r_spont_MT + b_pref*coherence*2;
    b_in = r_spont_MT + b_null*coherence*2;
    proba_in_temp = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t;
    
    %      max_rate_in = b_pref*coherence;
    %      proba_in = (max_rate_in/2*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t;
    
    %     proba_in_temp = (max_rate_in/2*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t; ;
    %     proba_in_temp2 = (max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1)))*delta_t; ;
    
    %    max_rate_in = r_spont_MT + b_pref*coherence;
    %    b_in = r_spont_MT + b_null*coherence;
    %    proba_in_temp = ((max_rate_in-b_in)*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))+b_in)*delta_t;
    
    %%% proba_in2: proba of firing of input neurons  in response to visual
    %%% targets
    
    % << Initialize firing rate of LIP neurons during "pretrial" by driving LIP with proba_in2. >>
    pos_targ = pos+[0:180/nu_targets:179];
    proba_in2__target = zeros(N_out,1);
    for m=1:nu_targets
        proba_in2__target = proba_in2__target + ((max_rate_in2-b_in2)*...
            exp(K_in2*(cos((pos_in'-pos_targ(m))/180*2*pi)-1))+b_in2)*delta_t;
    end
    proba_in2__target = proba_in2__target/(slope_norm*(nu_targets/2)+dc_norm);   % << Divisive normalization >>
    
    
    for j=1:nu_trial
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
        spikes_in21__target_train(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        spikes_in22__target_test(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        
        %%% proba_out01: probability of firing of output neuron during
        %%% pretrial in response to visual target alone. Set to proba_in2.
        %%% proba_out02: same for testing set
        proba_out01__pretrial = proba_in2__target*ones(1,pre_trial_dur)/delta_t; % << Why /delta_t here? Turn to firing rate? Yes. >>
        proba_out02__pretrial = proba_in2__target*ones(1,pre_trial_dur)/delta_t;
        
        %%% spikes_out01: pretrial output spike train for training set
        %%% spikes_out11: pretrial output spike train for testing set
        %%% during pretrial, output neurons simply respond like visual units
        % << Because we actually want to initialize LIP response, not MT response >>
        spikes_out01__pretrial_train = spikes_in21__target_train(:,1:pre_trial_dur);
        spikes_out02__pretrial_test = spikes_in22__target_test(:,1:pre_trial_dur);
        
        proba_out_av0__pretrial = proba_out_av0__pretrial + proba_out01__pretrial;  % << Averaged firing rate across nu_trials >>
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
        aux1_proba1 = zeros(N_out,trial_dur+2);
        aux1_proba1(:,1)= proba_out01__pretrial(:,pre_trial_dur);  % << Initial value of proba_out >>
        aux2_proba1 = zeros(N_out,trial_dur+2);
        
        aux1_proba2 = zeros(N_out,trial_dur+2);
        aux1_proba2(:,1)= proba_out02__pretrial(:,pre_trial_dur);
        aux2_proba2 = zeros(N_out,trial_dur+2);
        
        proba_out{1} = zeros(N_out,trial_dur);
        proba_out{2} = zeros(N_out,trial_dur);
        
        proba_out{1}(:,1)= proba_out01__pretrial(:,pre_trial_dur);  % << Initial value of proba_out. in Hz!! >>
        proba_out{2}(:,1)= proba_out02__pretrial(:,pre_trial_dur);
        
        spikes_out{1}(:,:)= zeros(N_out,trial_dur);
        spikes_out{2}(:,:)= zeros(N_out,trial_dur);
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
            %    However, it seems that not only the random seed fails to be updated, but even the noise is fixed! [randn(N_in,1)]
            % >>
            
%             aux_proba_in = proba_in_MTmotion + w_mt*randn(N_in,1);
%             spikes_in{m}(:,1:trial_dur)= rand(N_in,trial_dur)<(aux_proba_in*ones(1,trial_dur));  %#ok<*SAGROW>
             
            aux_proba_in = proba_in_MTmotion*ones(1,trial_dur) + w_mt*randn(N_in,trial_dur);            
            spikes_in{m}(:,1:trial_dur)= rand(N_in,trial_dur)<(aux_proba_in);

            %          aux_proba_in = proba_in_temp + w_mt*randn(N_in,1);
            %          spikes_in{m}(:,101:200)= rand(N_in,100)<(aux_proba_in*ones(1,100));
        end
        
        spikes_in21__target_train(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        spikes_in22__target_test(:,1:trial_dur) = rand(N_in,trial_dur)<(proba_in2__target*ones(1,trial_dur));
        
        count_time50__this_trial=1;
        count_time100__this_trial=1;
        k=1;
        
        while k<=trial_dur-1
            
            % << == Network dynamics == >>
            % << I don't understand yet why there are two variables aux1 and aux2 >>
            aux1_proba1(:,k+1) = (1-delta_t/time_const_in)*aux1_proba1(:,k)...   % << Self dynamics >>
                +1/time_const_in*(w_oi*spikes_in{1}(:,k)...     % << MT motion input >>
                +w_oi2*att_gain2*spikes_in21__target_train(:,k));             % << Target visual input (att_gain2 has been set to 0) >>
            
            aux2_proba1(:,k+1) = (1-delta_t/time_const_out)*aux2_proba1(:,k)+...
                +1/time_const_out*(w_oo*spikes_out{1}(:,k));  % << LIP recurrent >>
            
            proba_out{1}(:,k+1) = aux1_proba1(:,k+1) + aux2_proba1(:,k+1) + b_out;  % << Bias term: b_out >>
            
            aux1_proba2(:,k+1) = (1-delta_t/time_const_in)*aux1_proba2(:,k)...
                +1/time_const_in*(w_oi*spikes_in{2}(:,k)...
                +w_oi2*att_gain2*spikes_in22__target_test(:,k));
            aux2_proba2(:,k+1) = (1-delta_t/time_const_out)*aux2_proba2(:,k)...
                +1/time_const_out*(w_oo*spikes_out{2}(:,k));
            proba_out{2}(:,k+1) = aux1_proba2(:,k+1) + aux2_proba2(:,k+1) + b_out;
            
            % << Poisson with ReLU. Note that proba_out is in Hz!! >>
            spikes_out{1}(:,k+1)  = (rand(N_out,1) < delta_t*(proba_out{1}(:,k+1)-threshold));
            spikes_out{2}(:,k+1) = (rand(N_out,1) < delta_t*(proba_out{2}(:,k+1)-threshold));
            
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
        %% size nu_unit x 2 nu_trial (nu_trial for each of the 2 positions)
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
            for n=1:nu_window
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






%%%%%%%%%%%%%%%%%%%%%%%%
%% stats
%%%%%%%%%%%%%%%%%%%%%%%%
mean_in  = mean(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
mean_in2 = mean(spikes_in_trial11(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);
var_in = var(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);

mean_out = mean(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);
mean_out2 = mean(spikes_out_trial01(:,nu_trial+1:2*nu_trial)')/(trial_dur*delta_t);

mean_proba  = mean(proba_out_trial01(:,1:nu_trial)');
mean_proba2 = mean(proba_out_trial01(:,1:nu_trial)');

for m=1:2
    for n=1:nu_window
        mean_proba_count{m,n}=  mean(proba_count{m,n}(:,1:nu_trial)');
    end
end

var_out = var(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);
cov_in =  cov(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
corr_in =  cov_in./(diag(cov_in).^0.5*(diag(cov_in)').^0.5);

%compute the average correlation coefficients in MT
aux_mask = ones(N_in,1)*[1:N_in];
mask = ((1-cos(abs(aux_mask-aux_mask')/N_in*2*pi))/2)<.5;
mask = mask.*(1-eye(N_in));
cor90 = sum(sum(corr_in.*mask))/sum(sum(mask));

%%%%%%%%%%%%%%%%%%%%%%%% Debug the correlation problem %%%%%%%%%%%%%%%%%%%%%%
% HH20180503
fano = var_in ./ mean_in;

figure(1817); clf
subplot(1,3,1);
plot(mean_in); hold on; 
plot(proba_in_MTmotion/delta_t,'k--'); 
title('Firing rate'); ylim([0 max(ylim)*1.1])

subplot(1,3,2); 
plot(fano); hold on; 
plot(xlim,[1 1],'k--'); 
title('Fano'); ylim([0 max(ylim)*1.1])

subplot(1,3,3);
imagesc(corr_in); colorbar;
title(sprintf('cor90 = %g',cor90));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure();

mask = ((1-cos(abs(aux_mask-aux_mask')/N_in*2*pi))/2)>.5;
cor180 = sum(sum(corr_in.*mask))/sum(sum(mask));

cov_out =  cov(spikes_out_trial01(:,1:nu_trial)')/(trial_dur*delta_t);

difference = mean(spikes_out_trial01(:,1:nu_trial)')...
    -mean(spikes_out_trial01(:,nu_trial+1:nu_trial*2)');

proba_out_av0__pretrial = proba_out_av0__pretrial/(2*nu_trial);  % << Averaged firing rate. Divided by 2 because l = 1:2 >>
proba_out_av1 = proba_out_av1/(2*nu_trial);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if info_in_flag==1
    %input info
    [info_in info_train_in w_ole_in est_angle_in perc_right_test_in] = ...
        compute_info(spikes_in_trial11, spikes_in_trial12, step_test);
    
    %normalize for duration of trial
    info_in = info_in/(trial_dur*delta_t);
    info_train_in = info_train_in/(trial_dur*delta_t);
end

if info_out_flag==1
    % output info spikes
    [info_out info_train_out w_ole_out est_angle perc_right_test] = ...
        compute_info(spikes_out_trial01, spikes_out_trial02, step_test);
    
    %normalize for duration of trial
    info_out= info_out/(trial_dur*delta_t);
    info_train_out= info_train_out/(trial_dur*delta_t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info_outx=zeros(1,nu_window);
w_ole_outx = zeros(N_out+1,nu_window);
est_anglex= zeros(nu_window,nu_trial*2);

%%% information in MT
if info_in_flag==1
    for m=1:nu_window
        [info_inx(m) info_train_inx(m) w_ole_inx(:,m) est_angle_inx(m,:) perc_right_test_inx(m)] = ...
            compute_info(spikes_count_in{1,m}, spikes_count_in{2,m}, step_test);
        %normalize for duration of window
        info_inx(m) = info_inx(m) /(wind_size*delta_t);
        info_train_inx(m) = info_train_inx(m)/(wind_size*delta_t);
    end
end

%%%information in LIP
if info_out_flag==1
    for m=1:nu_window
        [info_outx(m) info_train_outx(m) w_ole_outx(:,m) est_anglex(m,:) perc_right_testx(m)] = ...
            compute_info(spikes_count{1,m}, spikes_count{2,m}, step_test);
        %normalize for duration of window
        info_outx(m) = info_outx(m) /(wind_size*delta_t);
        info_train_outx(m) = info_train_outx(m)/(wind_size*delta_t);
        
        [info_proba(m) info_train_proba(m) w_ole_proba(:,m) est_angle_proba(m,:) perc_right_test_proba(m)] = ...
            compute_info(proba_count{1,m}, proba_count{2,m}, step_test);
        %normalize for duration of window
        info_proba(m) = info_proba(m) /(wind_size*delta_t);
        info_train_proba(m) = info_train_proba(m)/(wind_size*delta_t);
    end
    
    
    
    if shuffle_flag ==1
        for j=1:N_out
            for m=1:2
                for n=1:nu_window
                    spikes_count_shuff{m,n}(j,1:nu_trial) = ...
                        spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+1);
                    spikes_count_shuff{m,n}(j,nu_trial+1:nu_trial*2) = ...
                        spikes_count{m,n}(j,mod(j-1:nu_trial-1+(j-1),nu_trial)+nu_trial+1);
                end
            end
        end
        
        for m=1:nu_window
            [info_shuff_out(m) info_train_shuff_out(m)...
                w_ole_shuff_out est_shuff_angle perc_right_shuff_test(m)] = ...
                compute_info(spikes_count_shuff{1,m}, spikes_count_shuff{2,m}, step_test);
            %normalize for duration of window
            info_shuff_out(m) = info_shuff_out(m)/(wind_size*delta_t);
            info_train_shuff_out(m) = info_train_shuff_out(m)/(wind_size*delta_t);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute info in small time windows using beck's grad descent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if info_out_flag==1
%    for m=1:nu_window
%       [binfo_train_out(m) binfo_out(m) bw_ole_out(:,m) bperc_right_train(m) bperc_right_test(m)] = ...
%          infobeck(spikes_count{1,m}, spikes_count{2,m}, step_test, 0.01);
%        %normalize for duration of window
%       binfo_out(m) = binfo_out(m) /(wind_size*delta_t);
%       binfo_train_out(m) = binfo_train_out(m)/(wind_size*delta_t);
%
%       [binfo_train_proba(m) binfo_proba(m) bw_ole_proba(:,m) bperc_right_train_proba(m) bperc_right_test_proba(m)] = ...
%          infobeck(proba_count{1,m}, proba_count{2,m}, step_test, 0.01);
%        %normalize for duration of window
%       binfo_proba(m) = binfo_proba(m) /(wind_size*delta_t);
%       binfo_train_proba(m) = binfo_train_proba(m)/(wind_size*delta_t);
%    end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate activation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out= mean_out;
in = mean_proba;

plot(mean_proba,2*log(1+exp(mean_proba/1000)),'-g');

parameters = fminsearch(@(param) dist2(param,in,out),[1000 2]);

alpha= parameters(1);
amp = parameters(2);

hold on
plot(mean_proba,mean_out,'o');
plot(mean_proba,amp*log(1+exp(mean_proba/alpha)),'-');
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute true info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if info_out_flag ==1
    ac_in = max_rate_in*exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1));
    ac_in_der = -max_rate_in * K_in * sin((pos_in'-pos)/180*2*pi).*...
        exp(K_in*(cos((pos_in'-pos)/180*2*pi)-1))/180*2*pi;
    
    info_in_true = sum(ac_in_der.^2./ac_in);
    
    true_cov_in =diag(ac_in);
    
    local_cov_out = diag(mean_out2);
    
    % 1e-6 is added to ensure invertibility in the next line when aux_info is
    % zero
    aux_info = w_oo*local_cov_out*w_oo' + diag(diag(ones(N_out)*1e-6));
    info_out_rate_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);
    
    D_mat = diag(amp/alpha*exp(mean_proba/alpha)./(1+exp(mean_proba/alpha)));
    D_mat = D_mat/max(max(D_mat));
    D_mat_inv = inv(D_mat);
    
    aux_info = D_mat_inv*local_cov_out*D_mat_inv + diag(diag(ones(N_out)*1e-6));
    
    %aux_info = local_cov_out + diag(diag(ones(N_out)*1e-6));
    info_out_true = ac_in_der'*w_oi'*inv(w_oi*true_cov_in*w_oi'+ aux_info)*(w_oi*ac_in_der);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate width of population activity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Estimate the width of the population activity \n');
data=mean_in;
[parameters val] = fminsearch(@(param) dist(param, data),[max(mean_in); K_in; peak; b_in]);

amp_fit= parameters(1);
K_fit = parameters(2);
peak_fit = parameters(3);
offset = parameters(4)
fit_in = amp_fit * exp(K_fit*(cos((pos_in'-peak_fit)/180*2*pi)-1)) + offset;

hwhh_in = acos(1/K_fit*log((amp_fit-exp(-2*K_fit))/(2*amp_fit))+1)/(2*pi)*180;


data=mean_out;
fit_out_init = max(mean_out)*exp(K_in*(cos((pos_in'-peak)/180*2*pi)-1));
[parameters val] = fminsearch(@(param) dist(param, data),[max(mean_out); K_in; peak; b_in]);

amp_fit= parameters(1);
K_fit  = parameters(2);
peak_fit = parameters(3);
offset = parameters(4)

fit_out = amp_fit*exp(K_fit*(cos((pos_in'-peak_fit)/180*2*pi)-1))+offset;

hwhh_out = acos(1/K_fit*log((amp_fit-exp(-2*K_fit))/(2*amp_fit))+1)/(2*pi)*180;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate PDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Estimate the pdfs \n');

for m=1:nu_window
    pdf{m} = exp(h_mat*spikes_count{1,m});
    pdf{m} = pdf{m}./(ones(N_out,1)*sum(pdf{m}));
    pdf_aver(m,:) = sum(pdf{m}')/sum(sum(pdf{m}'));
end

if shuffle_flag==1
    for m=1:nu_window
        pdf_shuff{m} = exp(h_mat*spikes_count_shuff{1,m});
        pdf_shuff{m} = pdf_shuff{m}./(ones(N_out,1)*sum(pdf_shuff{m}));
        pdf_shuff_aver{m} = sum(pdf_shuff{m}')/sum(sum(pdf_shuff{m}'));
    end
end

%%% compute distribution over ML estimates
ML_info(1) = 0;

for m=1:nu_window
    [aux_ML ML{m}] = max(pdf{m});
    pdf_ML{m} = hist(ML{m},[1:100]);
    [parameters val] =...
        fminsearch(@(param) dist(param, pdf_ML{m}),[max(pdf_ML{m}); K_in; peak; ,min(pdf_ML{m})]);
    K_fit = parameters(2);
    vardeg = 1/K_fit*(180/(2*pi))^2;
    ML_info(m+1) = 1/vardeg;
    ML_info(m+1) =ML_info(m+1)/(wind_size*delta_t);
end

%%% compute log odds over time
log_odds(1)=0;
for m=1:nu_window
    log_odds(m+1) = mean(log10(pdf{m}(50,:)./pdf{m}(100,:)));
end


%%% Compute information over time by fitting von mises functions to the pdf
%%% and using the K parameter to estimate information
pdf_info(1) = 0;
K_fit(1) = 10^-30;

for m=1:nu_window
    pdf_mean(m+1) = sum(pdf_aver(m,:).*pos_out);
    
    [parameters val] =...
        fminsearch(@(param) dist(param, pdf_aver(m,:)),[max(pdf_aver(m,:));...
        K_in; peak; min((pdf_aver(m,:)))]);
    K_fit(m+1) = parameters(2);
    vardeg = 1/K_fit(m+1)*(180/(2*pi))^2;
    pdf_est(m,:)= parameters(1)*exp(K_fit(m+1)*...
        (cos((pos_out-parameters(3))/180*2*pi)-1))+parameters(4);
    pdf_info(m+1) = 1/vardeg;
    pdf_info(m+1) = pdf_info(m+1)/(wind_size*delta_t);
end



%%% compute pdf at decision time based on spike counts over the last 100ms (the last Fisher window)
pdf_decis = exp(h_mat*last_spikes__last_Fisher_win);
pdf_decis = pdf_decis./(ones(N_out,1)*sum(pdf_decis));


%%% compute performance on each trial
%[junk pos_max] = max(pdf_decis);
[junk pos_max] = max(last_proba);
pos_max = pos_max/N_out*180;
bound_max = peak+step_test*(-0.5)+ 180/nu_targets/2;
bound_min = peak+step_test*(-0.5)- 180/nu_targets/2;
perf_trial = (bound_min<pos_max).*(pos_max<bound_max).*(RT'>0);


conf_level = log10(pdf_decis(50,:)./pdf_decis(100,:));
conf_level_proba = last_proba(50,:)-last_proba(100,:);

%%% compute overall performance
perf = sum(perf_trial)/(sum(RT>0));

RT_cor= RT.*(RT>0).*(perf_trial'>0);
RT_incor = RT.*(RT>0).*(perf_trial'<1);

mean_RT = sum(RT.*(RT>0).*(perf_trial'>0))/(sum((RT>0).*(perf_trial'>0)));

mean_RT_incor = sum(RT.*(RT>0).*(perf_trial'<1))/(sum((RT>0).*(perf_trial'<1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fano factor over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:nu_window
    size_fano = sum(mean(spikes_count{1,j}')~=0);
    aux_mean_fano = mean(spikes_count{1,j}');
    aux_mean_fano = aux_mean_fano+(aux_mean_fano==0);
    aux_fano = var(spikes_count{1,j}')./aux_mean_fano;
    
    fano(j) = sum(aux_fano)/size_fano;
    fano50(j) = aux_fano(1,50);
    fano100(j) = aux_fano(1,100);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% plots
%%%%%%%%%%%%%%%%%%%%%%%%
if info_in_flag==1
    fprintf('\n');
    fprintf('Info In test           : %2.3f\n',info_in);
    fprintf('Info In train          : %2.3f\n',info_train_in);
end

if info_out_flag==1
    fprintf('True Info rate Out test : %2.3f\n',info_out_rate_true);
    fprintf('True Info Out test      : %2.3f\n',info_out_true);
    fprintf('\n');
    
    fprintf('Percentage classification correct: %2.2f\n',perc_right_test_in);
    fprintf('\n');
    
    fprintf('Info Out test          : %2.3f\n',info_out);
    fprintf('Info Out train         : %2.3f\n',info_train_out);
    fprintf('Percentage classification correct: %2.2f\n\n',perc_right_test);
    fprintf('\n');
end

fprintf('Max firing rate: %2.3f\n',max(mean_out));
fprintf('\n');

fprintf('Norm of W_oo : %2.6f\n',w_oo(:,1)'*w_oo(:,1));
fprintf('Norm of W_oi: %2.6f\n',w_oi(:,1)'*w_oi(:,1));
fprintf('\n');
fprintf('\n');

fprintf('Width in      : %2.3f\n',hwhh_in);
fprintf('Width out      : %2.3f\n',hwhh_out);
fprintf('\n');

fprintf('Performance: %f\n', perf);
fprintf('Mean RT correct  : %f ms\n', mean_RT);
fprintf('Mean RT incorrect: %f ms\n', mean_RT_incor);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Move contents of plot_results here for a better experience of debugging. HH20161210
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_results;

subplot(221)
wind_times = wind_offset+[0:nu_window]*wind_size;

if nu_targets==2
    lik1 = [1 pdf_aver(1,50) pdf_aver(2,50) pdf_aver(3,50) pdf_aver(4,50)];
    lik2 = [1 pdf_aver(1,100) pdf_aver(2,100) pdf_aver(3,100) pdf_aver(4,100)];
    loglik= log10(lik1./lik2);
    %    if shuffle_flag ==1
    %       lik1_shuff = [1 pdf_shuff_aver{1}(50) pdf_shuff_aver{2}(50) pdf_shuff_aver{3}(50) pdf_shuff_aver{4}(50)];
    %       lik2_shuff = [1 pdf11_shuff_aver{1}(100) pdf_shuff_aver{2}(100) pdf_shuff_aver{4}(100) pdf_shuff_aver{4}(100)];
    %       loglik_shuff= log10(lik1_shuff./lik2_shuff);
    %    end
    
    % size   plot(wind_times, loglik,'*-');
    plot(wind_times, loglik,'*-');
    %    if shuffle_flag==1
    %        hold on
    %        plot([1 100 200 300 trial_dur], loglik_shuff,'r*-');
    %        hold off
    %    end
    title('Log odds');
elseif info_out_flag==1
    plot(wind_times,[0 info_outx],'r-o')
    hold on
    if shuffle_flag==1
        plot(wind_times, [0 info_shuff_out],'b-o')
    end
    plot(wind_times,[0 info_train_outx],'r:o');
    plot(wind_times,[0 (info_outx+info_train_outx)/2],'g:o');
    plot(wind_times,ML_info,'c:*');
    hold off
    title('Fisher information');
end


%%% plot information over time in input layer
%    plot(wind_times,[0 info_inx],'r-o')
%    hold on
%    plot(wind_times,[0 info_train_inx],'r:o');
%    plot(wind_times,[0 (info_inx+info_train_inx)/2],'g:o');
%    hold off
%    title('Fisher information in input layer');
%



subplot(222)
plot(pdf{1}(:,1),'.r');
hold on
plot(pdf{2}(:,1),'.b');
plot(pdf{3}(:,1),'.g');
plot(pdf{4}(:,1),'.c');
plot(pdf_aver(1,:),'r')
plot(pdf_aver(2,:),'b')
plot(pdf_aver(3,:),'g')
plot(pdf_aver(4,:),'c')
hold off

axis tight
ylabel('P(s|r)')
xlabel('Stimulus')





subplot(223)



for m=1:nu_window
    
    % mean_in  = mean(spikes_in_trial11(:,1:nu_trial)')/(trial_dur*delta_t);
    
    plot(pos_in,mean_in,'o')
    hold on
    plot(pos_out,mean_out,'or')
    plot(pos_in,fit_in,'-b')
    plot(pos_in,fit_out,'-r')
    hold off
    
    max_plot=max(max(mean_out),max(mean_in));
    axis([0 180 0 max_plot*1.2]);
    
    aux11{m} = mean(proba_count{1,m}');
    aux11{m} = aux11{m}.*(aux11{m}>0);
    %     aux11{m} = aux11{m}-min(aux11{m});
    %     aux11{m} = aux11{m}/max(aux11{m});
end

plot(aux11{1},'or')
hold on
plot(aux11{2},'ob')
plot(aux11{3},'og')
plot(aux11{4},'oc')
hold off




subplot(224)
%surf([0:180/(N_out-1):180], [0:180/(N_out-1):180],cov_out-diag(diag(cov_out)))
%shading interp
aux_proba= [proba_out_av0__pretrial(:,pre_trial_dur-49:10:pre_trial_dur) proba_out_av1(:,1:10:trial_dur)];
%surfl(aux_proba)

plot([-49:10:trial_dur],aux_proba(50,:))
hold on
plot([-49:10:trial_dur],aux_proba(100,:),'g')
plot([-49:10:trial_dur],aux_proba(25,:),'c')
plot([-49:10:trial_dur],aux_proba(50,:)-(aux_proba(100,:).*(aux_proba(100,:)>0)),'r')
hold off
axis([-49 trial_dur 0 60] )
title('Conditioned on neuron')
 

%%% compute and plot Tin and Tout trajectories conditioned on response
proba_out50_cor = proba_out50.*(perf_trial'*ones(1,trial_dur));
proba_out100_cor = proba_out100.*(perf_trial'*ones(1,trial_dur));

neuron_traj = sum(proba_out50_cor)/(sum(perf_trial));
antineuron_traj = sum(proba_out100_cor)/(sum(perf_trial));
%
plot(neuron_traj)
hold on
plot(mean(proba_out50),'b--')
plot(antineuron_traj,'g');
plot(mean(proba_out100),'g--')
hold off
axis([0 200 20 60])
ylabel('Activity (Hz)');
xlabel('Time (ms)');
title('Conditioned on correct')


%save exp_data spike_time50 spike_time100 perf_trial coherence RT
save datalast RT conf_level perf_trial perf mean_RT mean_RT_incor coherence neuron_traj antineuron_traj

beep
toc
