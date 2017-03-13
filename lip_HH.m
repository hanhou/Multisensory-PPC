%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LNP simulations for decision making
% Modified by HH 2017 @ UNIGE
% Adapted for the vestibular-visual multisensory heading discrimintation task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
clear
clear global

load h

global Y_targ X_train;
rand('state',sum(100*clock))

% %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === Switches ===
if_debug = 1;
if_test_set = 0;

if_bounded = 0; %% decis_flag=1 means that the trial stops at the bound (reaction time version)
decis_thres = 58; %% bound height

% === Sizes ===
% Input layers
N_vis = 100; % Visual motion signal
prefs_vis = linspace(-180,180,N_vis);

N_vest = 100; % Vestibular motion signal
prefs_vest = linspace(-180,180,N_vis);

N_targets = 2; % Target input

% Output layers
N_lip = 100;
prefs_lip = linspace(-180,180,N_lip);

% === Time ===
dt = 1e-3; % Size of time bin in seconds
trial_dur_total = 1700; %% in time bins (including pre_trial_dur)
stim_on_time = 200; % Time of motion start, in time bins.
N_trial = 100; % For each condition

ts = (0:trial_dur_total-1)*dt*1000; % in ms

if stim_on_time>trial_dur_total
    error('pre_trial_dur cannot be longer than trial_dur');
end

% === Stimuli ===
headings = [-8 -4 -2 -1 0 1 2 4 8];

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
t_motion = 0:trial_dur_total-stim_on_time; % in s
duration = 1500; num_of_sigma = 3.5; amp = 0.2;
miu = 750; sigma = duration/2/num_of_sigma; 
vel = exp(-(t_motion-miu).^2/(2*sigma^2));
vel = vel*amp/sum(vel*dt) ; % in m/s. Normalized to distance
acc = diff(vel)/dt; % in m/s^2

% To make sure t_motion, vel, and acc have the same length
t_motion(end) = []; 
vel(end) = []; 

% Gains
gain_vel_vis = 5; % (m/s)^-1
gain_acc_vest = 1.5; % (m^2/s)^-1

if if_debug
    figure(111);clf
    set(gcf,'name','Motion profile','uni','norm','pos',[0.632       0.381       0.358       0.403]);
    plotyy(t_motion,vel,t_motion,acc);
end


% === Network configuration ===

% -- Time constant for integration
time_const_lip = 1000e-3; % in s
time_const_vis  = time_const_lip;
time_const_vest = time_const_lip;

% -- Visual to LIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_lip_vis = 30;
K_lip_vis = 5;
dc_w_lip_vis = -5/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- Vestibular to LIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_lip_vest = 30;
K_lip_vest = 5;
dc_w_lip_vest = -5/100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Targets to LIP
g_w_lip_targ= 40;
K_lip_targ= 5;
att_gain2_targ = 0; % Drop in attention to visual target once motion stimulus appears.

% --- Recurrent connectivity in LIP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_w_lip_lip = 25;
K_lip_lip = 10;
dc_w_lip_lip = -10/100;

amp_I = 0;  % Mexico hat shape
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
%% Variable initialization (I moved here)
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

decision_ac = zeros(N_lip,trial_dur_total+2);

spikes_target = cell(2,1); spikes_vis = cell(2,1); spikes_vest = cell(2,1); rate_lip = cell(2,1); spikes_lip = cell(2,1);
spikes_count_vis = cell(2,1); spikes_count_vest = cell(2,1); proba_count_lip = cell(2,1); spikes_count_lip = cell(2,1);

for m = 1:if_test_set+1 % Train and Test
    
    % Raw data
    spikes_target{m} = nan(N_lip,trial_dur_total,N_trial,length(headings));
    spikes_vis{m} = nan(N_vis,trial_dur_total,N_trial,length(headings));
    spikes_vest{m} = nan(N_vest,trial_dur_total,N_trial,length(headings));

    rate_lip{m} = nan(N_lip,trial_dur_total,N_trial,length(headings));
    spikes_lip{m} = nan(N_lip,trial_dur_total,N_trial,length(headings));
    
    % Spike count in sliding windows for Fisher information
    spikes_count_vis{m} = nan(N_vis, fisher_N_window, N_trial,length(headings));
    spikes_count_vest{m} = nan(N_vest, fisher_N_window, N_trial,length(headings));
    proba_count_lip{m} = nan(N_lip, fisher_N_window, N_trial, length(headings));
    spikes_count_lip{m} = nan(N_lip, fisher_N_window, N_trial, length(headings));
    if shuffle_flag==1
        spikes_count_lip_shuff{m} = zeros(N_lip, fisher_N_window, N_trial, length(headings));
    end
    
end

% Output variables
RT = zeros(N_trial,length(headings));

%  Network connections 
w_lip_vis = zeros(N_lip,N_vis);
w_lip_vest = zeros(N_lip,N_vest);
w_lip_targ = zeros(N_lip,N_lip); % Not N_target, but N_lip (supposing the same number as the lip neurons)
w_lip_lip = zeros(N_lip,N_lip);

for nn=1:N_lip
    % -- Original coarse task --- 
    %     w_lip_vis(nn,:) = g_w_lip_vis/N_vis *(exp(K_lip_vis*(cos((prefs_vis-prefs_lip(nn))/360 *2*pi)-1)));  %  MT input 
    
    % -- Fine discrimination task --
    w_lip_vis(nn,:) = g_w_lip_vis/N_vis *(exp(K_lip_vis * (cos((prefs_lip(nn)-(-90*(prefs_vis<0)+90*(prefs_vis>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vis/180*pi))... % Gaussian(theta_lip - +/-90) * Sin(heading)
        + dc_w_lip_vis;   % Offset
    w_lip_vest(nn,:) = g_w_lip_vest/N_vest *(exp(K_lip_vest * (cos((prefs_lip(nn)-(-90*(prefs_vest<0)+90*(prefs_vest>0)))/180*pi)-1)))...
        .* abs(sin(prefs_vest/180*pi))... % Gaussian(theta_lip - +/-90) * Sin(heading)
        + dc_w_lip_vest;   % Offset

    % -- Targ->LIP and LIP->LIP. These are the same for fine/coarse task --
    w_lip_targ(nn,:) = g_w_lip_targ/N_lip *(exp(K_lip_targ*(cos((prefs_lip-prefs_lip(nn))/360 *2*pi)-1)));  %  Target input 
    w_lip_lip(nn,:) = g_w_lip_lip/N_lip*...   %  LIP recurrent 
        ((exp(K_lip_lip*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1)))-...
        amp_I*(exp(K_lip_I*(cos((prefs_lip-prefs_lip(nn))/360*2*pi)-1))))...
        + dc_w_lip_lip;
end

if if_debug
    figure(90); clf;
    set(gcf,'uni','norm','pos',[0.006       0.099       0.779       0.766],'name','Weights & Raster plot');
    subplot(3,3,7);
    imagesc(prefs_lip,prefs_vis,w_lip_vis'); colorbar; axis  xy; title('Vis->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{vis}'); ylim([-20 20]);
    set(gca,'xtick',-180:90:180);
    
    subplot(3,3,4);
    imagesc(prefs_lip,prefs_vest,w_lip_vest'); colorbar; axis  xy; title('Vest->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{vest}'); ylim([-20 20]);
    set(gca,'xtick',-180:90:180);
   
    subplot(3,3,1);
    % surf(prefs_lip,prefs_lip,w_lip_lip'); 
    imagesc(prefs_lip,prefs_lip,w_lip_lip'); 
    colorbar; axis  xy; title('LIP->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{lip}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
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

for ss = 1:length(headings)  % Motion directions
    
    % ==== Not trial- or time- dependent stuffs ===
    stim_this = headings(ss);
    fprintf('%dth  pos: %4.2f\n  ',ss, stim_this);
    
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

    if if_debug
        figure(91);clf;
        plot(prefs_vis,proba_vis/dt,'r');
        plot(prefs_vest,proba_vest/dt,'b');
        title('Sensory population tuning');
    end
    
    % proba_targ: proba of firing of target neurons in response to visual targets
    
    % pos_targ = stim_this + [0:360/N_targets:359.9];
    pos_targ = [-90 90]; % Always these two for fine discrimination task. HH2017

    proba_target_tuning = zeros(N_lip,1);
    for nn=1:N_targets
        proba_target_tuning = proba_target_tuning + ((max_rate_target-b_target)*...
            exp(K_target*(cos((prefs_lip'-pos_targ(nn))/360*2*pi)-1))+b_target)*dt;
    end
    proba_target_tuning = proba_target_tuning/(slope_norm_target*(N_targets/2)+dc_norm_target);   %  Divisive normalization 
    
    
    % === Begin trials ===
    for mm = 1:if_test_set+1 % Train and test sets
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
            % aux_proba_vis = proba_vis*[zeros(1,stim_on_time) ones(1,trial_dur_total-stim_on_time)]...
            %    + w_cov_vis*randn(N_vis,trial_dur_total);
            
            % Temporal dynamics added. HH2017
            aux_proba_vis = proba_vis*[zeros(1,stim_on_time) vel*gain_vel_vis]...
                + w_cov_vis*randn(N_vis,trial_dur_total);
           
            spikes_vis{mm}(:,1:trial_dur_total,tt,ss)= rand(N_vis,trial_dur_total)<(aux_proba_vis);
            
            % -- Vestibular ---
            % With the temporal gain of abs(acc)
            aux_proba_vest = proba_vest*[zeros(1,stim_on_time) abs(acc)*gain_acc_vest]...
                + w_cov_vest*randn(N_vest,trial_dur_total);
           
            spikes_vest{mm}(:,1:trial_dur_total,tt,ss)= rand(N_vest,trial_dur_total)<(aux_proba_vest);
            
            % == Override LIP dynamics during pre-trial by setting the spikes manually (weird thing) ==
            rate_lip{mm}(:,1:stim_on_time,tt,ss) = aux_proba_target(:,1:stim_on_time) / dt;
            spikes_lip{mm}(:,1:stim_on_time,tt,ss) = spikes_target{mm}(:,1:stim_on_time,tt,ss);
            
            aux1_rate_lip(:,stim_on_time) = rate_lip{mm}(:,stim_on_time,tt,ss);
            aux2_rate_lip(:,stim_on_time) = 0;
            
            % == Real network dynamics after stimlus onset ==
            k = stim_on_time; % Time begins at stim_on_time
            
            while k<=trial_dur_total-1
                
                % -- I don't understand yet why there are two variables aux1 and aux2 --
                aux1_rate_lip(:,k+1) = (1-dt/time_const_vis)*aux1_rate_lip(:,k)...   %  Self dynamics.  in Hz!
                    + 1/time_const_vis*(w_lip_vis*spikes_vis{mm}(:,k,tt,ss)...     %  Visual input 
                    + w_lip_vest*spikes_vest{mm}(:,k,tt,ss)...     % Vestibular input
                    + w_lip_targ * att_gain2_targ * spikes_target{mm}(:,k,tt,ss));             %  Target visual input (att_gain2 has been set to 0) 
                
                aux2_rate_lip(:,k+1) = (1-dt/time_const_lip)*aux2_rate_lip(:,k)+...  %  Self dynamics.  in Hz!
                    +1/time_const_lip*(w_lip_lip*spikes_lip{mm}(:,k,tt,ss));  %  LIP recurrent. in Hz!
                
                rate_lip{mm}(:,k+1,tt,ss) = aux1_rate_lip(:,k+1) + aux2_rate_lip(:,k+1) + bias_lip;  %  Bias term: b_out 
                
                % -- Note that proba_out is in Hz!! --
                spikes_lip{mm}(:,k+1,tt,ss)  = (rand(N_lip,1) < dt*(rate_lip{mm}(:,k+1,tt,ss)-threshold_lip));
                
                % -- Variable used for stopping the integration --
                
                if mm == 1 % Only apply to train set
                    decision_ac(:,k+1) = rate_lip{mm}(:,k+1,tt,ss);
                    %          decision_ac(:,k+1) = (1-delta_t/time_const_out)*decision_ac(:,k+1)+...
                    %                               +1/time_const_out*((w_oo-dc_w_oo)*spikes_out{1}(:,k));
                    
                    % -- Termination --
                    if ((range(if_bounded*decision_ac(:,k+1)) > decis_thres) && (RT(tt,ss)==0)) || (k==trial_dur_total-1)
                        RT(tt,ss) = k;
                        if k>fisher_wind_size-1   %  Enough data for Fisher info of at least one window 
                            last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,k-fisher_wind_size+1:k)')';
                        else
                            last_spikes__last_Fisher_win(:,count) = sum(spikes_lip{1}(:,1:k)')';
                        end
                        last_proba(:,count) = rate_lip{mm}(:,k,tt,ss);
                        k=trial_dur_total;  %  Terminate trial immediately
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

rate_real_vis_aver = nanmean(spikes_vis{1}(:,:,:,1),3)/dt;
rate_real_vest_aver = nanmean(spikes_vest{1}(:,:,:,1),3)/dt;
rate_expected_lip_aver = nanmean(rate_lip{1}(:,:,:,1),3);
rate_real_lip_aver = nanmean(spikes_lip{1}(:,:,:,1),3)/dt;

%{
%%  Animation
figure(1001); clf

for ttt = 1:round(trial_dur_total/200):trial_dur_total
    % LIP layer
    subplot(2,1,1);
    plot(prefs_lip,rate_expected_lip_aver(:,ttt),'r');  hold on;
    plot(prefs_lip,rate_real_lip_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',min(prefs_vis):90:max(prefs_vis));
    xlim([min(prefs_vis) max(prefs_vis)]); ylim([-20 150]);
    title(ttt);
    
    % Visual layer
    subplot(2,1,2);
    plot(prefs_vis,rate_real_vis_aver(:,ttt),'k');  hold off;
    set(gca,'xtick',min(prefs_vis):90:max(prefs_vis));
    xlim([min(prefs_vis) max(prefs_vis)]); ylim([-20 150]);
    
    drawnow;
end


%}

%%{
%% Raster plot
if if_debug
    figure(90);
    
    subplot(3,3,[2 3]);
    imagesc(ts,prefs_lip,rate_real_lip_aver*dt); axis xy;
    ylabel('LIP');  title(sprintf('Firing prob. for each bin, averaged over %g trials',N_trial));  colorbar
    
    subplot(3,3,[5 6]);
    imagesc(ts,prefs_vest,rate_real_vest_aver*dt); axis xy; colorbar
    ylabel('Vest');
    
    subplot(3,3,[8 9]);
    imagesc(ts,prefs_vis,rate_real_vis_aver*dt); axis xy; colorbar
    ylabel('Vis');
    
    colormap hot;
end
%}

%%{
%%
set(figure(1002),'name','PSTH','uni','norm','pos',[0.005 0.044 0.767 0.438]); clf;
subplot(1,2,1);

to_plot_heading = length(headings);
[~,pref_ind] = min(abs(prefs_lip - 90)); % Left and right for the heading task
[~,null_ind] = min(abs(prefs_lip - -90));

for cc = 1:round(N_trial/10):N_trial
    plot(ts,rate_lip{1}(pref_ind,:,cc,to_plot_heading),'r','linewid',2);
    hold on; plot(rate_lip{1}(null_ind,:,cc,to_plot_heading),'k--','linewid',1);
end

plot(t_motion+stim_on_time,vel/max(vel)*max(ylim)/3,'k--');

if if_bounded
    plot(xlim,[decis_thres decis_thres],'k--');
end

title(sprintf('rate\\_lip, heading = %g, coh = %g',headings(to_plot_heading),coherence));

subplot(1,2,2);
plot(ts,nanmean(rate_lip{1}(pref_ind,:,:,to_plot_heading),3)-nanmean(rate_lip{1}(null_ind,:,:,to_plot_heading),3),'r');
hold on; plot(t_motion+stim_on_time,vel/max(vel)*max(ylim)/3,'k--');
title(sprintf('averaged of %g trials, pref-null',N_trial));

%}

%% Behavior

rate_lip_at_decision = squeeze(rate_lip{1}(:,end,:,:));
[~, pos_max_rate_lip_at_decision] = max(rate_lip_at_decision,[],1);
choices = prefs_lip(squeeze(pos_max_rate_lip_at_decision)) >= 0; % 1 = rightward, 0 = leftward
figure(99); hold on;
plot(headings,sum(choices,1)'/N_trial,'go','markersize',15); % Psychometric
xxx = min(headings):0.1:max(headings);

[bias, threshold] = cum_gaussfit_max1([headings' sum(choices,1)'/N_trial ]);
plot(xxx,normcdf(xxx,bias,threshold),'g-','linewid',4);

set(text(-8,0.6,sprintf('threshold = %g\n',threshold)),'color','g');                    
   
                  






