%----------------------------------------------------------------------------------------------------------------------
%-- PSTH.m -- Plots Post Stimulus Time Histogram for MOOG 3D tuning expt
%--	Yong, 6/27/03
% Hou, Han, To show that MSTd's visual and vestibular responses obey ilPPC. 20180223
%-----------------------------------------------------------------------------------------------------------------------

function Azimuth_PSTH_HH_for_ilPPC(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, batch_flag);

TEMPO_Defs;
Path_Defs;

%% ========== Get data =============

trials = 1:size(data.moog_params,2);		% a vector of trial indices

% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
select_trials = false(length(trials),1);
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    keyboard;
end

NULL_TRIAL = -9999;

stim_type_per_trial = data.moog_params(STIM_TYPE,select_trials,MOOG)';
azimuth_per_trial   = data.moog_params(AZIMUTH, select_trials, MOOG)';
elevation_per_trial =  data.moog_params(ELEVATION, select_trials, MOOG)';

unique_stim_type = setxor(NULL_TRIAL, munique(stim_type_per_trial));
unique_azimuth = setxor(NULL_TRIAL, munique(azimuth_per_trial));
unique_elevation = setxor(NULL_TRIAL, munique(elevation_per_trial));

% -- Time information
eye_timeWin = 1000/(data.htb_header{EYE_DB}.speed_units/data.htb_header{EYE_DB}.speed/(data.htb_header{EYE_DB}.skip+1)); % in ms
spike_timeWin = 1000/(data.htb_header{SPIKE_DB}.speed_units/data.htb_header{SPIKE_DB}.speed/(data.htb_header{SPIKE_DB}.skip+1)); % in ms
event_timeWin = 1000/(data.htb_header{EVENT_DB}.speed_units/data.htb_header{EVENT_DB}.speed/(data.htb_header{EVENT_DB}.skip+1)); % in ms

% -- Spike data
spike_in_bin = squeeze(data.spike_data(SpikeChan,:,select_trials))';   % TrialNum * 5000
spike_in_bin( spike_in_bin > 100 ) = 1; % something is absolutely wrong

% -- Event data
event_in_bin = squeeze(data.event_data(:,:,select_trials))';  % TrialNum * 5000

% Get Coherence   20180603
coherence = mean(data.moog_params(COHERENCE,:,1));

%% ========== 1. Calculate time-sliding tuning curves ===========

% ------- Time-related -------
trial_begin = mode(mod(find(event_in_bin'==4),5000));
trial_end = mode(mod(find(event_in_bin'==5),5000));
ROI = [ 
    -300 -100;
    -100 100;
    100 300;
    300 500;
    500 700;  % Region of interests, in ms
    700 900;
    900 1100;
    1100 1300;
    1300 1500;
    1500 1700;
    1700 1900;
    1900 2100;
    2100 2300;
    (trial_end-trial_begin)/2-750 (trial_end-trial_begin)/2+750];  % Classical defination as Gu

% ------- Condition related --------
elevation_included = [0]; % [-45 0 45]

mean_firing_matrix = nan(length(unique_azimuth), size(ROI,1), length(unique_stim_type));
se_firing_matrix = nan(length(unique_azimuth), size(ROI,1), length(unique_stim_type));
p_value = nan(3,1); % Only classical time range

%---------------------------------------------
for k = 1:length(unique_stim_type)
    for tt = 1:length(ROI)
        ROI_bins = trial_begin + round(ROI(tt,1)/spike_timeWin) : trial_begin + round(ROI(tt,2)/spike_timeWin);
        for aa = 1:length(unique_azimuth)
            select_trials = (stim_type_per_trial == unique_stim_type(k)) & ...
                (azimuth_per_trial == unique_azimuth(aa)) & ...
                any(elevation_per_trial == elevation_included,2);
            
            if tt == length(ROI)
                raw_firing_for_p_value(aa,1:sum(select_trials),k) = sum(spike_in_bin(select_trials,ROI_bins),2);
            end
            
            mean_firing_matrix(aa,tt,k) = mean(sum(spike_in_bin(select_trials,ROI_bins),2),1)/(range(ROI(tt,:))/1000);
            se_firing_matrix(aa,tt,k) = std(sum(spike_in_bin(select_trials,ROI_bins),2),1)/(range(ROI(tt,:))/1000)/sqrt(sum(select_trials));
        end
        
        % Calculate preferred direction
        [pref_az(k,tt), el, amp] = vectorsumAngle(mean_firing_matrix(:,tt,k), unique_azimuth, unique_azimuth * 0);
        
    end
    %     pref_az
    %     figure(223); subplot(1,3,k);
    %     errorbar(repmat(unique_azimuth,1,size(ROI,1)),mean_firing_matrix(:,:,k),se_firing_matrix(:,:,k));
    %     polarplot([unique_azimuth/180*pi; unique_azimuth(1)/180*pi],...
    %         [mean_firing_matrix(:,:,k); mean_firing_matrix(1,:,k);]);
    
    % P value
    p_value(k) = anova1(raw_firing_for_p_value(:,:,k)','','off');
    
    
    % Wrap tuning curve so that pref is at the center
    % The caveats is that the sampling should be uniform here!! (or using curve fitting first which I'm not willing to do so)
    pref_this_k = pref_az(k,end);
    uniform_unique_azimuth = 0:45:315; % Downsample azimuths    
    down_sampled_azimuth = ismember(unique_azimuth, uniform_unique_azimuth);
    uniform_mean_firing_matrix = mean_firing_matrix(down_sampled_azimuth,:,k);
    uniform_se_firing_matrix = se_firing_matrix(down_sampled_azimuth,:,k);
        
    [~,ind(k)] = min(abs(pref_this_k - uniform_unique_azimuth));
    mean_firing_matrix_wrap(:,:,k) = circshift(uniform_mean_firing_matrix, length(uniform_unique_azimuth)/2+1-ind(k));
    se_firing_matrix_wrap(:,:,k) = circshift(uniform_se_firing_matrix, length(uniform_unique_azimuth)/2+1-ind(k));
end

mean_firing_matrix_wrap(end+1,:,:) = mean_firing_matrix_wrap(1,:,:);
se_firing_matrix_wrap(end+1,:,:) = se_firing_matrix_wrap (1,:,:);

% Congruency index
r = corrcoef(mean_firing_matrix(:,end,1),mean_firing_matrix(:,end,2));
CI.all = r(2); % Simple one

r = corrcoef(sum(mean_firing_matrix(:,3:7,1),2), sum(mean_firing_matrix(:,3:7,2),2));
CI.first_half = r(2);

r = corrcoef(sum(mean_firing_matrix(:,8:12,1),2), sum(mean_firing_matrix(:,8:12,2),2));
CI.second_half = r(2);


color_order = colormap(jet);
color_order = color_order(round(linspace(1,size(color_order,1),length(ROI))),:);

set(0,'defaultaxescolororder',color_order);

figure(223); clf;      set(gcf,'uni','norm','pos',[0.023       0.079       0.776       0.763]);
for k = 1:length(unique_stim_type)
    subplot(2,3,k);
    errorbar(mean_firing_matrix_wrap(:,:,k),se_firing_matrix_wrap(:,:,k));
    hold on; 
    plot(mean_firing_matrix_wrap(:,end,k),'k','linew',2);
    title(sprintf('stim type = %g, p=%g', unique_stim_type(k), p_value(k)));
    
    % Plot "90 heading"
    pos_90 = 3 + (length(uniform_unique_azimuth)/2+1-ind(k));
    plot([pos_90 pos_90], ylim, 'k--', 'linew',2)
    
    if unique_stim_type(k)==2,  ylabel(sprintf('Coherence = %g',coherence)), end
   
end
set(gcf,'name',FILE);

%% ========== 2. 10-ms window PSTH of Pref-Null to fitting Model M1 (20180916) ===========
% Since only temporal dynamics matters in fitting the traces, I only care about the relative responses.
% Given that, I 
%   1) Find out the pref direction (using Gu's center 1.5 s window)
%   2) Output the PSTH along this direction (PREF) as well as the deltaPSTH between pref and null (PREF-NULL) as a function of time.

% Same as for LIP data
binSize_rate = 10;  % in ms
stepSize_rate = 10; % in ms
smoothFactor = 50; % in ms !!

PSTH_ts = -150 : stepSize_rate : 2150; % Only 0~2000 ms
PSTH_pref = nan(length(PSTH_ts),length(unique_stim_type));
PSTH_null = nan(length(PSTH_ts),length(unique_stim_type));

for k = 1:length(unique_stim_type)
    % 1) Pref direction (simply using max instead of vector sum)
    [~,pref_this] = max(mean_firing_matrix(:,end,k)); % "End" means center 1.5s window
    null_this = unique_azimuth == mod(unique_azimuth(pref_this) + 180, 360);
    
    % 2) Get PSTH from raw spike train
    select_pref_trials = (stim_type_per_trial == unique_stim_type(k)) & ...
        (azimuth_per_trial == unique_azimuth(pref_this)) & ...
        any(elevation_per_trial == elevation_included,2);
    select_null_trials = (stim_type_per_trial == unique_stim_type(k)) & ...
        (azimuth_per_trial == unique_azimuth(null_this)) & ...
        any(elevation_per_trial == elevation_included,2);
    
    % --- Calculate PSTH --  HH20180916
    for tt = 1:length(PSTH_ts)
        this_t = trial_begin + (PSTH_ts(tt) - binSize_rate/2)/spike_timeWin + 1 : trial_begin + (PSTH_ts(tt) + binSize_rate/2)/spike_timeWin;
        PSTH_pref(tt,k) = mean(sum(spike_in_bin(select_pref_trials,this_t),2),1)/(binSize_rate/1000);
        PSTH_null(tt,k) = mean(sum(spike_in_bin(select_null_trials,this_t),2),1)/(binSize_rate/1000);
    end
    
    PSTH_pref_smoothed(:,k) =  GaussSmooth(PSTH_ts,PSTH_pref(:,k),smoothFactor);
    PSTH_null_smoothed(:,k) =  GaussSmooth(PSTH_ts,PSTH_null(:,k),smoothFactor);
    
    % --- Calculate FI --- (simple Gu) HH20180917
    
end

set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);
figure(223); subplot(2,3,4)
plot(PSTH_ts,PSTH_pref_smoothed - PSTH_null_smoothed)

%% ===== 3. Fisher info (20180917) ======

binSize_FI = 250;  % in ms  % According to Shadlen 2001
stepSize_FI = 50; % in ms  % Match PSTH, CD, and others

FI_ts = -150 : stepSize_FI : 2150; 
FI = nan(length(FI_ts),length(unique_stim_type));

mean_firing_matrix_for_FI = nan(length(unique_azimuth), length(FI_ts), length(unique_stim_type));
se_firing_matrix_for_FI = nan(length(unique_azimuth), length(FI_ts), length(unique_stim_type));

xx = 0:0.1:360; % Gu 2010
xx_left = find(xx<90,1,'last');
xx_right = find(xx>90,1);

for k = 1:length(unique_stim_type)
    
    % --- Calculate PSTH --  HH20180917
    for tt = 1:length(FI_ts)
        this_t = trial_begin + (FI_ts(tt) - binSize_FI/2)/spike_timeWin + 1 : trial_begin + (FI_ts(tt) + binSize_FI/2)/spike_timeWin;
        
        % -- Get spatial tuning --
        for aa = 1:length(unique_azimuth)
            select_trials = (stim_type_per_trial == unique_stim_type(k)) & ...
                (azimuth_per_trial == unique_azimuth(aa)) & ...
                any(elevation_per_trial == elevation_included,2);
            
            mean_firing_matrix_for_FI(aa,tt,k) = mean(sum(spike_in_bin(select_trials,this_t),2),1)/(binSize_FI/1000);
            se_firing_matrix_for_FI(aa,tt,k) = std(sum(spike_in_bin(select_trials,this_t),2),1)/(binSize_FI/1000)/sqrt(sum(select_trials));
        end
        
        % -- Calculate FI --
        spline_this = spline(unique_azimuth, mean_firing_matrix_for_FI(:,tt,k), xx); % Spline fit
        slope_at_90 = (spline_this(xx_left)-spline_this(xx_right))/(xx(xx_left)-xx(xx_right)); % Local slope
        variance_at_90 = max(spline_this((xx_left+xx_right)/2), 0.5); % Poisson assumption with 0.5 Hz minimal variance, Gu 2010
        
        % --- Save data ---
        FI(tt,k) = slope_at_90^2/variance_at_90;
        
    end
        
end

figure(223); subplot(2,3,5)
plot(FI_ts, FI)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Saving

% Reorganized. HH20141124
config.batch_flag = batch_flag;

% Output information for test. HH20160415
if isempty(batch_flag)
    config.batch_flag = 'test.m';
    disp('Saving results to \batch\test\ ');
end

%%%%%%%%%%%%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = PackResult(FILE, PATH, SpikeChan, unique_stim_type, Protocol, ... % Obligatory!!
    unique_azimuth, unique_elevation,...
    mean_firing_matrix_wrap, se_firing_matrix_wrap, ...
    mean_firing_matrix, se_firing_matrix, ...
    pref_az, ROI, p_value, coherence, ...
    PSTH_pref_smoothed, PSTH_null_smoothed, PSTH_ts, ...
    FI_ts, FI, CI ...
    ); % model info

config.suffix = 'ilPPC';

% figures to save
config.save_figures = [223];

% Only once
config.sprint_once_marker = '';
config.sprint_once_contents = '';

% loop across stim_type
config.sprint_loop_marker = '';
config.sprint_loop_contents = '';

config.append = 1; % Overwrite or append

SaveResult(config, result);

return;

