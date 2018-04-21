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

%% ========== Calculate time-sliding tuning curves ===========

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
    pref_this_k = pref_az(k,end);
    [~,ind] = min(abs(pref_this_k - unique_azimuth));
    mean_firing_matrix_wrap(:,:,k) = circshift(mean_firing_matrix(:,:,k),length(unique_azimuth)/2+1-ind);
    se_firing_matrix_wrap(:,:,k) = circshift(se_firing_matrix(:,:,k),length(unique_azimuth)/2+1-ind);
end

mean_firing_matrix_wrap(end+1,:,:) = mean_firing_matrix_wrap(1,:,:);
se_firing_matrix_wrap(end+1,:,:) = se_firing_matrix_wrap (1,:,:);

figure(223); clf;      set(gcf,'uni','norm','pos',[0.023        0.34       0.855       0.365]);
for k = 1:length(unique_stim_type)
    subplot(1,3,k);
    errorbar(mean_firing_matrix_wrap(:,:,k),se_firing_matrix_wrap(:,:,k));
    hold on; 
    plot(mean_firing_matrix_wrap(:,end,k),'k','linew',2);
    title(['p=' num2str(p_value(k))]);
end
set(gcf,'name',FILE);

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
    pref_az, ROI, p_value ...
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

