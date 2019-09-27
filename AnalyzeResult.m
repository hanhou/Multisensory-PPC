%% Analysis switch
if ION_cluster
analysis_switch = [0;  % 0p1 Fig2a of Beck 2008 and noise correlation calculation 
                   1;  % 0p5 Decoding method 
                       % 1 Overview (mandatary)
                   1;  % 1p5 Overview (normalized)
                   1;  % 2 Example
                   1;  % 3 Cells deltaPSTH
                   1;  % 3p5 Cells rawPSTH
                   1;  % 4 Heterogeneity
                   1;  % 5 Linear regression (Gu 2008)
                   1;  % 6 Information in sensory areas
                   ];
else
analysis_switch = [0;%1;  % 0p1 Fig2a of Beck 2008 and noise correlation calculation
                   0;  % 0p5 Decoding method 
                       % 1 Overview (mandatary)
                   0;  % 1p5 Overview (normalized)
                   0;  % 2 Example
                   0;%1;  % 3 Cells deltaPSTH
                   0;% 1;  % 3p5 Cells rawPSTH
                   0;%1;  % 4 Heterogeneity
                   0;%1;  % 5 Linear regression (Gu 2008)
                   1;%1;  % 6 Information in sensory areas
                   ];
end

%% Get data
% if ~exist('RT','var')  % Used in offline analysis
%     clear;
%     close all
%     disp('Offline...');
%     load('./result/last_result.mat');
%     
%     if_debug = 0;
% end

if ~exist('trial_dur_total_in_bins','var')  % Used in offline analysis
    % Unpack all paras
    para_list = fieldnames(paras);
    for i = 1:length(para_list)
        eval([para_list{i} '= paras.' para_list{i} ';']);
    end

    if_debug = 0; % Override
end

unique_abs_heading = unique(abs(unique_heading));
analysis_time = tic;
colors = [0 0 1; 1 0 0; 0 0.8 0.4];
LEFT = -1; RIGHT = 1;

if ~exist('para_override_txt','var')
    para_override_txt = '';
end

%% === Rectify rate_lip ===
% I use rectified rate_lip because
% 1. The mean of which is equal to the mean firing rate calculated from the real spike train.
% 2. The variance of which corresponds to VarCE of the real spike trains, without the need to get rid of the Poisson variability

rate_lip(rate_lip<0) = 0;
rate_int(rate_int<0) = 0;
% rate_vest(rate_vest<0) = 0;
% rate_vis(rate_vis<0) = 0;
ts = ((0:trial_dur_total_in_bins-1)-stim_on_time_in_bins)*dt; % in s

% I don't use these two any more. Just for backup. HH20180621
% %{
rate_int_Ann = rate_int;
rate_lip_Ann = rate_lip;
ts_Ann = ts;
win_step = dt;
%}

%% == Actual shifting window spike count for LIP as in experimental data (~ rectified rate_lip etc.) ===
%{ 
% I don't use this. See above.  % No, I decide to use this. HH20180621
fprintf('Generating PSTH (same as experiments)...');

% %{
% if if_debug

% I decide to use EXACTLY the same way as in the experiments:
%    binSize 10ms, stepSize 10ms sliding windows + 50ms Gaussian smoother
% GaussSmooth has the "dropping at edges" problem, so I use conventional smooth methods, which, in turn, 
% equvalent to increase the win_wid to ~50e-3 s without smooth (NO... there would be small oscillations. Well, add a small smoother.)
win_wid = 50e-3; % in s    % Should be larger enough than dt, otherwise there is another "time-window-induced-low-firing-rate" problem...
win_step = 10e-3;
smooth_win = 2*win_step; % GaussSmooth has the "dropping at edges" problem, so I use conventional smooth

ts_Ann = ts(1) + win_wid/2: win_step: ts(end) - win_wid/2; % New time center for rate
rate_lip_Ann = nan(N_lip,length(ts_Ann),N_rep,length(unique_heading),length(unique_stim_type)); % Ann Churchland's version
rate_int_Ann = rate_lip_Ann;

% -- Sliding window --
for tt = 1:length(ts_Ann)
    % fprintf('%.2g%%\n',tt/length(ts_Ann)*100);
    this_range = ts_Ann(tt)-win_wid/2 <= ts & ts <= ts_Ann(tt)+win_wid/2;
    rate_lip_Ann(:,tt,:,:,:) = sum(spikes_lip(:,this_range,:,:,:),2)/win_wid;
    rate_int_Ann(:,tt,:,:,:) = sum(spikes_int(:,this_range,:,:,:),2)/win_wid;
end

oriSize = size(rate_lip_Ann);
tmp_rate_lip_Ann = reshape(shiftdim(rate_lip_Ann,1),length(ts_Ann),[]);

% -- Gauss smooth --
nSmooth = round(smooth_win/win_step);
parfor tt = 1:size(tmp_rate_lip_Ann,2)
    tmp_rate_lip_Ann(:,tt) = smooth(tmp_rate_lip_Ann(:,tt),nSmooth);
end

rate_lip_Ann_smooth = shiftdim(reshape(tmp_rate_lip_Ann,oriSize(2),oriSize(3),oriSize(4),oriSize(5),oriSize(1)),4);

fprintf(' Done!\n');
%} 

if if_debug
    
    figure(1041); clf
    plot(ts_Ann,mean(rate_lip_Ann(left_targ_ind,:,:,1,2),3),'b'); hold on;
    plot(ts_Ann,mean(rate_lip_Ann(right_targ_ind,:,:,1,2),3),'b--')

%     plot(ts_Ann,mean(rate_lip_Ann_smooth(left_targ_ind,:,:,1,2),3),'k');
%     plot(ts_Ann,mean(rate_lip_Ann_smooth(right_targ_ind,:,:,1,2),3),'k--')

    % plot(ts,mean(rectified_rate_lip(left_targ_ind,:,:,end,1),3),'r--');
    % plot(ts,mean(rectified_rate_lip(right_targ_ind,:,:,end,1),3),'r');
    plot(ts,mean(rate_lip(left_targ_ind,:,:,1,2),3),'g');
    plot(ts,mean(rate_lip(right_targ_ind,:,:,1,2),3),'g--');
    % export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/test.png'));
    
    
    plot(ts_Ann,mean(rate_int_Ann(left_targ_ind,:,:,1,2),3),'b'); hold on;
    plot(ts_Ann,mean(rate_int_Ann(right_targ_ind,:,:,1,2),3),'b--')

    plot(ts,mean(rate_int(left_targ_ind,:,:,1,2),3),'k');
    plot(ts,mean(rate_int(right_targ_ind,:,:,1,2),3),'k--')


    
    %% Example snapshot of population activity (for demo)
    spkCntCent = [0.8];
    spkCntWin = 0.2; % 200 ms
    areas = {'lip','int','target','vest','vis'};
    cols = {'g','k','k','b','r'};
    figure(15); clf
    prefs_target = prefs_lip;
    
    for aa = 1:length(areas) % Draw the last condition
        eval(sprintf('real_count_%s = sum(spikes_%s(:,spkCntCent-spkCntWin/2 <= ts & ts < spkCntCent+spkCntWin/2,end,end,end),2);',areas{aa},areas{aa}));
        subplot(2,3,aa); 
        eval(sprintf('plot(prefs_%s,real_count_%s,''o%s''); title(''%s'');',areas{aa},areas{aa},cols{aa},areas{aa}));
        if aa >= 4
            hold on; plot([0 0],[0 max(ylim)],'k--');
        end
        ylabel('Spike count');
    end
end
%}

%% == Plotting example network activities ==
this_heading_ind = length(unique_heading);
to_plot_cond = length(unique_stim_type);

if if_debug
%%  ====== Animation ======
    rate_real_vis_aver = nanmean(spikes_vis(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_vest_aver = nanmean(spikes_vest(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_expected_int_aver = nanmean(rate_int_Ann(:,:,:,this_heading_ind,to_plot_cond),3);
    rate_expected_lip_aver = nanmean(rate_lip_Ann(:,:,:,this_heading_ind,to_plot_cond),3);
    rate_real_int_aver = nanmean(spikes_int(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_lip_aver = nanmean(spikes_lip(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_targ_aver = nanmean(spikes_target(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    
    figure(1001); clf
    
    for ttt = 1:10:length(ts_Ann)
        
        % INTEGRATOR layer
        subplot(2,2,1);
        plot(prefs_int,rate_expected_int_aver(:,ttt)); hold on;
        plot(prefs_int,rate_real_int_aver(:,ttt),'r');  hold off;
        axis(gca,[min(prefs_int) max(prefs_int) min(rate_expected_int_aver(:)) max(rate_expected_int_aver(:))*2]);
        title(sprintf('Integrator, heading = %g, cond = %g, aver %g trials',unique_heading(this_heading_ind),unique_stim_type(to_plot_cond),N_rep));
        
        % LIP layer
        subplot(2,2,2);
        plot(prefs_lip,rate_expected_lip_aver(:,ttt)); hold on;
        plot(prefs_lip,rate_real_lip_aver(:,ttt),'r');  hold off;
        axis(gca,[min(prefs_lip) max(prefs_lip) min(rate_expected_lip_aver(:)) max(rate_expected_lip_aver(:))*2+0.1]);
        title(sprintf('LIP, t = %g',ttt*dt*1000));
        
        % Visual layer
        subplot(2,2,3);
%         plot(pref_vis,aux_proba_vis(:,ttt)/dt); hold on;
        plot(prefs_vis,rate_real_vis_aver(:,ttt),'r'); hold off;
%         axis(gca,[min(prefs_int) max(prefs_int) min(aux_proba_vis(:)/dt) max(aux_proba_vis(:)/dt)*2]);
        title('Visual');
        
        % Target layer
        subplot(2,2,4);
        plot(prefs_lip,proba_target_tuning/dt); hold on;
        plot(prefs_lip,rate_real_targ_aver(:,ttt),'r'); hold off;
        axis(gca,[min(prefs_lip) max(prefs_lip) min(proba_target_tuning(:)/dt) max(proba_target_tuning(:)/dt)*2]);
        title('Target');
        
        drawnow;
    end
else % Not debug, Fig.2a in Beck 2008
  %% Averaged population activity for the to_plot_cond condition
  %   %{
  if analysis_switch(1)
      rate_expected_lip_aver = nanmean(rate_lip_Ann(:,:,:,this_heading_ind,to_plot_cond),3);
      figure(1647);  set(gcf,'uni','norm','pos',[0.252       0.431       0.669       0.435]);    clf;
      
      subplot(1,2,1);
      to_snap_time = [0.25 0.5 0.75 1.0 1.25];
      snapshots =  arrayfun(@(x)find(ts_Ann>x,1),to_snap_time);
      set(0,'defaultaxescolororder',lines)
      plot(prefs_lip,rate_expected_lip_aver(:,snapshots),'o');
      set(gca,'xtick',[-180:90:180]);
      title('PSTH');
      legend(num2str(to_snap_time'));
      
%       subplot(1,3,2);
%       plot(prefs_lip(1:end/2),rate_expected_lip_aver(end/2+1:end,snapshots)-fliplr(rate_expected_lip_aver(1:end/2,snapshots)),'o');
%       title('Diff PSTH');
%       set(gca,'xtick',-180:90:180);
      
      %% Verifying noise correlation. HH20171127
      spikecount_vest_each_trial = squeeze(reshape(sum(spikes_vest(:,ts_Ann>stim_on_time,:,:,[1 3]),2),N_vest,N_rep,[])); % All spikes
      noise_vest = reshape(spikecount_vest_each_trial - repmat(mean(spikecount_vest_each_trial,2),1,N_rep,1),N_vest,[]); % Get noise
      noise_correlation_vest = corrcoef(noise_vest');
            
      figure(1647);       subplot(1,2,2); 
      noise_correlation_vest(logical(eye(N_vest))) = nan;
      imagesc(prefs_vest, prefs_vest, noise_correlation_vest);
      
      % Follow Beck, 2008
      aux_mask = ones(N_vest,1)*[1:N_vest];
      mask = ((1-cos(abs(aux_mask-aux_mask')/N_vest*2*pi))/2)<.5;
      mask = mask.*(1-eye(N_vest));

      mean_corr_less_than_90 = nanmean(noise_correlation_vest(logical(mask)));
      title(['during stimuli, mean corr = ' num2str(mean_corr_less_than_90)]);

      %       % Check noise correlation before stimuli onset
      %       subplot(1,2,1);
      %       spikecount_vest_before_onset = squeeze(reshape(sum(spikes_vest(:,ts<=stim_on_time,:,:,[1 3]),2), N_vest,N_rep,[])); % All spikes
      %       noise_vest_before_onset = reshape(spikecount_vest_before_onset,N_vest,[]); % Noise = spikes count itself
      %       noise_correlation_vest_before_onset = corrcoef(noise_vest_before_onset');
      %
      %       noise_correlation_vest_before_onset(logical(eye(N_vest)))= nan;
      %       imagesc(prefs_vest, prefs_vest, noise_correlation_vest_before_onset);
      %
      %       mean_corr_less_than_90_before_onset = nanmean(noise_correlation_vest_before_onset(logical(mask)));
      %       title(['before stimuli, mean corr = ' num2str(mean_corr_less_than_90_before_onset)]);

      
      if ION_cluster
          fprintf('...Saving fig 0p1...');
          export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s0p1_Fig2a%s.png',save_folder,para_override_txt));
          saveas(gcf,sprintf('./result/%s0p1_Fig2a%s.fig',save_folder,para_override_txt),'fig');
      end
      
      
  end
  %}
end

%% ====== Connection and Raster plot ======
if if_debug
    figure(90);
    
    subplot(4,3,[2 3]);
    imagesc(ts_Ann,prefs_lip,rate_real_lip_aver*dt); axis xy;
    ylabel('LIP');
    title(sprintf('Firing prob. for each bin, averaged over %g trials, cond = %g, heading = %g',...
        N_rep,unique_stim_type(to_plot_cond),unique_heading(this_heading_ind)));  colorbar
    
    subplot(4,3,[5 6]);
    imagesc(ts_Ann,prefs_int,rate_real_int_aver*dt); axis xy; colorbar
    ylabel('INTEGRATOR');
    
    subplot(4,3,[8 9]);
    imagesc(ts_Ann,prefs_vest,rate_real_vest_aver*dt); axis xy; colorbar
    ylabel('Vest');
    
    subplot(4,3,[11 12]);
    imagesc(ts_Ann,prefs_vis,rate_real_vis_aver*dt); axis xy; colorbar
    ylabel('Vis');
end
%}

%% ====== Readout LIP activity and make decisions =====
disp('Decoding...');

if read_out_at_the_RT == 1 % Freeze population activity at RT for each trial
    rate_lip_at_decision = nan(N_lip,N_rep,length(unique_heading),length(unique_stim_type));
    for aa = 1:N_rep
        for bb = 1:length(unique_heading)
            for cc = 1:length(unique_stim_type)
                try
                    bin_in_Ann = find(ts_Ann >= RT(aa,bb,cc) - dt, 1, 'first');
                    rate_lip_at_decision(:,aa,bb,cc) = rate_lip_Ann(:,bin_in_Ann,aa,bb,cc);
                catch
                    keyboard
                end
            end
        end
    end
else
    rate_lip_at_decision = squeeze(rate_lip_Ann(:,end,:,:,:)); % Or, get the population activity at the end of the trials
end

% -------- Simple readout: Location of max rate -------------
[~, readout_lip_at_decision] = max(rate_lip_at_decision,[],1); % Readout population activity    
choices_maxpos = sign((prefs_lip(shiftdim(readout_lip_at_decision)) >= 0) -0.5); % choices(reps,headings,stimtypes):  1 = rightward, -1 = leftward
correct_rate_maxpos = sum(choices_maxpos(:) == correct_ans_all(:)) / numel(correct_ans_all);

%     % I just flip the psychometric curve to the negative headings
%     psychometric = [unique_heading' sum(reshape(choices{ss},[],length(unique_heading)),1)'/N_rep];
%     fake_psy = flipud(psychometric(unique_heading>0,:));
%     fake_psy(:,1) = -fake_psy(:,1);
%     fake_psy(:,2) = 1-fake_psy(:,2);
%     psychometric = [fake_psy ; psychometric];

% ------------ SVM decoder --------------
train_ratio = 0.3;

% To make sure each condition always has the same number of trials.
n_train_each_condition = round(N_rep * train_ratio); 
n_test_each_condition = N_rep - n_train_each_condition;
n_train = n_train_each_condition * length(unique_heading) * length(unique_stim_type); 
n_test = n_test_each_condition * length(unique_heading) * length(unique_stim_type);

rate_lip_to_train = nan(n_train,N_lip);  correct_ans_train = nan(n_train,1);
rate_lip_to_test = nan(n_test,N_lip);  correct_ans_test = nan(n_test,1);

% Generate train and test sets
n = 0;
for ss = 1:length(unique_stim_type)
    for hh = 1:length(unique_heading)
        n = n + 1;
        train_this = randperm(N_rep,n_train_each_condition);
        test_this = setdiff(1:N_rep,train_this);
        
        where_train_this = (n-1)*n_train_each_condition + 1 : n*n_train_each_condition;
        rate_lip_to_train(where_train_this,:) = rate_lip_at_decision(:,train_this,hh,ss)';
        correct_ans_train(where_train_this,1) = correct_ans_all(train_this,hh,ss);
        
        where_test_this = (n-1)*n_test_each_condition + 1 : n*n_test_each_condition;
        rate_lip_to_test(where_test_this,:) = rate_lip_at_decision(:,test_this,hh,ss)';
        correct_ans_test(where_test_this,1) = correct_ans_all(test_this,hh,ss);        
    end
end

% Train linear SVM, try different Cs
% Cs = 10.^(-3:0.5:0);
Cs = 1e-5;
correct_rate_svm_test_Cs = nan(1,length(Cs));

for ccs = 1:length(Cs)
%     try
        svm_model_Cs{ccs} = svmtrain(rate_lip_to_train,correct_ans_train,'boxconstraint',Cs(ccs),'tol',1e-7);
        choices_svm_test_Cs{ccs} = svmclassify(svm_model_Cs{ccs},rate_lip_to_test);
        correct_rate_svm_test_Cs(ccs) = sum(choices_svm_test_Cs{ccs} == correct_ans_test) / length(correct_ans_test);
%     end
end

% Select the best C
fprintf('Correct rate of SVM with different Cs:\n%s\n',num2str(correct_rate_svm_test_Cs));
[~,bestC] = max(correct_rate_svm_test_Cs);
svm_model = svm_model_Cs{bestC};
choices_svm_test = choices_svm_test_Cs{bestC};
correct_rate_svm_test = correct_rate_svm_test_Cs(bestC);

svm_weight = -svm_model.SupportVectors'*svm_model.Alpha;
svm_weight = svm_weight./svm_model.ScaleData.scaleFactor'; % MUST ADD THIS!!!!

% SVM train performance
choices_svm_train = svmclassify(svm_model,rate_lip_to_train);
correct_rate_svm_train = sum(choices_svm_train == correct_ans_train) / length(correct_ans_train);

choices_svm_all = svmclassify(svm_model,reshape(rate_lip_at_decision,N_lip,[])');   % Simply combine them (should do bootstrapping here?)
correct_rate_svm_all = sum(choices_svm_train == correct_ans_train) / length(correct_ans_train);

% Turn to the same format
choices_svm_train = reshape(choices_svm_train,n_train_each_condition,length(unique_heading),length(unique_stim_type));    
choices_svm_test = reshape(choices_svm_test,n_test_each_condition,length(unique_heading),length(unique_stim_type));    
choices_svm_all = reshape(choices_svm_all,N_rep,length(unique_heading),length(unique_stim_type));  

% --- Make final decisions here ---
% choices = choices_maxpos;
choices = choices_svm_all;

%% ======= Make choices at different TIME for calculating information ====== HH20171106
info_ts = 0:0.1:1.5;
% rate_lip_at_info_ts = squeeze(rate_lip_Ann(:,round((info_ts + stim_on_time)/dt),:,:,:)); % Or, get the population activity at the end of the trials
rate_lip_at_info_ts = squeeze(rate_lip_Ann(:,min(length(ts_Ann),round((info_ts + stim_on_time)/win_step)),:,:,:)); % Or, get the population activity at the end of the trials
info_choices_svm_all = svmclassify(svm_model,reshape(rate_lip_at_info_ts,N_lip,[])');   % Simply combine them (should do bootstrapping here?)
info_choices_svm_all = reshape(info_choices_svm_all,length(info_ts),N_rep,length(unique_heading),length(unique_stim_type));  

info = [];
for ss = 1:length(unique_stim_type)
    for itt = 1:length(info_ts)
        
        psychometric = [unique_heading' sum(reshape(info_choices_svm_all(itt,:,:,ss)==RIGHT,[],length(unique_heading)),1)'/N_rep];
        xxx = -max(unique_heading):0.1:max(unique_heading);
        
        [bias, threshold] = cum_gaussfit_max1(psychometric);
        info(itt,ss) = 2 * 1/threshold^2; % Fisher info = (d'/ds)^2 = (ds/sigma_r/ds)^2 = (1/sigma_r)^2 = (sqrt(2)/sigma_psycho)^2 = 2/threshold^2.
        
        % plot(psychometric(:,1),psychometric(:,2),'o','color',colors(unique_stim_type(ss),:),'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',10); % Psychometric
        % set(text(min(xlim),0.6+0.06*ss,sprintf('%g\n',threshold)),'color',colors(unique_stim_type(ss),:));
        % plot(xxx,normcdf(xxx,bias,threshold),'-','color',colors(unique_stim_type(ss),:),'linewid',4);
        %  axis([-8.5 8.5 0 1]);
        
    end
end

% test_svm

%% ====== SVM toy test: Generate uncorrelated population activity to test SVM decoder (after talk 20170609) ======
%{
rate_lip_toy_to_train = nan(n_train,N_lip);
rate_lip_toy_to_test = nan(n_test,N_lip);
aver_population_act = nan(N_lip,length(unique_heading),length(unique_stim_type));

% Generate toy train and test sets
n = 0;
for ss = 1:length(unique_stim_type)
    for hh = 1:length(unique_heading)
        n = n + 1;
        train_this = randperm(N_rep,n_train_each_condition);
        
        where_train_this = (n-1)*n_train_each_condition + 1 : n*n_train_each_condition;
        
        % Uncorrelated Poisson noise over the averaged population activity for this condition
        aver_this = mean(rate_lip(:,end,:,hh,ss),3)';
%         rate_lip_toy_to_train(where_train_this,:) = poissrnd(repmat(aver_this,length(where_train_this),1));
        rate_lip_toy_to_train(where_train_this,:) = (1+randn(length(where_train_this),N_lip)*4)...
                                                    .*repmat(aver_this,length(where_train_this),1);
        % Correct answer is the same as real train set
        
        aver_population_act(:,hh,ss) = aver_this;
    end
end

rate_lip_toy_to_train(rate_lip_toy_to_train<0) = 0; % Rectify

% %{
% Train linear SVM, try different Cs
% Cs = 10.^(-2.5:0.5:0);
% Cs = 1e-5;
correct_rate_toy_svm_train_Cs = nan(1,length(Cs));
svm_toy_model_Cs = [];

for ccs = 1:length(Cs)
   try
        svm_toy_model_Cs{ccs} = svmtrain(rate_lip_toy_to_train,correct_ans_train,'boxconstraint',Cs(ccs),'autoscale',1,'tol',1e-7);
        choices_svm_toy_train_Cs{ccs} = svmclassify(svm_toy_model_Cs{ccs},rate_lip_toy_to_train); % Just training set for simplification
        correct_rate_toy_svm_train_Cs(ccs) = sum(choices_svm_toy_train_Cs{ccs} == correct_ans_train) / length(correct_ans_train);
   end
end

% Select the best C
fprintf('Correct rate of (toy) SVM with different Cs: \n%s\n',num2str(correct_rate_toy_svm_train_Cs));
[~,bestToyC] = max(correct_rate_toy_svm_train_Cs);
svm_toy_model = svm_toy_model_Cs{bestToyC};

svm_toy_weight = -svm_toy_model.SupportVectors'*svm_toy_model.Alpha; 
svm_toy_weight = svm_toy_weight./svm_toy_model.ScaleData.scaleFactor'; % MUST ADD THIS!!!!

% Turn to the same format
choices_svm_toy_train = reshape(choices_svm_toy_train_Cs{bestToyC},n_train_each_condition,length(unique_heading),length(unique_stim_type));    
%}

%% ====== Fig.0.5 Compare MaxPos and SVM (training, test) ======
% %{
% Compare different decoding methods (MaxPos, SVM_all, SVM_test)
if analysis_switch(2)
    figure(2135); clf;
    set(gcf,'name','Decoding methods');
    set(gcf,'uni','norm','pos',[0.013        0.06       0.972       0.464]);
    
    decoding_methods = {choices_maxpos, 'MaxPos';
        choices_svm_all, 'SVM all';
        choices_svm_train, 'SVM train';
        choices_svm_test, 'SVM test';
        % choices_svm_toy_train, 'SVM toy train';
        };
    hs = tight_subplot(2,size(decoding_methods,1),0.05,0.1);
    
    % ======= Behavior ========
    
    % --- Bootstrapping ---
    if ION_cluster
        N_bootstrapping = 200;
        N_reps_boot = N_rep;
    else
        N_bootstrapping = 20;
        N_reps_boot = 20;
    end
    
    for dm = 1:size(decoding_methods,1)
        thres_boot = nan(length(unique_stim_type),N_bootstrapping);
        
        for ss = 1:length(unique_stim_type)
            
            psychometric = [unique_heading' sum(reshape(decoding_methods{dm,1}(:,:,ss)==RIGHT,[],length(unique_heading)),1)'/size(decoding_methods{dm,1},1)];
            
            axes(hs(dm*2-1)); set(gca,'color','none');
            
            parfor bb = 1:N_bootstrapping
                psycho_boot = sum(bsxfun(@le,rand(size(psychometric,1),N_reps_boot),psychometric(:,2)),2)/N_reps_boot;
                [~, thres_boot_this] = cum_gaussfit_max1([unique_heading' psycho_boot]);
                thres_boot(ss,bb) = thres_boot_this;
            end
            
            set(bar(ss,mean(thres_boot(ss,:))),'facecolor',colors(unique_stim_type(ss),:)); hold on;
            errorbar(ss,mean(thres_boot(ss,:)),std(thres_boot(ss,:)),'color',colors(unique_stim_type(ss),:));
            
            axes(hs(dm*2));
            
            plot(psychometric(:,1),psychometric(:,2),'o','color',colors(unique_stim_type(ss),:),'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',10); % Psychometric
            hold on;
            xxx = -max(unique_heading):0.1:max(unique_heading);
            
            [bias, threshold] = cum_gaussfit_max1(psychometric);
            psycho(ss,:) = [bias,threshold];
            
            set(text(min(xlim),0.6+0.06*ss,sprintf('%g\n',threshold)),'color',colors(unique_stim_type(ss),:));
            
            plot(xxx,normcdf(xxx,bias,threshold),'-','color',colors(unique_stim_type(ss),:),'linewid',4);
            axis([-8.5 8.5 0 1]);
            
        end
        
        axes(hs(dm*2-1));
        thres_boot(4,:) = (thres_boot(1,:).^(-2) + thres_boot(2,:).^(-2)).^(-1/2);
        set(bar(4,mean(thres_boot(4,:))),'facecolor','k'); hold on;
        errorbar(4,mean(thres_boot(4,:)),std(thres_boot(4,:)),'k');
        
        title(sprintf('%s, N_{reps} = %g',decoding_methods{dm,2},size(decoding_methods{dm,1},1)));
        
        % Pred
        axes(hs(dm*2));
        if length(unique_stim_type) == 3
            pred_thres = (psycho(1,2)^(-2)+psycho(2,2)^(-2))^(-1/2);
            set(text(min(xlim),0.6+0.06*(ss+1),sprintf('pred = %g\n',pred_thres)),'color','k');
        end
        
    end
    
    if ION_cluster
        fprintf('...Saving fig 0p5...');
        file_name = sprintf('./result/%s0p5_DecodingMethods%s',save_folder,para_override_txt);
        export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end
    
    disp('Decoding done!');
end
%}

%% ====== Fig.1 Overview ======
figure(1000); clf;
set(gcf,'name','Overview');
set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);

hpsy = axes('position',[ 0.0540    0.6860    0.0890      0.1360]);

hs = tight_subplot(3,4,0.05);

% ======= Behavior ========
% Psychometric

% --- Bootstrapping ---
if ION_cluster
    N_bootstrapping = 1000;
    N_reps_boot = N_rep;
else
    N_bootstrapping = 20;
    N_reps_boot = 20; 
end
thres_boot = nan(length(unique_stim_type),N_bootstrapping);

for ss = 1:length(unique_stim_type)
    
    %     % I just flip the psychometric curve to the negative headings
%     psychometric = [unique_heading' sum(reshape(choices{ss},[],length(unique_heading)),1)'/N_rep];
%     fake_psy = flipud(psychometric(unique_heading>0,:));
%     fake_psy(:,1) = -fake_psy(:,1);
%     fake_psy(:,2) = 1-fake_psy(:,2);
%     psychometric = [fake_psy ; psychometric];

    psychometric = [unique_heading' sum(reshape(choices(:,:,ss)==RIGHT,[],length(unique_heading)),1)'/N_rep];

    axes(hs(1)); set(gca,'color','none');
    
    
    parfor bb = 1:N_bootstrapping
        psycho_boot = sum(bsxfun(@le,rand(size(psychometric,1),N_reps_boot),psychometric(:,2)),2)/N_reps_boot;
        [~, thres_boot_this] = cum_gaussfit_max1([unique_heading' psycho_boot]);
        thres_boot(ss,bb) = thres_boot_this;
    end
    
    set(bar(ss,mean(thres_boot(ss,:))),'facecolor',colors(unique_stim_type(ss),:)); hold on;
    errorbar(ss,mean(thres_boot(ss,:)),std(thres_boot(ss,:)),'color',colors(unique_stim_type(ss),:));
    
    axes(hpsy);
    plot(psychometric(:,1),psychometric(:,2),'o','color',colors(unique_stim_type(ss),:),'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',10); % Psychometric
    hold on;
    xxx = -max(unique_heading):0.1:max(unique_heading);
    
    [bias, threshold] = cum_gaussfit_max1(psychometric);
    psycho(ss,:) = [bias,threshold];
    
    set(text(min(xlim),0.6+0.06*ss,sprintf('%g\n',threshold)),'color',colors(unique_stim_type(ss),:));
    
    plot(xxx,normcdf(xxx,bias,threshold),'-','color',colors(unique_stim_type(ss),:),'linewid',4);
    axis([-8.5 8.5 0 1]);
    
    % Chronometric
    axes(hs(2));
    for hh = 1:length(unique_heading)
        select_correct = (sign(choices(:,hh,ss)) == sign(unique_heading(hh))) | unique_heading(hh)== 0 ; % Only correct trials should be plotted in RT
        
        RT_this_cc_mean(hh,1) = mean(RT(select_correct,hh,ss),1);
        RT_this_cc_se(hh,1) = std(RT(select_correct,hh,ss),[],1);
        
%         RT_this_cc_mean_wrong(hh,1) = mean(RT(~select_correct,hh,ss),1);
%         RT_this_cc_se_wrong(hh,1) = std(RT(~select_correct,hh,ss),[],1);
        
    end
    
    errorbar(unique_heading+(ss-2)*0.2,RT_this_cc_mean,RT_this_cc_se,'o-','color',colors(unique_stim_type(ss),:),...
        'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',10,'linew',2);
    hold on;
%     errorbar(unique_heading+(ss-2)*0.2,RT_this_cc_mean_wrong,RT_this_cc_se_wrong,'o--','color',colors(unique_stim_type(ss),:),...
%         'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',7,'linew',1);
    axis([-8.5 8.5 0 1.5]);

end

axes(hs(2));
plot(vel/max(vel)*max(xlim)/3-8,t_motion,'k--','linew',2);
if if_bound_RT
    title(sprintf('RT, correct trials, Bound = %s',num2str(decis_bound)));
else
    title(sprintf('Fixed duration'));
end
axes(hs(1));
thres_boot(4,:) = (thres_boot(1,:).^(-2) + thres_boot(2,:).^(-2)).^(-1/2);
set(bar(4,mean(thres_boot(4,:))),'facecolor','k'); hold on;
errorbar(4,mean(thres_boot(4,:)),std(thres_boot(4,:)),'k');

title(sprintf('N_{reps} = %g, Readout at RT = %g',N_rep,read_out_at_the_RT));

% Pred
axes(hpsy);
if length(unique_stim_type) == 3
    pred_thres = (psycho(1,2)^(-2)+psycho(2,2)^(-2))^(-1/2);
    set(text(min(xlim),0.6+0.06*(ss+1),sprintf('pred = %g\n',pred_thres)),'color','k');
end


% --- Calculate ps ---
xs = thres_boot';
ps_psycho = nan(size(xs,2));

% NOTE: ttest or ttest2 cannot be used because they're not "sampling data"!!
% Instead, use the definition of p-values (and different subsamplings can be pooled)
for g1 = 1:size(xs,2)
    for g2 = g1 + 1:size(xs,2)
        p = 2 * sum(xs(:,g1) > xs(:,g2))/size(xs,1);
        if p > 1, p = 2-p; end
        ps_psycho(g1,g2) = p;
        ps_psycho(g2,g1) = p;
    end
end

display(ps_psycho);

ps_psycho_ttest = nan(size(xs,2));


%}

% ===== Averaged PSTH, different angles =====
%%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
half_range = 30;
to_calculate_PSTH_cells_ind = find(abs(abs(prefs_lip)-90) <= half_range);   
% to_calculate_PSTH_cells_ind = find(prefs_lip>-117,1);
N_sample_cell = length(to_calculate_PSTH_cells_ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No more down sample of cells since cells are not enough already...
% if N_sample_cell > 60
%     to_calculate_PSTH_cells_ind = to_calculate_PSTH_cells_ind(round(1:N_sample_cell/30:N_sample_cell));
%     N_sample_cell = length(to_calculate_PSTH_cells_ind);
% end

% to_plot_grouped_PSTH = find(to_calculate_PSTH_cells_ind == right_targ_ind);

to_flip = prefs_lip(to_calculate_PSTH_cells_ind) < 0;

PSTH_correct_mean_headings = nan(N_sample_cell,length(ts_Ann),length(unique_stim_type),length(unique_abs_heading),2);
PSTH_correct_mean_headings_norm = nan(N_sample_cell,length(ts_Ann),length(unique_stim_type),length(unique_abs_heading),2);
PSTH_correct_mean_allheading = nan(N_sample_cell,length(ts_Ann),length(unique_stim_type),2);
PSTH_correct_mean_allheading_norm = nan(N_sample_cell,length(ts_Ann),length(unique_stim_type),2);

% For linear regression
PSTH_all_mean_headings = nan(N_sample_cell,length(ts_Ann),length(unique_stim_type),length(unique_abs_heading),2);

pref_null = [RIGHT LEFT];

% === Calculate cells' dynamic range for normalization HH20170714 ===
norm_time_range = 0 < ts_Ann;
dynamic_max = max(reshape(rate_lip_Ann(to_calculate_PSTH_cells_ind,norm_time_range,:,:,:),N_sample_cell,[]),[],2);
dynamic_min = min(reshape(rate_lip_Ann(to_calculate_PSTH_cells_ind,norm_time_range,:,:,:),N_sample_cell,[]),[],2);
% figure(); plot(dynamic_max); hold on; plot(dynamic_min)

% === Group PSTH data and Norm_PSTH data ===

for ss = 1:length(unique_stim_type)
    
%     zu_all{ss} = []; % Cache for VarCE
    
    for cc = 1:2 % Pref and Null
        PSTH_cache = [];
        PSTH_cache_norm = [];
         
        for abshh = 1:length(unique_abs_heading)         % For different abs(headings)
            % Assuming PREF = RIGHT for all the cells
            this_heading_ind = find(unique_heading == unique_abs_heading(abshh) * pref_null(cc),1);
            this_correct_ind = find(choices(:,this_heading_ind,ss) == pref_null(cc)); % Only correct trials
            
            % PSTH separated for each heading
            PSTH_correct_raw_this = rate_lip_Ann(to_calculate_PSTH_cells_ind,:,this_correct_ind,this_heading_ind,ss);
            PSTH_correct_mean_headings(:,:,ss,abshh,cc) = mean(PSTH_correct_raw_this,3);

            PSTH_all_raw_this = rate_lip_Ann(to_calculate_PSTH_cells_ind,:,:,this_heading_ind,ss);            
            PSTH_all_mean_headings(:,:,ss,abshh,cc) = mean(PSTH_all_raw_this,3); % Note that there's a duplication in zero heading
            
            PSTH_correct_raw_this_norm = bsxfun(@rdivide,bsxfun(@minus,PSTH_correct_raw_this,dynamic_min),(dynamic_max-dynamic_min));
            PSTH_correct_mean_headings_norm(:,:,ss,abshh,cc) = mean(PSTH_correct_raw_this_norm,3);
            
            % PSTH grouped cache for all headings
            PSTH_cache = cat(3,PSTH_cache,PSTH_correct_raw_this);
            PSTH_cache_norm = cat(3,PSTH_cache_norm,PSTH_correct_raw_this_norm);

            % For VarCE: all choice
%             if ~(unique_abs_heading(abshh)==0 && cc==2) % Skip abs(heading) = 0 and cc = 2 because cc = 1 already includes all trials for 0 heading
%                 PSTH_all_raw_this_forVarCE = rate_lip(to_calculate_PSTH_cells_ind,:,:,this_heading_ind,ss);
%                 zu_this = bsxfun(@minus, PSTH_all_raw_this_forVarCE, mean(PSTH_all_raw_this_forVarCE,3));
%                 zu_all{ss} = cat(3,zu_all{ss},zu_this);
%             end
        end
        
        % PSTH for all headings
        PSTH_correct_mean_allheading(:,:,ss,cc) = mean(PSTH_cache,3);
        PSTH_correct_mean_allheading_norm(:,:,ss,cc) = mean(PSTH_cache_norm,3);
    
    end
end
% Flip pref and null if PREF = LEFT for some of the cells
v = version;

if str2num(v(findstr(v,'R')+1:findstr(v,'R')+4)) >= 2014   % Matlab version issue
    PSTH_correct_mean_headings(to_flip,:,:,:,:) = flip(PSTH_correct_mean_headings(to_flip,:,:,:,:),5);
    PSTH_correct_mean_allheading(to_flip,:,:,:) = flip(PSTH_correct_mean_allheading(to_flip,:,:,:),4);
    PSTH_correct_mean_headings_norm(to_flip,:,:,:,:) = flip(PSTH_correct_mean_headings_norm(to_flip,:,:,:,:),5);
    PSTH_correct_mean_allheading_norm(to_flip,:,:,:) = flip(PSTH_correct_mean_allheading_norm(to_flip,:,:,:),4);
    PSTH_all_mean_headings(to_flip,:,:,:,:) = flip(PSTH_all_mean_headings(to_flip,:,:,:,:),5);
else
    PSTH_correct_mean_headings(to_flip,:,:,:,:) = flipdim(PSTH_correct_mean_headings(to_flip,:,:,:,:),5);
    PSTH_correct_mean_allheading(to_flip,:,:,:) = flipdim(PSTH_correct_mean_allheading(to_flip,:,:,:),4);    
    PSTH_correct_mean_headings_norm(to_flip,:,:,:,:) = flipdim(PSTH_correct_mean_headings_norm(to_flip,:,:,:,:),5);
    PSTH_correct_mean_allheading_norm(to_flip,:,:,:) = flipdim(PSTH_correct_mean_allheading_norm(to_flip,:,:,:),4);    
    PSTH_all_mean_headings(to_flip,:,:,:,:) = flipdim(PSTH_all_mean_headings(to_flip,:,:,:,:),5);
end

diff_PSTH_correct_mean_headings = - diff(PSTH_correct_mean_headings,[],5);
diff_PSTH_correct_mean_allheading = - diff(PSTH_correct_mean_allheading,[],4);
diff_PSTH_correct_mean_headings_norm = - diff(PSTH_correct_mean_headings_norm,[],5);
diff_PSTH_correct_mean_allheading_norm = - diff(PSTH_correct_mean_allheading_norm,[],4);

% -- Plotting --
y_max = -inf; y_min = inf;
colors_hue = [240 0 120]/360;

for ss = 1:length(unique_stim_type)
    yyy = squeeze(reshape(mean(PSTH_correct_mean_headings(:,:,ss,:,:),1),1,[],1,length(unique_abs_heading)*2));
    yyy = shiftdim(yyy,-1);
    
    colors_angles_hsv(:,2) = linspace(0.2,1,length(unique_abs_heading));
    colors_angles_hsv(:,1) = colors_hue(ss);
    colors_angles_hsv(:,3) = 0.9;
    colors_angles = hsv2rgb(colors_angles_hsv);
    colors_angles = mat2cell(colors_angles,ones(length(unique_abs_heading),1));
    
    h = SeriesComparison(yyy,ts_Ann,...
        'Colors',colors_angles,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
        'ErrorBar',0,'Xlabel',[],'axes',hs(3+ss));
    legend off;
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');

    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    xlim([min(ts_Ann),max(ts_Ann)]);
end

% title(hs(4),sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(to_plot_grouped_PSTH))));
title(hs(4),sprintf('Grouped %g cells',N_sample_cell));

set(hs(4:6),'ylim',[y_min y_max]);

% =====  Delta-PSTH (stim type), different heading, correct only ====
y_max = -inf; y_min = inf;

for abshh = 1:length(unique_abs_heading)
    
    axes(hs(6+abshh)); hold on;
    
    for ss = 1:length(unique_stim_type)
        plot(ts_Ann,mean(diff_PSTH_correct_mean_headings(:,:,ss,abshh),1),'color',colors(unique_stim_type(ss),:),'linew',2.5);
    end
    if length(unique_stim_type) == 3
        plot(ts_Ann,mean(sum(diff_PSTH_correct_mean_headings(:,:,[1 2],abshh),3),1),'k--');
    end
    title(sprintf('|heading| = %g',unique_abs_heading(abshh)));
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
end

% ======= Delta-PSTH, all heading, correct only =====
axes(hs(12)); hold on;    title('All headings');

mean_diff_PSTH_correct_allheading = mean(diff_PSTH_correct_mean_allheading,1); % Save for fitting real data. HH20170808

yyy = reshape(diff_PSTH_correct_mean_allheading(:,:,:),N_sample_cell,length(ts_Ann),[]);
try
h = SeriesComparison(yyy,ts_Ann,...
    'Colors',colors,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
    'ErrorBar',2,'Xlabel',[],'axes',hs(12));
catch
    keyboard
end
legend off;
plot(hs(12),t_motion,vel/max(vel)*max(ylim)/5,'k--'); axis tight;
title(hs(12),'raw');

% for ss = 1:length(unique_stim_type)
%     mean_diff_PSTH_correct_allheading(:,ss) = mean(diff_PSTH_correct_mean_allheading(:,:,ss),1);
%     plot(ts,mean_diff_PSTH_correct_allheading(:,ss),'color',colors(unique_stim_type(ss),:),'linew',2.5);
%     
% %     errorbar(ts,mean(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,ss),1),...
% %         std(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,ss),[],1)/sqrt(length(to_plot_grouped_PSTH)),'color',colors(unique_stim_type(ss),:),'linew',2.5);
% end

if length(unique_stim_type) == 3
    plot(ts_Ann,mean(sum(diff_PSTH_correct_mean_allheading(:,:,[1 2]),3),1),'k--');
end

axis tight;
y_min = min(y_min,min(ylim));
y_max = max(y_max,max(ylim));

set(hs(7:12),'ylim',[y_min y_max]);

% % ====== VarCE ======
% Info (1/variance^2)
axes(hs(3)); hold on; 
info(info>=90)=0; % Some weird values due to bad fitting
for ss = 1:length(unique_stim_type)
    plot(info_ts,info(:,ss),'o-','color',colors(unique_stim_type(ss),:),'markerfacecolor',colors(unique_stim_type(ss),:));
end

plot(info_ts,sum(info(:,[1 2]),2),'k--');

% title('VarCE');
% for ss = 1:length(unique_stim_type)
%     plot(ts,mean(var(zu_all{ss}(:,:,:)*0.06,[],3),1),...  % 60 ms spike counting window in Ann's paper
%         'color',colors(unique_stim_type(ss),:),'linew',2.5);
% end
% plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
% xlim([min(ts),max(ts)]);


if ION_cluster
    export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s1_Overview%s.png',save_folder,para_override_txt));
    saveas(gcf,sprintf('./result/%s1_Overview%s.fig',save_folder,para_override_txt),'fig');
    h_grouped = [h_grouped hs(1) hpsy hs(6) hs(12) hs(3)]; % Add h_grouped
end
disp('Overview done');

if analysis_switch(3)
%% ====== Fig.1.5 Overview_normalized ======
    figure(7141706); clf;
    set(gcf,'name','Overview_normalized');
    set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);
    
    hs = tight_subplot(3,4,0.05);
    
    % -- Plotting PSTH / norm PSTH, all in one --
    yyy = reshape(PSTH_correct_mean_allheading,N_sample_cell,length(ts_Ann),[]);
    h = SeriesComparison(yyy,ts_Ann,...
        'Colors',colors,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
        'ErrorBar',2,'Xlabel',[],'axes',hs(1));
    legend off;
    plot(hs(1),t_motion,vel/max(vel)*max(ylim)/5,'k--'); axis tight;
    title(hs(1),'raw');
    
    yyy = reshape(PSTH_correct_mean_allheading_norm,N_sample_cell,length(ts_Ann),[]);
    h = SeriesComparison(yyy,ts_Ann,...
        'Colors',colors,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
        'ErrorBar',2,'Xlabel',[],'axes',hs(2));
    legend off;
    plot(hs(2),t_motion,vel/max(vel)*max(ylim)/5,'k--'); axis tight;
    title(hs(2),'norm');
    
    % -- Plotting norm PSTHs, grouped by stim type, different angles --
    y_max = -inf; y_min = inf;
    colors_hue = [240 0 120]/360;
    
    for ss = 1:length(unique_stim_type)
        yyy = squeeze(reshape(mean(PSTH_correct_mean_headings_norm(:,:,ss,:,:),1),1,[],1,length(unique_abs_heading)*2));
        yyy = shiftdim(yyy,-1);
        
        colors_angles_hsv(:,2) = linspace(0.2,1,length(unique_abs_heading));
        colors_angles_hsv(:,1) = colors_hue(ss);
        colors_angles_hsv(:,3) = 0.9;
        colors_angles = hsv2rgb(colors_angles_hsv);
        colors_angles = mat2cell(colors_angles,ones(length(unique_abs_heading),1));
        
        h = SeriesComparison(yyy,ts_Ann,...
            'Colors',colors_angles,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
            'ErrorBar',0,'Xlabel',[],'axes',hs(3+ss));
        legend off;
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
        axis tight;
        y_min = min(y_min,min(ylim));
        y_max = max(y_max,max(ylim));
        xlim([min(ts_Ann),max(ts_Ann)]);
    end
    
    % title(hs(4),sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(to_plot_grouped_PSTH))));
    title(hs(4),sprintf('Grouped %g cells',N_sample_cell));
    
    set(hs(4:6),'ylim',[y_min y_max]);
    
    % =====  Delta-PSTH (stim type), different heading, correct only ====
    y_max = -inf; y_min = inf;
    
    for abshh = 1:length(unique_abs_heading)
        
        axes(hs(6+abshh)); hold on;
        
        for ss = 1:length(unique_stim_type)
            plot(ts_Ann,mean(diff_PSTH_correct_mean_headings_norm(:,:,ss,abshh),1),'color',colors(unique_stim_type(ss),:),'linew',2.5);
        end
        if length(unique_stim_type) == 3
            plot(ts_Ann,mean(sum(diff_PSTH_correct_mean_headings_norm(:,:,[1 2],abshh),3),1),'k--');
        end
        title(sprintf('|heading| = %g',unique_abs_heading(abshh)));
        
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        axis tight;
        y_min = min(y_min,min(ylim));
        y_max = max(y_max,max(ylim));
    end
    
    % ======= Delta-PSTH, all heading, correct only =====
    axes(hs(12)); hold on;    title('All headings');
    
    for ss = 1:length(unique_stim_type)
        plot(ts_Ann,mean(diff_PSTH_correct_mean_allheading_norm(:,:,ss),1),'color',colors(unique_stim_type(ss),:),'linew',2.5);
        %     errorbar(ts,mean(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,ss),1),...
        %         std(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,ss),[],1)/sqrt(length(to_plot_grouped_PSTH)),'color',colors(unique_stim_type(ss),:),'linew',2.5);
    end
    if length(unique_stim_type) == 3
        plot(ts_Ann,mean(sum(diff_PSTH_correct_mean_allheading_norm(:,:,[1 2]),3),1),'k--');
    end
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    
    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    
    set(hs(7:12),'ylim',[y_min y_max]);
    
    
    if ION_cluster
        export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s1p5_Overview_norm%s.png',save_folder,para_override_txt));
        saveas(gcf,sprintf('./result/%s1p5_Overview_norm%s.fig',save_folder,para_override_txt),'fig');
        h_grouped = [h_grouped hs(1) hpsy hs(6)]; % Add h_grouped
    end
    disp('Overview done');
end

%% ====== Fig.2 Example Pref and Null traces (correct only) ======
%%{
if analysis_switch(4)
    to_plot_abs_headings = [0 4 8];
    % to_plot_cell_ind = right_targ_ind;  % Cell-based plotting
    [~, ind] = max(mean(mean(diff_PSTH_correct_mean_allheading(:,:,:),2),3),[],1);
    to_plot_cell_ind = to_calculate_PSTH_cells_ind(ind);
    n_to_plot_trace = 5;
    
    set(figure(1002),'name','Example PSTHs (correct only)'); clf;
    set(gcf,'uni','norm','pos',[0.005       0.056       0.33*length(to_plot_abs_headings)       0.832]);
    hs = tight_subplot(3,2*length(to_plot_abs_headings),[0.04 0.03]);
    
    pref_null = [RIGHT LEFT]; % Default: PREF = RIGHT, NULL = LEFT
    if prefs_lip(to_plot_cell_ind)<0  % If this target cell has a negative pref heading
        pref_null = fliplr(pref_null);
    end
    
    y_max = -inf; y_min = inf; to_sync = [];
    
    for tph = 1:length(to_plot_abs_headings)
        
        for ss = 1:length(unique_stim_type)
            
            for cc = 1:2 % Pref and Null
                
                this_heading_ind = find(unique_heading == to_plot_abs_headings(tph) * pref_null(cc),1);
                % this_correct_ind = find(choices(:,this_heading_ind,ss) == pref_null(cc)); % Only correct trials
                this_correct_ind = find(choices(:,this_heading_ind,ss) > -inf); % All trials
                
                for tt = 1:ceil(length(this_correct_ind)/n_to_plot_trace):length(this_correct_ind)
                    
                    RT_this = RT(this_correct_ind(tt),this_heading_ind,ss);
                    rate_at_RT = rate_lip_Ann(to_plot_cell_ind,find(ts_Ann<=RT_this,1,'last'),this_correct_ind(tt),this_heading_ind,ss);
                    
                    % --- LIP ---
                    axes(hs(6*(tph-1)+ss));
                    
                    if cc == 1 % Pref
                        plot(ts_Ann,rate_lip_Ann(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'color',colors(unique_stim_type(ss),:),'linewid',2); hold on;
                        plot([RT_this RT_this],[rate_at_RT+5 rate_at_RT-5],'m','linew',5);
                    else
                        plot(ts_Ann,rate_lip_Ann(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'k--','linewid',1);
                        plot([RT_this RT_this],[rate_at_RT+5 rate_at_RT-5],'k','linew',5);
                    end
                    
                    axis tight; y_min = min(y_min,min(ylim)); y_max = max(y_max,max(ylim)); to_sync = [to_sync gca];
                    
                    % --- Int ---
                    axes(hs(6*(tph-1)+ss+3));
                    if cc == 1 % Pref
                        plot(ts_Ann,rate_int_Ann(round(to_plot_cell_ind/N_lip*N_int),:,this_correct_ind(tt),this_heading_ind,ss),'color',colors(unique_stim_type(ss),:),'linewid',2); hold on;
                    else
                        plot(ts_Ann,rate_int_Ann(round(to_plot_cell_ind/N_lip*N_int),:,this_correct_ind(tt),this_heading_ind,ss),'k--','linewid',1); hold on;
                    end
                    %     if if_bound_RT
                    %         plot(xlim,[decis_bound(unique_stim_type(ss)) decis_bound(unique_stim_type(ss))],'c--');
                    % %         ylim([min(ylim),decis_bound*1.1]);
                    %     end
                    
                end
            end
            
            axes(hs(6*(tph-1)+ss));
            title(sprintf('rate\\_lip, pref = %g, |heading| = %g',prefs_lip(to_plot_cell_ind), to_plot_abs_headings(tph)));
            plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
            axis tight;
            if if_bound_RT
                plot(xlim,[decis_bound(unique_stim_type(ss)) decis_bound(unique_stim_type(ss))],'c--');
                %         ylim([min(ylim),decis_bound(k)*1.1]);
            end
            
            axes(hs(6*(tph-1)+ss+3));
            title(sprintf('rate\\_int'));
            plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
            axis tight;
            
        end
        
    end
    
    set(to_sync,'ylim',[y_min y_max]);
    
    if ION_cluster
        export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s2_Example%s.png',save_folder,para_override_txt));
        saveas(gcf,sprintf('./result/%s2_Example%s.fig',save_folder,para_override_txt),'fig');
    end
    
    disp('Example done');
end
%}

%% ====== Fig.3 Different cells, delta PSTH ======
if analysis_switch(5)
    figure(1002); clf;
    set(gcf,'name','Different cells, delta PSTH');
    set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);
    
    % -- diff_PSTH for all headings --
    % Only plot diff_PSTH and raw PSTH for part of cells (Max number = 30)
    to_plot_single_cell_PSTH = unique(round(1:N_sample_cell/30:N_sample_cell));
    
    m = fix(sqrt(length(to_plot_single_cell_PSTH)+5));
    n = ceil((length(to_plot_single_cell_PSTH)+5)/m);
    [~,hs] = tight_subplot(m,n,0.05);
    y_max = -inf; y_min = inf;
    
    for cc = 1:length(to_plot_single_cell_PSTH)
        axes(hs(cc)); hold on;
        for ss = 1:length(unique_stim_type)
            plot(ts_Ann,diff_PSTH_correct_mean_allheading(to_plot_single_cell_PSTH(cc),:,ss),'color',colors(unique_stim_type(ss),:),'linew',2.5);
        end
        if length(unique_stim_type) == 3
            plot(ts_Ann,sum(diff_PSTH_correct_mean_allheading(to_plot_single_cell_PSTH(cc),:,[1 2]),3),'k--');
        end
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
        
        title(sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc)))));
        
        if abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc))-left_targ_ind) == min(abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH)-left_targ_ind)) ...
                || abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc)) - right_targ_ind) == min(abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH) - right_targ_ind))
            set(gca,'color','y');
        end
        
        axis tight;
        y_min = min(y_min,min(ylim));
        y_max = max(y_max,max(ylim));
    end
    
    set(hs,'ylim',[y_min y_max]);
    
    % -- Weights --
    axes(hs(end-4));
    h1 = imagesc(prefs_lip,prefs_lip,w_lip_lip'); axis  xy; title('LIP->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{lip}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    % set(h1,'ButtonDownFcn',{@plot_weight,prefs_lip,prefs_lip,w_lip_lip'});
    
    axes(hs(end-3));
    h2 = imagesc(prefs_int,prefs_lip,w_lip_int');  axis xy; title('Int->LIP');
    xlabel('\theta_{lip}'); ylabel('\theta_{int}');
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    % set(h2,'ButtonDownFcn',{@plot_weight,prefs_lip,prefs_int,w_lip_int'});
    
    axes(hs(end-2));
    h3 = imagesc(prefs_vis,prefs_int,w_int_vis'); axis xy; title('Vis->Int');
    xlabel('\theta_{int}'); ylim([-30 30]);    ylabel('\theta_{vis}');
    set(gca,'xtick',-180:90:180);
    % set(h3,'ButtonDownFcn',{@plot_weight,prefs_int,prefs_vis,w_int_vis'});
    
    axes(hs(end-1));
    h3 = imagesc(prefs_vest,prefs_int,w_int_vest'); axis xy; title('Vest->Int');
    xlabel('\theta_{int}'); ylim([-30 30]);    ylabel('\theta_{vest}');
    set(gca,'xtick',-180:90:180);
    % set(h3,'ButtonDownFcn',{@plot_weight,prefs_int,prefs_vis,w_int_vis'});
    
    colormap hot;
    
    if ION_cluster
        export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s3_Cells%s.png',save_folder,para_override_txt));
        saveas(gcf,sprintf('./result/%s3_Cells%s.fig',save_folder,para_override_txt),'fig');
        h_grouped = [h_grouped]; % Add h_grouped
    end
    disp('Different cells done');
end

%% ====== Fig.3.5 Different cells, PSTH ======
if analysis_switch(6)
    
    figure(1003); clf;
    set(gcf,'name','Different cells, PSTH');
    set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);
    
    % -- PSTH for all headings --
    
    m = fix(sqrt(length(to_plot_single_cell_PSTH)+5));
    n = ceil((length(to_plot_single_cell_PSTH)+5)/m);
    [~,hs] = tight_subplot(m,n,0.05);
    y_max = -inf; y_min = inf;
    
    for cc = 1:length(to_plot_single_cell_PSTH)
        axes(hs(cc)); hold on;
        
        for ss = 1:length(unique_stim_type)
            plot(ts_Ann,PSTH_correct_mean_allheading(to_plot_single_cell_PSTH(cc),:,ss,1),'linestyle','-','color',colors(unique_stim_type(ss),:),'linew',2.5);
            plot(ts_Ann,PSTH_correct_mean_allheading(to_plot_single_cell_PSTH(cc),:,ss,2),'linestyle','--','color',colors(unique_stim_type(ss),:),'linew',2.5);
        end
        
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
        title(sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc)))));
        
        % Mark the nearest cells to the targets
        if abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc))-left_targ_ind) == min(abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH)-left_targ_ind)) ...
                || abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH(cc)) - right_targ_ind) == min(abs(to_calculate_PSTH_cells_ind(to_plot_single_cell_PSTH) - right_targ_ind))
            set(gca,'color','y');
        end
        
        axis tight;
        y_min = min(y_min,min(ylim));
        y_max = max(y_max,max(ylim));
    end
    
    set(hs,'ylim',[y_min y_max]);
    
    if ION_cluster
        export_fig('-painters','-nocrop','-m1.5' ,sprintf('./result/%s3p5_Cells_PSTH%s.png',save_folder,para_override_txt));
        saveas(gcf,sprintf('./result/%s3p5_Cells_PSTH%s.fig',save_folder,para_override_txt),'fig');
    end
    disp('Different cells done');
end

%% ====== Fig.4 Heterogeneity and Decoding Weights ========
if analysis_switch(7)
    
    if length(unique_stim_type) == 3
        
        figure(1418); clf
        set(gcf,'name','Heterogeneity and Decoding Weights');
        set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);
        
        hs = tight_subplot(2,3,0.1);
        
        % -- Multisensory enhancement --
        axes(hs(1));
        div_vest = mean(diff_PSTH_correct_mean_allheading(:,:,1),2);
        div_vis = mean(diff_PSTH_correct_mean_allheading(:,:,2),2);
        div_comb = mean(diff_PSTH_correct_mean_allheading(:,:,3),2);
        div_comb_minus_vest = div_comb - div_vest;
        div_comb_minus_vis = div_comb - div_vis;
        
        cols = colormap(jet);
        
        % Fit color to to N_sample_cell/4
        cols = interp1(1:length(cols),cols,1:length(cols)/(N_sample_cell/4+4):length(cols));
        
        for cc = 1:N_sample_cell
            dis_prop = min(abs(cc-N_sample_cell/4),abs(cc-N_sample_cell/4*3));
            this_col = cols(1+round(dis_prop),:);
            plot(div_comb_minus_vest(cc),div_comb_minus_vis(cc),'og','markersize',7,'color',this_col,'linewid',2); hold on;
        end
        
        xlabel('Comb - Vest'); ylabel('Comb - Vis');
        max_axis = max(abs([ylim xlim]));
        plot([0 0],[-max_axis max_axis],'--k');
        plot([-max_axis max_axis],[0 0],'--k');
        plot([-max_axis max_axis],[-max_axis max_axis],'--k');
        axis tight square;
        
        % -- SVM weights --
        axes(hs(2));
        plot(prefs_lip,svm_weight,'k');
        xlim([-180 180]);
        set(gca,'xtick',-180:90:180);
        xx = prefs_lip(to_calculate_PSTH_cells_ind);
        yy = svm_weight(to_calculate_PSTH_cells_ind);
        hold on;
        for cc = 1:N_sample_cell
            dis_prop = min(abs(cc-N_sample_cell/4),abs(cc-N_sample_cell/4*3));
            this_col = cols(1+round(dis_prop),:);
            plot(xx(cc),yy(cc),'og','markersize',7,'color',this_col,'linewid',2); hold on;
        end
        
        % plot(prefs_lip,svm_toy_weight,'m','linew',2);
        
        xlabel('Prefs_{LIP}'); ylabel('SVM weight');
        title(sprintf('Best C = %g, CR = %g',Cs(bestC),correct_rate_svm_test_Cs(bestC)));
        
        % -- SVM weights vs Div_vest --
        xxs = {div_vest, 'Divergence vest';
            div_vis,  'Divergence vis';
            div_comb, 'Divergence comb';
            div_comb_minus_vest + div_comb_minus_vis, '(comb-vest)+(comb-vis)'};
        
        for xxxx = 1:length(xxs)
            xx = xxs{xxxx,1};
            yy = abs(svm_weight(to_calculate_PSTH_cells_ind));
            hl = LinearCorrelation(xx,yy,'Axes',hs(2+xxxx));
            hold on;
            delete([hl.leg hl.group.dots]);
            for cc = 1:N_sample_cell
                dis_prop = min(abs(cc-N_sample_cell/4),abs(cc-N_sample_cell/4*3));
                this_col = cols(1+round(dis_prop),:);
                plot(xx(cc),yy(cc),'og','markersize',7,'color',this_col,'linewid',2); hold on;
            end
            text(min(xlim),min(ylim)+range(ylim)*0.2,sprintf('r^2=%g\n p=%g',hl.group.r_square,hl.group.p))
            xlabel(xxs{xxxx,2}); ylabel('abs(svm weight)');
        end
        
        if ION_cluster
            file_name = sprintf('./result/%s4_Heterogeneity%s',save_folder,para_override_txt);
            export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
            saveas(gcf,[file_name '.fig'],'fig');
            h_grouped = [h_grouped hs(1)]; % Add h_grouped
        end
        
    end
end

%% ===== Fig.5 Linear regression of Rcomb = Wvest * Rvest + Wvis * Rvis
if analysis_switch(8)
    
    % ======  Part I: Small window steps =======
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_span = 0.2;
    %         t_span = 1000;
    t_step = 0.05; % For time evolution
    toi = [0.95 1.3] ;   % Time of interests for weights illustration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    to_plot_linear = 1:N_sample_cell;
    
    % PSTH_correct_mean_headings = [cell_for_regression time stim_type abs_heading choice];
    
    t_centers = ts_Ann(1) + t_span / 2 : t_step : ts_Ann(end)-t_span / 2;
    
    weight_vest_vis = nan(length(to_plot_linear),2 * 2 + 2,length(t_centers)); % Weights (3), p_values (3), r^2, p-value of r^2
    
    
    % measure_pred_all = nan(length(to_plot_linear) * 2 * length(unique_abs_heading),2,length(t_centers));
    
    for i = 1:length(to_plot_linear)  % For each cell
        
        for tcc = 1:length(t_centers) % For each time bin
            
            t_range = t_centers(tcc) - t_span / 2 <= ts_Ann & ts_Ann <= t_centers(tcc) + t_span / 2;
            
            r = [];
            for k = 1:3
                
                % Use the raw spike count including correct AND incorrect trials. HH20171121
                r(:,k) = squeeze(mean(sum(spikes_lip(to_calculate_PSTH_cells_ind(to_plot_linear(i)),t_range,:,:,k),2),3));
                
                % Use correct only firing rate
                % r(:,k) = mean(reshape(squeeze(PSTH_correct_mean_headings(to_plot_linear(i),t_range,k,:,:)),sum(t_range),[]),1);
%                 r(:,k) = mean(reshape(squeeze(PSTH_all_mean_headings(to_plot_linear(i),t_range,k,:,:)),sum(t_range),[]),1);

            end
            
            
            % if tcc == 17    keyboard; end
            
            % GLM fit
%             [b,~,stat] = glmfit(r(:,1:2) ,r(:,3),'normal','link','identity','constant','on');
%             [~,~,~,~,stat_reg] = regress(r(:,3),[ones(size(r,1),1) r(:,1:2)]);
%             
%             weight_vest_vis(i,:,tcc) = [b(2:3)' stat.p(2:3)' stat_reg([1 3])]; % Weights, r^2 and p-value

            r = r - repmat(nanmean(r,1),size(r,1),1); % Mean removed for each
            X = r(:,1:2); Y = r(:,3);
            [b,~,stat] = glmfit(X, Y,'normal','link','identity','constant','off');
            
            [~,~,~,~,stat_reg] = regress(Y,[ones(size(r,1),1) X]);
            weight_vest_vis(i,:,tcc) = [b' stat.p' stat_reg([1 3])]; % Weights, r^2 and p-value
           
            % yfit = r(:,1:2) * b;
            % measure_pred_all((i-1)*10 + 1:(i-1)*10+size(r,1),:,tcc) = [r(:,3) yfit];
            
            % r = r([end:-2:2 1:2:end-1],:);
            % yfit = glmval(b,r(:,1:2),'identity');
            % figure(348); hold on; plot(1:10,r(:,3),'o',1:10,yfit,'r-');
            
        end
    end
    
    %%
    figure(7252243);
    set(gcf,'uni','norm','pos',[0.005       0.057       0.958       0.822]);clf; hold on;
    subplot(2,5,[1]);
    mean_paras = squeeze(mean(weight_vest_vis(:,[1 2 end-1],:),1))';
    sem_paras = squeeze(std(weight_vest_vis(:,[1 2 end-1],:)))'/sqrt(size(weight_vest_vis,1));
    
    temp_col = {'b','r','k'};
    temp_marker = {'s','o','o'}; % {'','',''};
    for pp = 1 : 3
        h = shadedErrorBar(repmat(t_centers',1,1),mean_paras(:,pp),sem_paras(:,pp),'lineProp',{[temp_marker{pp} temp_col{pp} '-'],'linew',2,'markersize',10});
        hold on;
    end
    axis tight;
    
    plot(repmat(t_centers',1,1),mean_paras(:,1)+mean_paras(:,2),'g','linew',2);
    plot(xlim,[1 1],'g:','linew',3);
    
    set(legend('w_{vest}','w_{vis}','R^2'),'location','best'); axis tight; %ylim([0 max(ylim)])
    xlabel(sprintf('Center of %g ms window',t_span)); ylabel('Median');
    
    % Gaussian vel
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    
    % ======  Part II: Large window, R^2 and weight. Figure 5 of Gu 2008  =====
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toi = [1 1.4] ;   % Time of interests for weights illustration
    t_span_large = 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_centers = toi;
    
    to_plot_linear = 1:N_sample_cell;
    % PSTH_correct_mean_headings = [cell_for_regression time stim_type abs_heading choice];
    
    weight_vest_vis = nan(length(to_plot_linear),2 * 2 + 2,length(t_centers)); % Weights (3), p_values (3), r^2, p-value of r^2
   
    weight_vest_vis_test = nan(length(to_plot_linear),2,length(t_centers));
    measure_pred_all = nan(length(to_plot_linear) * 2 * length(unique_abs_heading),2,length(t_centers));

%     weight_vest_vis_lsqlin = nan(length(to_plot_linear),2,length(t_centers));
    
    for i = 1:length(to_plot_linear)  % For each cell
        
        for tcc = 1:length(t_centers) % For each time bin
            
            t_range = t_centers(tcc) - t_span_large / 2 <= ts_Ann & ts_Ann <= t_centers(tcc) + t_span_large / 2;
            
            r = [];
            for k = 1:3
                
                % Use the raw spike count including correct AND incorrect trials. HH20171121
                r(:,k) = squeeze(mean(sum(spikes_lip(to_calculate_PSTH_cells_ind(to_plot_linear(i)),t_range,:,:,k),2),3));
                
                % Use correct only firing rate
                % r(:,k) = mean(reshape(squeeze(PSTH_correct_mean_headings(to_plot_linear(i),t_range,k,:,:)),sum(t_range),[]),1);
%                 r(:,k) = mean(reshape(squeeze(PSTH_all_mean_headings(to_plot_linear(i),t_range,k,:,:)),sum(t_range),[]),1);
            end
            
            
            % if tcc == 17    keyboard; end
            
            % GLM fit

%             [b,~,stat] = glmfit(r(:,1:2) ,r(:,3),'normal','link','identity','constant','on');
%             b = b(2:3);   stat.p = stat.p(2:3);            
            
            r = r - repmat(nanmean(r,1),size(r,1),1); % Mean removed for each
            X = r(:,1:2); Y = r(:,3);
            [b,~,stat] = glmfit(X ,Y,'normal','link','identity','constant','off');

            % Other linear fittings
%             b_test = pinv(X)*Y;
%             weight_vest_vis_test(i,:,tcc) = b_test';  % Pseudoinverse
            
            lamda = 10;
            b_test = inv(X'*X + lamda * eye(size(X,2)))*X'*Y; % Ridge regularization
            weight_vest_vis_test(i,:,tcc) = b_test';  % The same result for sure.
    
            [~,~,~,~,stat_reg] = regress(r(:,3),[ones(size(r,1),1) r(:,1:2)]);
            weight_vest_vis(i,:,tcc) = [b' stat.p' stat_reg([1 3])]; % Weights, r^2 and p-value
            
             
            
%             b_lsqlin = lsqlin([r(:,1:2)],r(:,3),[],[],[],[],[0 0]');
%             lin = lsqlin([ones(size(r(:,1))) r(:,1:2)],r(:,3),[],[],[],[],[0 0 0]'); b_lsqlin = b_lsqlin(2:3);
            
%             weight_vest_vis_lsqlin(i,:,tcc) = b_lsqlin;
            
            yfit = r(:,1:2) * b;
            measure_pred_all((i-1)*10 + 1:(i-1)*10+size(r,1),:,tcc) = [r(:,3) yfit];
            
            % r = r([end:-2:2 1:2:end-1],:);
            % yfit = glmval(b,r(:,1:2),'identity');
            % figure(348); hold on; plot(1:10,r(:,3),'o',1:10,yfit,'r-');
            
        end
    end
    
    % --- Plotting ---
    for toii = 1:2
        [~,toi_ind] = min(abs(toi(toii)-t_centers));
        
        sig_ind = weight_vest_vis(:,end,toi_ind)<0.05;
        nsig_ind = weight_vest_vis(:,end,toi_ind)>=0.05;
        
        % Predicted v.s. measured
        ax0 = subplot(2,5,2 + (toii-1)*5);
        h = LinearCorrelation( measure_pred_all(:,1,toi_ind), measure_pred_all(:,2,toi_ind),...
            'CombinedIndex',[],...
            'Xlabel','Measured response','Ylabel','Predicted response',...
            'FaceColors',{[0.9 0.9 0.9]},'Markers',{'.'},...
            'LineStyles',{'k-'},'MarkerSize',10,...
            'XHist',0,'YHist',0,...
            'XHistStyle','grouped','YHistStyle','grouped','SameScale',1,'Method','Spearman','axes',ax0);
        text(min(xlim)+2,min(ylim)+2,sprintf('r^2 = %g, p = %g',h.group(1).r_square,h.group(1).p),'FontSize',13);
        legend off;
        
        % Distribution of R^2
        ax1 = subplot(2,5,3 + (toii-1)*5);
        h = HistComparison({weight_vest_vis(sig_ind,end-1,toi_ind),weight_vest_vis(nsig_ind,end-1,toi_ind)},...
            'EdgeColors',{'k','k'},'FaceColors',{[.3 .3 .3 ],[0.9 0.9 0.9]},'XCenters',0.05:0.1:1,'MeanType','Median','Axes',ax1);
        xlabel('R^2'); ylabel('Number of neurons'); xlim([0 1]); ylim([0 23]);
        title(sprintf('t_{center} = %g (%g s window)',t_centers(toi_ind),t_span_large));
        
        % Vestibular / visual weights
        ax2(toii) = subplot(2,5,4 + (toii-1)*5); 
        %    plot(weight_vest_vis(:,1,toi_ind),weight_vest_vis(:,2,toi_ind),'o');
        
        %     h = LinearCorrelation({
        %         weight_vest_vis(nsig_ind,1,toi_ind);
        %         weight_vest_vis(sig_ind,1,toi_ind);
        %         },...
        %         {
        %         weight_vest_vis(nsig_ind,2,toi_ind) ;
        %         weight_vest_vis(sig_ind,2,toi_ind) ;
        %         },...
        %         'CombinedIndex',[],...
        %         'Xlabel','Vestibular weight','Ylabel','Visual weight',...
        %         'FaceColors',{[0.9 0.9 0.9],[.3 .3 .3 ]},'Markers',{'o'},...
        %         'LineStyles',{'k--','k-'},'MarkerSize',9,...
        %         'XHist',20,'YHist',10,...
        %         'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,'Method','Spearman','axes',ax2);
        % text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(2).r_square,h.group(2).p),'FontSize',13);
        
        % Fit color to to N_sample_cell/4
        cols = colormap(jet);
        cols = interp1(1:length(cols),cols,1:length(cols)/(N_sample_cell/4+4):length(cols));
        for cc = 1:N_sample_cell
            dis_prop = min(abs(cc-N_sample_cell/4),abs(cc-N_sample_cell/4*3));
            this_col = cols(1+round(dis_prop),:);
            plot(weight_vest_vis(cc,1,toi_ind),weight_vest_vis(cc,2,toi_ind),'o','markersize',7,'color',this_col,'linewid',2); hold on;
        end
        xlims = xlim; ylims = ylim;
        xlims_new = [min(xlims(1),ylims(1)) max(xlims(2),ylims(2))];
        xlims_new = [xlims_new(1)-range(xlims_new)/20 xlims_new(2)+range(xlims_new)/20];
        
        axis([xlims_new xlims_new]);
        h.diag = plot([xlims_new(1) xlims_new(2)],[xlims_new(1) xlims_new(2)],'k:');
        
        % Annotate significant regressions
        sig_ind = weight_vest_vis(:,end,toi_ind)<0.05;
        plot(weight_vest_vis(sig_ind,1,toi_ind),weight_vest_vis(sig_ind,2,toi_ind),'o','markersize',7,'markerfacecolor','k');
        
        if toii == 1 ; xlabel(''); end
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
        legend off;
        set(gca,'xtick',-10:1:10,'ytick',-10:1:10); axis square
        
%         ax3(toii) = subplot(2,5,5 + (toii-1)*5);
%         plot(ax3(toii),weight_vest_vis_lsqlin(:,1,toi_ind),weight_vest_vis_lsqlin(:,2,toi_ind),'o'); hold on
%         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
%         set(gca,'xtick',-10:1:10,'ytick',-10:1:10);
        
        % ------  Test plot ----------

        ax3(toii) = subplot(2,5,5 + (toii-1)*5);
        
        % Fit color to to N_sample_cell/4
        cols = colormap(jet);
        cols = interp1(1:length(cols),cols,1:length(cols)/(N_sample_cell/4+4):length(cols));
        for cc = 1:N_sample_cell
            dis_prop = min(abs(cc-N_sample_cell/4),abs(cc-N_sample_cell/4*3));
            this_col = cols(1+round(dis_prop),:);
            plot(weight_vest_vis_test(cc,1,toi_ind),weight_vest_vis_test(cc,2,toi_ind),'o','markersize',7,'color',this_col,'linewid',2); hold on;
        end
        xlims = xlim; ylims = ylim;
        xlims_new = [min(xlims(1),ylims(1)) max(xlims(2),ylims(2))];
        xlims_new = [xlims_new(1)-range(xlims_new)/20 xlims_new(2)+range(xlims_new)/20];
        
        axis([xlims_new xlims_new]);
        h.diag = plot([xlims_new(1) xlims_new(2)],[xlims_new(1) xlims_new(2)],'k:');

        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--'); plot(xlim,[1 1],'k--'); plot([1 1],ylim,'k--');
        legend off;
        set(gca,'xtick',-10:1:10,'ytick',-10:1:10); axis square
        title(sprintf('Ridge lamda = %g', lamda))


    end
    
    if ION_cluster
        file_name = sprintf('./result/%s5_LinearRegression%s',save_folder,para_override_txt);
        export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
        h_grouped = [h_grouped ax2(1)]; % Add h_grouped
    end
end

%% Information of sensory areas
if analysis_switch(9)
    
    %% --------------- Get spike counts (in Hz) -----------------
    info_win_wid = 200e-3; % 100e-3; % in s
    info_win_step = 20e-3; % 100e-3; 
    n_info_win = round((trial_dur_total-info_win_wid)/info_win_step)+1;
    
    spikecount_dt_lip = nan(N_lip,n_info_win,N_rep,length(unique_heading),length(unique_stim_type));
    spikecount_t_vest = nan(N_vest,n_info_win,N_rep,length(unique_heading),length(unique_stim_type));
    spikecount_t_vis = nan(N_vis,n_info_win,N_rep,length(unique_heading),length(unique_stim_type));
    info_ts = nan(1,n_info_win);
    
    for tt = 1:n_info_win
        this_win_beg = (tt-1)*info_win_step + ts(1);
        this_win_end = this_win_beg + info_win_wid;
        info_ts(tt) =  mean([this_win_beg this_win_end]);
        
        % Info in [t,t+dt] of LIP
        spikecount_dt_lip(:,tt,:,:,:) = sum(spikes_lip(:,this_win_beg <= ts & ts < this_win_end,:,:,:),2);
        
        % Try rate instead (mean firing rate over dt)
%         spikecount_dt_lip(:,tt,:,:,:) = sum(rate_lip(:,this_win_beg <= ts & ts < this_win_end,:,:,:),2);
        
        % Info in [0,t+dt] of vest/vis
        spikecount_t_vest(:,tt,:,:,:) = sum(spikes_vest(:,0 <= ts & ts < this_win_end,:,:,:),2); 
        spikecount_t_vis(:,tt,:,:,:) = sum(spikes_vis(:,0 <= ts & ts < this_win_end,:,:,:),2);
       
        % Try rate instead (mean firing rate over 0->t)
%         spikecount_t_vest(:,tt,:,:,:) = sum(rate_vest(:,0 <= ts & ts < this_win_end,:,:,:),2); 
%         spikecount_t_vis(:,tt,:,:,:) = sum(rate_vis(:,0 <= ts & ts < this_win_end,:,:,:),2);
    end
    
    train_ratio = 0.5;
    
    infos_dt_lip = zeros(n_info_win,3);
    infos_t_vest = zeros(n_info_win,3);
    infos_t_vis = zeros(n_info_win,3);
    infos_dt_lip_simpleGu = zeros(n_info_win,3); infos_dt_lip_simpleGu_SE = infos_dt_lip_simpleGu;
    infos_dt_lip_partialSensoryFI = zeros(n_info_win,3);   infos_dt_lip_partialSensoryFI_SE = zeros(n_info_win,3);   
    infos_dt_lip_partialChoiceFI = zeros(n_info_win,3);   infos_dt_lip_partialChoiceFI_SE = zeros(n_info_win,3);   
    
    disp('Calculate information...')
    headings = reshape(repmat(unique_heading,N_rep,1),[],1);

    % ------------ Train the decoder once using the most informative (the last one) spike counts -----------------
%     [~,svm_model_lip_count] = fisher_HH(reshape(spikecount_dt_lip(:,end,:,:,3),N_lip,[])',headings,train_ratio);
%     [~,svm_model_vest_count] = fisher_HH(reshape(spikecount_t_vest(:,end,:,:,3),N_vest,[])',headings,train_ratio);
%     [~,svm_model_vis_count] = fisher_HH(reshape(spikecount_t_vis(:,end,:,:,3),N_vis,[])',headings,train_ratio);
    
    %% ------------ Info in each small window --------------
    for tt = 1:n_info_win
        for kk = 1:3
            
            activities = reshape(spikecount_dt_lip(:,tt,:,:,kk),N_lip,[])';
%             infos_dt_lip(tt,kk) = fisher_HH(activities,headings,train_ratio,svm_model_lip_count);
            
            % ---- Add FisherSimpleGu as in the experiment data. HH20180622 03:18 Argentina 0:3 Croatia... ----
            [infos_dt_lip_simpleGu(tt,kk), infos_dt_lip_simpleGu_SE(tt,kk)] = ...
                        fisher_HH_simpleGu(activities(:,to_calculate_PSTH_cells_ind),headings);
                    
            % ---- Add Partial FisherInfo. HH20180824
            this_result = fisher_HH_partialFI(activities(:,to_calculate_PSTH_cells_ind),headings,reshape(choices(:,:,kk),[],1));
            
            infos_dt_lip_partialSensoryFI(tt,kk) = this_result.infoPartialSensory;
            infos_dt_lip_partialSensoryFI_SE(tt,kk) = this_result.bootSESensory;
            infos_dt_lip_partialChoiceFI(tt,kk) = this_result.infoPartialChoice;
            infos_dt_lip_partialChoiceFI_SE(tt,kk) = this_result.bootSEChoice;
                        
            
            if kk ~= 2
                activities = reshape(spikecount_t_vest(:,tt,:,:,kk),N_vest,[])';
%                 infos_t_vest(tt,kk) = fisher_HH(activities,headings,train_ratio,svm_model_vest_count);
            end
            
            if kk ~= 1
                activities = reshape(spikecount_t_vis(:,tt,:,:,kk),N_vis,[])';
%                 infos_t_vis(tt,kk) = fisher_HH(activities,headings,train_ratio,svm_model_vis_count);
            end
        end
    end
    
    disp('Calculate information... Done')
    
    figure(1446);   clf;
    set(gcf,'uni','norm','pos',[0.013       0.079       0.875       0.733]);    
    
    subplot(2,3,1);
    plot(info_ts,infos_dt_lip(:,1),'bo-', info_ts, infos_t_vest(:,1),'b--'); hold on; title('vest');
    
    subplot(2,3,2); 
    plot(info_ts,infos_dt_lip(:,2),'ro-', info_ts, infos_t_vis(:,2),'r--'); hold on; title('vis');

    h_3 = subplot(2,3,3);
    plot(info_ts,infos_dt_lip(:,3),'go-', info_ts, infos_t_vest(:,3),'b--', ...
            info_ts, infos_t_vis(:,3),'r--', ...
            info_ts, infos_t_vis(:,3) + infos_t_vest(:,3),'k--',...
            info_ts, infos_dt_lip(:,1) + infos_dt_lip(:,2),'k-'); hold on ; title('comb');
   
    linkaxes(findobj(gcf,'type','axes'),'xy');
    
    try
        h_4 = subplot(2,3,4);
        SeriesComparison(shiftdim(infos_dt_lip_simpleGu,-1),info_ts,'OverrideError',infos_dt_lip_simpleGu_SE,'axes',h_4);
        hold on; plot(info_ts,sum(infos_dt_lip_simpleGu(:,1:2),2),'m','linew',2);
        xlim([0 1.5]); legend off; title('Standard FI (Gu 2010)')
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
        h_5 = subplot(2,3,5);
        SeriesComparison(shiftdim(infos_dt_lip_partialSensoryFI,-1),info_ts,'OverrideError',infos_dt_lip_partialSensoryFI_SE,'axes',h_5);
        hold on; plot(info_ts,sum(infos_dt_lip_partialSensoryFI(:,1:2),2),'m','linew',2);
        xlim([0 1.5]); legend off; title('Partial Sensory FI (Gu)')
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
        linkaxes([h_4 h_5],'xy')
        
        h_6 = subplot(2,3,6);
        SeriesComparison(shiftdim(infos_dt_lip_partialChoiceFI,-1),info_ts,'OverrideError',infos_dt_lip_partialChoiceFI_SE,'axes',h_6);
        hold on; plot(info_ts,sum(infos_dt_lip_partialChoiceFI(:,1:2),2),'m','linew',2);
        xlim([0 1.5]); legend off; title('Partial Choice FI (Gu)')
        plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
        
    end
    
    SetFigure(20)    
    
    if ION_cluster
        file_name = sprintf('./result/%s6_Info%s',save_folder,para_override_txt);
        export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
        h_grouped = [h_grouped h_3]; % Add h_grouped
    end
  
end

%%
analysis_time = toc(analysis_time)
disp('---- All done ----')
