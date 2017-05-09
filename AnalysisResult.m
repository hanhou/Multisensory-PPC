
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
ts = ((0:trial_dur_total_in_bins-1)-stim_on_time_in_bins)*dt; % in s
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

%% == Actual shifting window spike count for LIP as in experimental data (~ rectified rate_lip etc.) ===
% I don't use this. See above.
%{
win_wid = 100e-3; % in s
win_step = 10e-3;
real_rate_len = fix((trial_dur_total-win_wid)/win_step)+1;

real_rate_lip = nan(N_lip,real_rate_len,N_rep,length(unique_heading),length(unique_stim_type));
real_ts = nan(1,real_rate_len);

for tt = 1:real_rate_len
    this_win_beg = (tt-1)*win_step + ts(1); 
    this_win_end = this_win_beg + win_wid;
    real_ts(tt) =  mean([this_win_beg this_win_end]);
    
    real_rate_lip(:,tt,:,:,:) = sum(spikes_lip(:,this_win_beg <= ts & ts < this_win_end,:,:,:),2)/win_wid;
end

figure();
plot(real_ts,mean(real_rate_lip(left_targ_ind,:,:,end,1),3),'b--'); hold on;
plot(real_ts,mean(real_rate_lip(right_targ_ind,:,:,end,1),3),'b')
plot(ts,mean(rectified_rate_lip(left_targ_ind,:,:,end,1),3),'r--');
plot(ts,mean(rectified_rate_lip(right_targ_ind,:,:,end,1),3),'r');
plot(ts,mean(rate_lip(left_targ_ind,:,:,end,1),3),'g--');
plot(ts,mean(rate_lip(right_targ_ind,:,:,end,1),3),'g');
export_fig('-painters','-nocrop','-m2' ,sprintf('./result/test.png'));
%}

%% == Plotting example network activities ==
this_heading_ind = length(unique_heading);
to_plot_cond = length(unique_stim_type);

if if_debug
%%  ====== Animation ======
    rate_real_vis_aver = nanmean(spikes_vis(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_vest_aver = nanmean(spikes_vest(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_expected_int_aver = nanmean(rate_int(:,:,:,this_heading_ind,to_plot_cond),3);
    rate_expected_lip_aver = nanmean(rate_lip(:,:,:,this_heading_ind,to_plot_cond),3);
    rate_real_int_aver = nanmean(spikes_int(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_lip_aver = nanmean(spikes_lip(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    rate_real_targ_aver = nanmean(spikes_target(:,:,:,this_heading_ind,to_plot_cond),3)/dt;
    
    figure(1001); clf
    
    for ttt = 1:10:length(ts)
        
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
        axis(gca,[min(prefs_lip) max(prefs_lip) min(rate_expected_lip_aver(:)) max(rate_expected_lip_aver(:))*2]);
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
end

%% ====== Connection and Raster plot ======
if if_debug
    figure(90);
    
    subplot(4,3,[2 3]);
    imagesc(ts,prefs_lip,rate_real_lip_aver*dt); axis xy;
    ylabel('LIP');
    title(sprintf('Firing prob. for each bin, averaged over %g trials, cond = %g, heading = %g',...
        N_rep,unique_stim_type(to_plot_cond),unique_heading(this_heading_ind)));  colorbar
    
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

%% ====== Readout LIP activity and make decisions =====

rate_lip_at_decision = squeeze(rate_lip(:,end,:,:,:)); % Get the population activity at the end of the trials
[~, readout_lip_at_decision] = max(rate_lip_at_decision,[],1); % Readout population activity    
choices = sign((prefs_lip(shiftdim(readout_lip_at_decision)) >= 0) -0.5); % choices(reps,headings,stimtypes):  1 = rightward, -1 = leftward
    
%     % I just flip the psychometric curve to the negative headings
%     psychometric = [unique_heading' sum(reshape(choices{ss},[],length(unique_heading)),1)'/N_rep];
%     fake_psy = flipud(psychometric(unique_heading>0,:));
%     fake_psy(:,1) = -fake_psy(:,1);
%     fake_psy(:,2) = 1-fake_psy(:,2);
%     psychometric = [fake_psy ; psychometric];

    
%% ====== Example Pref and Null traces (correct only) ======
%%{

to_plot_abs_headings = [0 4 8];
to_plot_cell_ind = right_targ_ind;  % Cell-based plotting
n_to_plot_trace = 5;

set(figure(1001),'name','Example PSTHs (correct only)'); clf;
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
            this_correct_ind = find(choices(:,this_heading_ind,ss) == pref_null(cc)); % Only correct trials
            
            for tt = 1:ceil(length(this_correct_ind)/n_to_plot_trace):length(this_correct_ind)
                
                RT_this = RT(this_correct_ind(tt),this_heading_ind,ss);
                rate_at_RT = rate_lip(to_plot_cell_ind,find(ts<=RT_this,1,'last'),this_correct_ind(tt),this_heading_ind,ss);
                
                % --- LIP ---
                axes(hs(6*(tph-1)+ss));
                
                if cc == 1 % Pref
                    plot(ts,rate_lip(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'color',colors(unique_stim_type(ss),:),'linewid',2); hold on;
                    plot([RT_this RT_this],[rate_at_RT+5 rate_at_RT-5],'m','linew',5);
                else
                    plot(ts,rate_lip(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'k--','linewid',1);
                    plot([RT_this RT_this],[rate_at_RT+5 rate_at_RT-5],'k','linew',5);
                end
                
                axis tight; y_min = min(y_min,min(ylim)); y_max = max(y_max,max(ylim)); to_sync = [to_sync gca];
                
                % --- Int ---
                axes(hs(6*(tph-1)+ss)+3);
                if cc == 1 % Pref
                    plot(ts,rate_int(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'color',colors(unique_stim_type(ss),:),'linewid',2); hold on;
                else
                    plot(ts,rate_int(to_plot_cell_ind,:,this_correct_ind(tt),this_heading_ind,ss),'k--','linewid',1); hold on;
                end
                %     if if_bounded
                %         plot(xlim,[decis_thres(unique_stim_type(ss)) decis_thres(unique_stim_type(ss))],'c--');
                % %         ylim([min(ylim),decis_thres*1.1]);
                %     end
                
            end
        end
        
        axes(hs(6*(tph-1)+ss));
        title(sprintf('rate\\_lip, pref = %g, |heading| = %g',prefs_lip(to_plot_cell_ind), to_plot_abs_headings(tph)));
        plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
        axis tight;
        if if_bounded
            plot(xlim,[decis_thres(unique_stim_type(ss)) decis_thres(unique_stim_type(ss))],'c--');
            %         ylim([min(ylim),decis_thres(k)*1.1]);
        end
        
        axes(hs(6*(tph-1)+ss)+3);
        title(sprintf('rate\\_int'));
        plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
        axis tight;
        
    end
    
end

set(to_sync,'ylim',[y_min y_max]);

if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/1_Example%s.png',para_override_txt));
    saveas(gcf,sprintf('./result/1_Example%s.fig',para_override_txt),'fig');
end

disp('Example done');

%}

%% ====== Fig.2 Overview ======
figure(1000); clf;
set(gcf,'name','Overview');
set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);

hs = tight_subplot(3,4,0.05);

% ======= Behavior ========
% Psychometric
for ss = 1:length(unique_stim_type)
        
%     % I just flip the psychometric curve to the negative headings
%     psychometric = [unique_heading' sum(reshape(choices{ss},[],length(unique_heading)),1)'/N_rep];
%     fake_psy = flipud(psychometric(unique_heading>0,:));
%     fake_psy(:,1) = -fake_psy(:,1);
%     fake_psy(:,2) = 1-fake_psy(:,2);
%     psychometric = [fake_psy ; psychometric];

    psychometric = [unique_heading' sum(reshape(choices(:,:,ss)==RIGHT,[],length(unique_heading)),1)'/N_rep];

    axes(hs(1));
    plot(psychometric(:,1),psychometric(:,2),'o','color',colors(unique_stim_type(ss),:),'markerfacecolor',colors(unique_stim_type(ss),:),'markersize',10); % Psychometric
    hold on;
    xxx = -max(unique_heading):0.1:max(unique_heading);
    
    [bias, threshold] = cum_gaussfit_max1(psychometric);
    psycho(ss,:) = [bias,threshold];
    plot(xxx,normcdf(xxx,bias,threshold),'-','color',colors(unique_stim_type(ss),:),'linewid',4);
    axis([-8.5 8.5 0 1]);

    set(text(min(xlim),0.6+0.06*ss,sprintf('threshold = %g\n',threshold)),'color',colors(unique_stim_type(ss),:));
    title(sprintf('N_{reps} = %g',N_rep));
    
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
title(sprintf('RT, correct trials, Bound = %s',num2str(decis_thres)));

axes(hs(1));
if length(unique_stim_type) == 3
    pred_thres = (psycho(1,2)^(-2)+psycho(2,2)^(-2))^(-1/2);
    set(text(min(xlim),0.6+0.06*(ss+1),sprintf('pred = %g\n',pred_thres)),'color','k');
end

%}

% ===== Averaged PSTH, different angles =====
%%{

to_calculate_PSTH_cells_ind = [right_targ_ind-20:5:right_targ_ind+20 left_targ_ind-20:5:left_targ_ind+20];   
N_sample_cell = length(to_calculate_PSTH_cells_ind);
% to_plot_grouped_PSTH = find(to_calculate_PSTH_cells_ind == right_targ_ind);
to_plot_grouped_PSTH = 1:N_sample_cell; % Group all calculated cells

to_flip = prefs_lip(to_calculate_PSTH_cells_ind) < 0;
PSTH_correct_mean_headings = nan(N_sample_cell,length(ts),length(unique_stim_type),length(unique_abs_heading),2);
PSTH_correct_mean_allheading = nan(N_sample_cell,length(ts),length(unique_stim_type),2);
pref_null = [RIGHT LEFT];

% === Group PSTH data ===

for ss = 1:length(unique_stim_type)
    
    zu_all{ss} = []; % Cache for VarCE
    
    for cc = 1:2 % Pref and Null
        PSTH_cache = [];
        
        for abshh = 1:length(unique_abs_heading)         % For different abs(headings)
            % Assuming PREF = RIGHT for all the cells
            this_heading_ind = find(unique_heading == unique_abs_heading(abshh) * pref_null(cc),1);
            this_correct_ind = find(choices(:,this_heading_ind,ss) == pref_null(cc)); % Only correct trials
            
            % PSTH separated for each heading
            PSTH_correct_raw_this = rate_lip(to_calculate_PSTH_cells_ind,:,this_correct_ind,this_heading_ind,ss);
            PSTH_correct_mean_headings(:,:,ss,abshh,cc) = mean(PSTH_correct_raw_this,3);
            
            % PSTH grouped cache for all headings
            PSTH_cache = cat(3,PSTH_cache,PSTH_correct_raw_this);

            % For VarCE: all choice
            if ~(unique_abs_heading(abshh)==0 && cc==2) % Skip abs(heading) = 0 and cc = 2 because cc = 1 already includes all trials for 0 heading
                PSTH_all_raw_this_forVarCE = rate_lip(to_calculate_PSTH_cells_ind,:,:,this_heading_ind,ss);
                zu_this = bsxfun(@minus, PSTH_all_raw_this_forVarCE, mean(PSTH_all_raw_this_forVarCE,3));
                zu_all{ss} = cat(3,zu_all{ss},zu_this);
            end
        end
        
        % PSTH for all headings
        PSTH_correct_mean_allheading(:,:,ss,cc) = mean(PSTH_cache,3);
    
    end
end
% Flip pref and null if PREF = LEFT for some of the cells
PSTH_correct_mean_headings(to_flip,:,:,:,:) = flip(PSTH_correct_mean_headings(to_flip,:,:,:,:),5);
PSTH_correct_mean_allheading(to_flip,:,:,:) = flip(PSTH_correct_mean_allheading(to_flip,:,:,:),4);

diff_PSTH_correct_mean_headings = - diff(PSTH_correct_mean_headings,[],5);
diff_PSTH_correct_mean_allheading = - diff(PSTH_correct_mean_allheading,[],4);

% -- Plotting --
y_max = -inf; y_min = inf;
colors_hue = [240 0 120]/360;

for ss = 1:length(unique_stim_type)
    yyy = squeeze(reshape(mean(PSTH_correct_mean_headings(to_plot_grouped_PSTH,:,ss,:,:),1),1,[],1,length(unique_abs_heading)*2));
    yyy = shiftdim(yyy,-1);
    
    colors_angles_hsv(:,2) = linspace(0.2,1,length(unique_abs_heading));
    colors_angles_hsv(:,1) = colors_hue(ss);
    colors_angles_hsv(:,3) = 0.9;
    colors_angles = hsv2rgb(colors_angles_hsv);
    colors_angles = mat2cell(colors_angles,ones(length(unique_abs_heading),1));
    
    h = SeriesComparison(yyy,ts,...
        'Colors',colors_angles,'LineStyles',[repmat({'-'},1,length(unique_abs_heading)) repmat({'--'},1,length(unique_abs_heading))],...
        'ErrorBar',0,'Xlabel',[],'axes',hs(3+ss));
    legend off;
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');

    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    xlim([min(ts),max(ts)]);
end

% title(hs(4),sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(to_plot_grouped_PSTH))));
title(hs(4),sprintf('Grouped %g cells',length(to_plot_grouped_PSTH)));

set(hs(4:6),'ylim',[y_min y_max]);

% =====  Delta-PSTH (stim type), different heading, correct only ====
y_max = -inf; y_min = inf;

for abshh = 1:length(unique_abs_heading)
    
    axes(hs(6+abshh)); hold on;
    
    for ss = 1:length(unique_stim_type)
        plot(ts,mean(diff_PSTH_correct_mean_headings(to_plot_grouped_PSTH,:,ss,abshh),1),'color',colors(unique_stim_type(ss),:),'linew',2.5);
    end
    plot(ts,mean(sum(diff_PSTH_correct_mean_headings(to_plot_grouped_PSTH,:,[1 2],abshh),3),1),'k--');
    title(sprintf('|heading| = %g',unique_abs_heading(abshh)));
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
end

% ======= Delta-PSTH, all heading, correct only =====
axes(hs(12)); hold on;    title('All headings');

for ss = 1:length(unique_stim_type)
    plot(ts,mean(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,ss),1),'color',colors(unique_stim_type(ss),:),'linew',2.5);
end
plot(ts,mean(sum(diff_PSTH_correct_mean_allheading(to_plot_grouped_PSTH,:,[1 2]),3),1),'k--');
plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');

axis tight;
y_min = min(y_min,min(ylim));
y_max = max(y_max,max(ylim));

set(hs(7:12),'ylim',[y_min y_max]);


% ====== VarCE ======
axes(hs(3)); hold on; title('VarCE');

for ss = 1:length(unique_stim_type)
    plot(ts,mean(var(zu_all{ss}(to_plot_grouped_PSTH,:,:)*0.06,[],3),1),...  % 60 ms spike counting window in Ann's paper
        'color',colors(unique_stim_type(ss),:),'linew',2.5);
end
plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
xlim([min(ts),max(ts)]);


if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/2_Overview%s.png',para_override_txt));
    saveas(gcf,sprintf('./result/2_Overview%s.fig',para_override_txt),'fig');
end
disp('Overview done');

%% ====== Fig.3 Different cells ======

figure(1002); clf;
set(gcf,'name','Different cells');
set(gcf,'uni','norm','pos',[0.014        0.06       0.895       0.829]);

% -- diff_PSTH for all headings --
m = fix(sqrt(N_sample_cell+1));
n = ceil((N_sample_cell+1)/m);
[~,hs] = tight_subplot(m,n,0.05);
y_max = -inf; y_min = inf;

for cc = 1:N_sample_cell
    axes(hs(cc)); hold on;
    for ss = 1:length(unique_stim_type)
        plot(ts,diff_PSTH_correct_mean_allheading(cc,:,ss),'color',colors(unique_stim_type(ss),:),'linew',2.5);
    end
    plot(ts,sum(diff_PSTH_correct_mean_allheading(cc,:,[1 2]),3),'k--');
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    
    axis tight;
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    
    title(sprintf('pref = %g',prefs_lip(to_calculate_PSTH_cells_ind(cc))));
end

set(hs,'ylim',[y_min y_max]);

% -- Raster plot of multisensory enhancement --
axes(hs(end));
comb_minus_vest = mean(diff_PSTH_correct_mean_allheading(:,:,3)-diff_PSTH_correct_mean_allheading(:,:,1),2);
comb_minus_vis = mean(diff_PSTH_correct_mean_allheading(:,:,3)-diff_PSTH_correct_mean_allheading(:,:,2),2);
plot(comb_minus_vest,comb_minus_vis,'og','markersize',7,'markerfacecol','g'); hold on;
xlabel('Comb - Vest'); ylabel('Comb - Vis');
max_axis = max(abs([ylim xlim]));
plot([0 0],[-max_axis max_axis],'--k');
plot([-max_axis max_axis],[0 0],'--k');
plot([-max_axis max_axis],[-max_axis max_axis],'--k');
axis tight square;


if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/3_Cells%s.png',para_override_txt));
    saveas(gcf,sprintf('./result/3_Cells%s.fig',para_override_txt),'fig');
end
disp('Different cells done');


%%

analysis_time = toc(analysis_time)
disp('---- All done ----')
