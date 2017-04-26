%% Get data
colors = [0 0 1; 1 0 0; 0 0.8 0.4];

%% == Plotting example network activities ==
to_plot_heading = length(unique_heading);
to_plot_cond = length(unique_condition);

% %{
%%  ====== Animation ======
if if_debug
    
    rate_real_vis_aver = nanmean(spikes_vis(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
    rate_real_vest_aver = nanmean(spikes_vest(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
    rate_expected_int_aver = nanmean(rate_int(:,:,:,to_plot_heading,to_plot_cond),3);
    rate_expected_lip_aver = nanmean(rate_lip(:,:,:,to_plot_heading,to_plot_cond),3);
    rate_real_int_aver = nanmean(spikes_int(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
    rate_real_lip_aver = nanmean(spikes_lip(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
    rate_real_targ_aver = nanmean(spikes_target(:,:,:,to_plot_heading,to_plot_cond),3)/dt;
    
    figure(1001); clf
    
    for ttt = 1:10:length(ts)
        
        % INTEGRATOR layer
        subplot(2,2,1);
        plot(prefs_int,rate_expected_int_aver(:,ttt)); hold on;
        plot(prefs_int,rate_real_int_aver(:,ttt),'r');  hold off;
        axis(gca,[min(prefs_int) max(prefs_int) min(rate_expected_int_aver(:)) max(rate_expected_int_aver(:))*2]);
        title(sprintf('Integrator, heading = %g, cond = %g, aver %g trials',unique_heading(to_plot_heading),unique_condition(to_plot_cond),N_rep));
        
        % LIP layer
        subplot(2,2,2);
        plot(prefs_lip,rate_expected_lip_aver(:,ttt)); hold on;
        plot(prefs_lip,rate_real_lip_aver(:,ttt),'r');  hold off;
        axis(gca,[min(prefs_lip) max(prefs_lip) min(rate_expected_lip_aver(:)) max(rate_expected_lip_aver(:))*2]);
        title(sprintf('LIP, t = %g',ttt*dt*1000));
        
        % Visual layer
        subplot(2,2,3);
        %     plot(prefs_int,aux_proba_vis(:,ttt)/dt); hold on;
        %     plot(prefs_int,rate_real_vis_aver(:,ttt),'r'); hold off;
        %     axis(gca,[min(prefs_int) max(prefs_int) min(aux_proba_vis(:)/dt) max(aux_proba_vis(:)/dt)*2]);
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

%}

%% ====== Raster plot ======
if if_debug
    figure(90);
    
    subplot(4,3,[2 3]);
    imagesc(ts,prefs_lip,rate_real_lip_aver*dt); axis xy;
    ylabel('LIP');
    title(sprintf('Firing prob. for each bin, averaged over %g trials, cond = %g, heading = %g',...
        N_rep,unique_condition(to_plot_cond),unique_heading(to_plot_heading)));  colorbar
    
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
if if_debug || 1
    
    to_plot_headings = [0 8];
    n_to_plot_trace = 5;
    
    for tph = 1:length(to_plot_headings)
        
        to_plot_heading  = find(to_plot_headings(tph) == unique_heading,1);

        set(figure(1002+tph),'name','Example PSTHs'); clf;
        set(gcf,'uni','norm','pos',[0.003       0.039       0.555       0.878]);
        
        for cc = 1:length(unique_condition)
            % --- LIP ---
            subplot(3,2,2*unique_condition(cc)-1);
            
            for trialtoplot = 1:ceil(N_rep/n_to_plot_trace):N_rep
                plot(ts,rate_lip(right_targ_ind,:,trialtoplot,to_plot_heading,cc),'color',colors(unique_condition(cc),:),'linewid',2); hold on;
                plot(ts,rate_lip(left_targ_ind,:,trialtoplot,to_plot_heading,cc),'--k','linewid',1);
            end
            plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
            axis tight;
            
            if if_bounded
                plot(xlim,[decis_thres(unique_condition(cc)) decis_thres(unique_condition(cc))],'c--');
                %         ylim([min(ylim),decis_thres(k)*1.1]);
            end
            
            title(sprintf('rate\\_lip, heading = %g, coh = %g',unique_heading(to_plot_heading),coherence));
            
            % --- Int ---
            subplot(3,2,2*unique_condition(cc));
            
            for trialtoplot = 1:ceil(N_rep/10):N_rep
                plot(ts,rate_int(right_targ_ind,:,trialtoplot,to_plot_heading,cc),'color',colors(unique_condition(cc),:),'linewid',2); hold on;
                plot(ts,rate_int(left_targ_ind,:,trialtoplot,to_plot_heading,cc),'--k','linewid',1);
            end
            plot(t_motion,vel/max(vel)*max(ylim)/3,'k--');
            axis tight;
            
            %     if if_bounded
            %         plot(xlim,[decis_thres(unique_condition(cc)) decis_thres(unique_condition(cc))],'c--');
            % %         ylim([min(ylim),decis_thres*1.1]);
            %     end
            
            title(sprintf('rate\\_int, heading = %g, coh = %g',unique_heading(to_plot_heading),coherence));
            
            
            %     subplot(3,2,2);
            %     plot(ts,nanmean(rate_lip(right_targ_ind,:,:,to_plot_heading,cc),3)...
            %         -nanmean(rate_lip(left_targ_ind,:,:,to_plot_heading,cc),3),'color',colors(unique_condition(cc),:),'linew',2);
            %     hold on;
            %     title(sprintf('averaged of all %g trials (correct + wrong), pref-null',N_rep));
            %     axis tight;
        end
        
        if ION_cluster
            export_fig('-painters','-nocrop','-m2' ,sprintf('./result/example_%gdegree_%s.png',to_plot_headings(tph),hostname));
        end
        
    end
end

%}

% ====== Behavior performance ======
%%{
figure(99); clf; hold on; axis([-8 8 0 1]);
set(gcf,'name','Behavior','uni','norm','pos',[0.32       0.071       0.357        0.83]);

% Psychometric
for cc = 1:length(unique_condition)
    subplot(2,1,1);
    
    % Make decisions
    rate_int_at_decision = squeeze(rate_lip(:,end,:,:,cc)); % Using lip
    
    [~, pos_max_rate_int_at_decision] = max(rate_int_at_decision,[],1);
    choices{cc} = prefs_int(shiftdim(pos_max_rate_int_at_decision)) >= 0; % 1 = rightward, 0 = leftward
    
    % I just flip the psychometric curve to the negative headings
    psychometric = [unique_heading' sum(reshape(choices{cc},[],length(unique_heading)),1)'/N_rep];
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
    
    % Chronometric
    subplot(2,1,2);
    for hh = 1:length(unique_heading)
        select_correct = (choices{cc}(:,hh) == (unique_heading(hh)>=0)) | unique_heading(hh)== 0 ; % Only correct trials should be plotted in RT
        
        RT_this_cc_mean(hh,1) = mean(RT(select_correct,hh,cc),1);
        RT_this_cc_se(hh,1) = std(RT(select_correct,hh,cc),[],1);
        
        RT_this_cc_mean_wrong(hh,1) = mean(RT(~select_correct,hh,cc),1);
        RT_this_cc_se_wrong(hh,1) = std(RT(~select_correct,hh,cc),[],1);
        
    end
    
    errorbar(unique_heading+(cc-2)*0.2,RT_this_cc_mean,RT_this_cc_se,'o-','color',colors(unique_condition(cc),:),...
        'markerfacecolor',colors(unique_condition(cc),:),'markersize',10,'linew',2);
    hold on;
    errorbar(unique_heading+(cc-2)*0.2,RT_this_cc_mean_wrong,RT_this_cc_se_wrong,'o--','color',colors(unique_condition(cc),:),...
        'markerfacecolor',colors(unique_condition(cc),:),'markersize',7,'linew',1);
end

subplot(2,1,2);
plot(vel/max(vel)*max(xlim)/5,t_motion,'k--');
title(sprintf('RT, all trials (correct), Bound = %s',num2str(decis_thres)));

subplot(2,1,1);
if length(unique_condition) == 3
    pred_thres = (psycho(1,2)^(-2)+psycho(2,2)^(-2))^(-1/2);
    set(text(min(xlim),0.4+0.1*(cc+1),sprintf('pred = %g\n',pred_thres)),'color','k');
end

if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/behavior_%s.png',hostname));
end
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
        PSTH_condition_absheading_choice(:,cc,hh,1) = squeeze(nanmean(rate_lip(right_targ_ind,:,select_correct,hh,cc),3));
        PSTH_condition_absheading_choice(:,cc,hh,2) = squeeze(nanmean(rate_lip(left_targ_ind,:,select_correct,hh,cc),3));
        
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,1),'color',colors(unique_condition(cc),:),'linew',2.5);
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,2),'--','color',colors(unique_condition(cc),:),'linew',2.5);
        title(sprintf('abs(heading) = %g, %g trials, correct only',abs(unique_heading(hh)),N_rep));
    end
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    axis tight;
    
    if if_bounded
        plot(xlim,[decis_thres(unique_condition(cc)) decis_thres(unique_condition(cc))],'c--');
        %         ylim([min(ylim),decis_thres*1.1]);
    end
    
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    
end

set(findobj(435,'type','axes'),'ylim',[y_min y_max]);

if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/PSTH_%s.png',hostname));
end

% =====  Delta-PSTH (stim type), heading, correct only ====
set(figure(436),'name','Delta PSTH'); clf;
set(gcf,'uni','norm','pos',[0.326       0.041       0.668       0.681]);
y_max = -inf; y_min = inf;

for hh = 1:length(unique_heading)
    
    figure(436);  subplot(sub_x,sub_y,hh); hold on;
    
    for cc = 1:length(unique_condition)
        plot(ts,PSTH_condition_absheading_choice(:,cc,hh,1)-PSTH_condition_absheading_choice(:,cc,hh,2),'color',colors(unique_condition(cc),:),'linew',2.5);
    end
    title(sprintf('abs(heading) = %g, %g trials, correct only',abs(unique_heading(hh)),N_rep));
    
    plot(t_motion,vel/max(vel)*max(ylim)/5,'k--');
    y_min = min(y_min,min(ylim));
    y_max = max(y_max,max(ylim));
    axis tight;
end
set(findobj(436,'type','axes'),'ylim',[y_min y_max]);

if ION_cluster
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/deltaPSTH_%s.png',hostname));
end
%}

total_time = toc(total_time)

