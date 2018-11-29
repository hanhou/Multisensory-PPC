%% Parameter scan


%% 20170508_g_w_lip_lip_dc_w_lip_lip
% %%
% for g_w_lip_lip = 0:0.5:10
%     for dc_w_lip_lip = 0:-0.2:-g_w_lip_lip/1.5
%         lip_HH({'g_w_lip_lip',g_w_lip_lip;'dc_w_lip_lip',dc_w_lip_lip});
%     end
% end

%% 20170508_int_vis_gAndK_1
% for g_w_int_vis = 1:3:30
%     for k_int_vis = [4 6 8 10 15 30]
%         lip_HH({'g_w_int_vis',g_w_int_vis;'k_int_vis',k_int_vis});
%     end
% end

%% 20170509_bounds_scan_2_good
% n=0;
% group_result = [];
% for thr_single = [28:2:58 inf]
%     for thr_comb = [28:2:58 inf]
%         n=n+1;
%         result = lip_HH({'decis_thres',[thr_single thr_single thr_comb]},{'psycho','RT'});
%
%         group_result(n,:) = [thr_single thr_comb result.psycho(:)']
%
%     end
% end

%% 20170509_bounds_scan_RT
% group_result = [];
% bounds = [28:1:55 inf];
%
% for bb = 1:length(bounds)
%     result = lip_HH({'decis_thres',[bounds(bb) bounds(bb) nan];'unique_stim_type',[1 2];'save_folder','20170509_bounds_scan_RT/'},{'psycho','RT'});
%     group_result(bb,[1 2 3 5 6]) = [bounds(bb) result.psycho(:,2)' squeeze(mean(mean(result.RT,1),2))'];
%
%     result = lip_HH({'decis_thres',[nan nan bounds(bb)];'unique_stim_type',[3];'save_folder','20170509_bounds_scan_RT/'},{'psycho','RT'});
%     group_result(bb,[4 7]) = [result.psycho(:,2)' squeeze(mean(mean(result.RT,1),2))']
% end

%% 20170510_lipTolip_fixedRT
% n=0;
% for g_w_lip_lip = 3:2:10
%     for dc_w_lip_lip = -linspace(g_w_lip_lip/3,g_w_lip_lip-1,10)
%         n=n+1;
%         fprintf('%g, %g, %g\n',n,g_w_lip_lip,dc_w_lip_lip);
%         lip_HH({'g_w_lip_lip',g_w_lip_lip;'dc_w_lip_lip',dc_w_lip_lip;'save_folder','20170510_lipTolip_fixedRT/'});
%     end
% end

%% 20170510_intTolip_fixedRT
% n=0;
% for g_w_lip_int = 5:2:20
%     for dc_w_lip_int= -linspace(g_w_lip_int/3, g_w_lip_int-1,10)
%         n=n+1;
%         fprintf('%g, %g, %g\n',n,g_w_lip_int,dc_w_lip_int);
%         lip_HH({'g_w_lip_int',g_w_lip_int;'dc_w_lip_int',dc_w_lip_int;'save_folder','20170510_intTolip_fixedRT/'});
%     end
% end

%% 20170528_dropout (non-scaled)
% drop_outs = 0.2:0.1:0.9;
% for dd = 1:length(drop_outs)
%     heter_dropout = ones(1,)*drop_outs(dd);
%     lip_HH({'heter_dropout',heter_dropout;'save_folder','20170528_dropout/'});
% end

%% 20170529_dropout (non-scaled) and log normal
% drop_outs = 0:0.1:0.5;
% lognormal_ratio = 0.1:0.3:2;
%
% for dd = 1:length(drop_outs)
%     for ll = 1:length(lognormal_ratio)
%         heter_dropout = ones(1,4)*drop_outs(dd);
%         heter_lognormal = ones(1,4)*lognormal_ratio(ll);
%         lip_HH({'heter_dropout',heter_dropout;'heter_lognormal',heter_lognormal;'save_folder','20170529_dropout_and_lognormal/'});
%     end
% end

%% 20170530_distribution of log normal diagonal weight
%{
lognormal_ratio = 0.1:0.3:2;
figure(1317); clf; hold on;
cc = lines;
delta_theta = 0;
h=[];

for ll = 1:length(lognormal_ratio)
    heter_lognormal = ones(1,4)*lognormal_ratio(ll);
    result = lip_HH({'heter_lognormal',heter_lognormal;'save_folder','20170529_dropout_and_lognormal/'});
    
    figure(1317);
    off_diag = round(delta_theta/(360/size(result.w_lip_lip,1)));
    x{ll} = diag(result.w_lip_lip,off_diag);
    [N,edges] = histcounts(x{ll},min(x{ll}):0.0003:max(x{ll}));
    plot([edges edges(end)+edges(2)-edges(1)],[0 ;smooth(N,5);0]/sum(N)/(edges(2)-edges(1)),'-','linew',2,'color',cc(ll,:));
end

legend(num2str(lognormal_ratio'));
ylabel('Prob density');
xlabel(['LIP --> LIP, \Delta\theta = ' num2str(delta_theta)]);

max_y = max(ylim);
for ll = 1:length(lognormal_ratio)
    plot(mean(x{ll})*ones(1,2),[max_y*1.05 max_y*1.15],'-','linew',2,'color',cc(ll,:))
end
%}

%% 20170530_scaled drop out and log normal with grouped figure
% %{
clear
scan = tic;
% save_folder = '20170602_drop0.6_lognorm1.2/';
% save_folder = '20170609_rescan_with_SVM_0.6_1.8/';
% save_folder = '20170613_SVM_rescaled_and_weight_noise_corr_-1_1.5/';

% Coarse
% drop_outs = 0:0.1:0.8;
% lognormal_ratio = 0:0.3:4;

% Fine
% drop_outs = 0.4:0.1:0.9;
% lognormal_ratio = 0:0.2:1.5;

% Fixed
% drop_outs = 0.6*ones(1,4);
% lognormal_ratio = 1.8*ones(1,5);

% Test
% xs = 0.6 ; x_name = 'drop_outs';
% ys = 1.2; y_name = 'lognormal_ratio';

% log-normal and vis/vest negative noise correlation
% xs = 0:-0.1:-1;  x_name = 'vis_vest_weight_noise_cor';
% ys = 0:0.3:4;  y_name = 'heter_lognormal';

% xs = 0:-0.5:-2;  x_name = 'vis_vest_weight_noise_cor';
% ys = 0:0.5:3;  y_name = 'heter_lognormal';

save_folder = '20171121_dropout_and_lognormal_all_raw_spike_corr/'
xs = 0:0.2:0.8;  x_name = 'heter_dropout';
ys = 0:0.5:3;  y_name = 'heter_lognormal';

initiated = 0;

total_n = length(xs)*length(ys)
errors = {};

for xx = 1:length(xs)
    for yy = 1:length(ys)
        try
            result = lip_HH({x_name, ones(1,4)*xs(xx); y_name, ones(1,4)*ys(yy);...
                'save_folder',save_folder},{'h_grouped'});
            
            finished = (xx-1)*length(ys)+yy;
            fprintf('=== Grouped progress: %g / %g ===\n',finished,total_n);
            fprintf('=== Estimated remaining time: %s ===\n', datestr(toc(scan)/finished*(total_n-finished)/24/3600,'HH:MM'));
            
            h_grouped = result.h_grouped;
            
            % Initiated once
            if ~initiated
                for ff = 1:length(h_grouped)
                    figure(1700+ff); clf;
                    set(gcf,'uni','norm','pos',[0.051       0.068       0.914        0.79]);
                    hs{ff} = tight_subplot(length(xs),length(ys),0.02); % Column-wise
                end
                initiated = 1;
            end
            
            % Copy grouped figure
            for ff = 1:length(h_grouped)
                this_f = figure(1700+ff);
                set(gcf,'name',sprintf('%s, %s',x_name,y_name));
                h_copyed = copyobj(h_grouped(ff),this_f);
                linkprop([hs{ff}(xx+(yy-1)*length(xs)),h_copyed],'pos'); % Put the figure into subfigure
                
                if ~(yy==1 && xx==length(xs))
                    set(h_copyed,'xtick',[],'ytick',[]);
                end
                
                xlabel(h_copyed,''); ylabel(h_copyed,''); title(h_copyed,'');
                
                if yy == 1;  ylabel(h_copyed,num2str(xs(xx))); end
                if xx == 1;  title(h_copyed,num2str(ys(yy))); end
                drawnow;
            end
        catch error
            errors = {errors, error};
        end
    end
end

for ff = 1:length(h_grouped)
    figure(1700+ff);
    
    % Clean up
    delete(hs{ff});
    
    % Adjust
    hhs = findobj(gcf,'type','axes');
    % axis(hhs,'tight');
    linkaxes(hhs,'xy');
    
    %{
        for i=1:length(hhs)
            axes(hhs(i));
            pos=get(gca,'pos');
            set(gca,'pos',[pos(1:2) pos(3)/1.32 pos(4)/1.14])
        end
    %}
    
    % Save figure
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%sGrouped_%g.png',save_folder,ff));
    saveas(gcf,sprintf('./result/%sGrouped_%g.fig',save_folder,ff),'fig');
    
end
disp('Grouped done');
toc(scan);
%}

%% 20170911_scan weights gain of sensory --> int
%{
clear
scan = tic;
% save_folder = '20170602_drop0.6_lognorm1.2/';
% save_folder = '20170609_rescan_with_SVM_0.6_1.8/';
save_folder = '20170914_Scan weights gain_no noise_conflict4/';

% Coarse

xs = 10*10.^(-0.6:0.2:0.6);  x_name = 'g_w_int_vest';
ys = 10*10.^(-0.6:0.2:0.6);  y_name = 'g_w_int_vis';

group_result = [];

initiated = 0;

total_n = length(xs)*length(ys);

for xx = 1:length(xs)
    for yy = 1:length(ys)
        
        result = lip_HH({x_name, xs(xx); y_name, ys(yy);...
            'save_folder',save_folder},{'h_grouped';'psycho'});
        
        group_result = [group_result; xs(xx) ys(yy) result.psycho(:)']

        
        finished = (xx-1)*length(ys)+yy;
        fprintf('=== Grouped progress: %g / %g ===\n',finished,total_n);
        fprintf('=== Estimated remaining time: %s ===\n', datestr(toc(scan)/finished*(total_n-finished)/24/3600,'HH:MM'));
        
        h_grouped = result.h_grouped;
        
        % Initiated once
        if ~initiated
            for ff = 1:length(h_grouped)
                figure(1700+ff); clf;
                set(gcf,'uni','norm','pos',[0.051       0.068       0.914        0.79]);
                hs{ff} = tight_subplot(length(xs),length(ys),0.02); % Column-wise
            end
            initiated = 1;
        end
        
        % Copy grouped figure
        for ff = 1:length(h_grouped)
            this_f = figure(1700+ff);
            set(gcf,'name',sprintf('%s, %s',x_name,y_name));
            h_copyed = copyobj(h_grouped(ff),this_f);
            linkprop([hs{ff}(xx+(yy-1)*length(xs)),h_copyed],'pos'); % Put the figure into subfigure
            
            if ~(yy==1 && xx==length(xs))
                set(h_copyed,'xtick',[],'ytick',[]);
            end
            
            xlabel(h_copyed,''); ylabel(h_copyed,''); title(h_copyed,'');
            
            if yy == 1;  ylabel(h_copyed,num2str(xs(xx))); end
            if xx == 1;  title(h_copyed,num2str(ys(yy))); end
            drawnow;
        end
    end
end

for ff = 1:length(h_grouped)
    figure(1700+ff);
    
    % Clean up
    delete(hs{ff});
    
    % Adjust
    hhs = findobj(gcf,'type','axes');
    axis(hhs,'tight');
    % linkaxes(hhs,'x');
    
%{
        for i=1:length(hhs)
            axes(hhs(i));
            pos=get(gca,'pos');
            set(gca,'pos',[pos(1:2) pos(3)/1.32 pos(4)/1.14])
        end
%}
    
    % Save figure
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%sGrouped_%g.png',save_folder,ff));
    saveas(gcf,sprintf('./result/%sGrouped_%g.fig',save_folder,ff),'fig');
    
end
save(sprintf('./result/%spsyho_scan_weight_gain_no_noise.mat',save_folder),'group_result');
disp('Grouped done');

toc(scan);

%% Plot result for this scan (threshold)
% load psyho_scan_weight_gain_no_noise;
t1 = reshape(group_result(:,end-2),length(ys),length(xs))';
t2 = reshape(group_result(:,end-1),length(ys),length(xs))';
t3 = reshape(group_result(:,end),length(ys),length(xs))';
p = t3.*sqrt(1./t1.^2+1./t2.^2); % Prediction ratio

figure(9131215); clf;
subplot(2,2,1); imagesc(t1); title('Vest threshold'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,2); imagesc(t2); title('Vis threshold'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,3); imagesc(t3); title('Comb threshold'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,4); imagesc(p); title('Comb/Optimal'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

% Save figure
export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%spsycho_scan_threshold.png',save_folder));
saveas(gcf,sprintf('./result/%spsycho_scan_threshold.fig',save_folder),'fig');

%% Plot result for this scan (bias)
% load psyho_scan_weight_gain_no_noise;
t1 = reshape(group_result(:,3),length(ys),length(xs))';
t2 = reshape(group_result(:,4),length(ys),length(xs))';
t3 = reshape(group_result(:,5),length(ys),length(xs))';
p = abs(t3-(t1+t2)/2); % Prediction bias

figure(9131215); clf;
subplot(2,2,1); imagesc(t1); title('Vest bias'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,2); imagesc(t2); title('Vis bias'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,3); imagesc(t3); title('Comb bias'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

subplot(2,2,4); imagesc(p); title('abs(Comb bias - Optimal bias)'); colorbar;
set(gca,'xaxisloc','top','xticklabel','','ytick',1:length(xs),'yticklabel',num2str(xs')); xlabel('vis weight'); ylabel('vest weight');

% Save figure
export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%spsycho_scan_bias.png',save_folder));
saveas(gcf,sprintf('./result/%spsycho_scan_bias.fig',save_folder),'fig');
%}

%% 201701103_Scan time-varying weight factor gamma
%{
clear
scan = tic;
% save_folder = '20171104_Time-varying gamma/';
% save_folder = '20171104_Time-varying gamma_fixed RT 1.00/';
save_folder = '20171109_Time-varying gamma_Max_the_same_Information/';
% Coarse

% xs = [0 1];   x_name = 'gamma'; % Rows
% ys = [5 10 20 40]; y_name = 'runs'; % Columns

xs = 1;   x_name = 'runs'; % Rows
ys = [0 0.5 1 1.5 2]; y_name = 'gamma'; % Columns


group_result = [];

initiated = 0;

total_n = length(xs)*length(ys);
errors = {};

for xx = 1:length(xs)
    for yy = 1:length(ys)
        try
            result = lip_HH({'gamma', ys(yy); 'N_rep',200; ...  % Fewer N_reps, different runs
                'save_folder',save_folder},{'h_grouped';'psycho';'thres_boot'});
            %         result = lip_HH({x_name, xs(xx); 'g_w_int_vest',ys(yy) ; 'g_w_int_vis',ys(yy);...  % Fewer N_reps, different runs
            %             'save_folder',save_folder},{'h_grouped';'psycho';'thres_boot'});
            
            finished = (xx-1)*length(ys) +yy;
            group_result = [group_result; xs(xx) result.psycho(:)']
            group_thres_boot{finished} = result.thres_boot;
            
            fprintf('=== Grouped progress: %g / %g ===\n',finished,total_n);
            fprintf('=== Estimated time: %s ===\n', datestr(toc(scan)/finished*(total_n-finished)/24/3600,'HH:MM'));
            
            h_grouped = result.h_grouped;
            
            % Initiated once
            if ~initiated
                for ff = 1:length(h_grouped)
                    figure(1700+ff); clf;
                    set(gcf,'uni','norm','pos',[0.051       0.068       0.914        0.79]);
                    hs{ff} = tight_subplot(length(xs),length(ys),0.02); % Column-wise
                end
                initiated = 1;
            end
            
            % Copy grouped figure
            for ff = 1:length(h_grouped)
                this_f = figure(1700+ff);
                set(gcf,'name',sprintf('%s, %s',x_name,y_name));
                h_copyed = copyobj(h_grouped(ff),this_f);
                linkprop([hs{ff}(xx+(yy-1)*length(xs)),h_copyed],'pos'); % Put the figure into subfigure
                
                if ~(yy==1 && xx==length(xs))
                    set(h_copyed,'xtick',[],'ytick',[]);
                end
                
                xlabel(h_copyed,''); ylabel(h_copyed,''); title(h_copyed,'');
                
                if yy == 1;  ylabel(h_copyed,num2str(xs(xx))); end
                if xx == 1;  title(h_copyed,num2str(ys(yy))); end
                drawnow;
            end
        catch error
            errors = {errors,error};
        end
    end
end

for ff = 1:length(h_grouped)
    figure(1700+ff);
    
    % Clean up
    delete(hs{ff});
    
    % Adjust
    hhs = findobj(gcf,'type','axes');
    axis(hhs,'tight');
    % linkaxes(hhs,'x');
    
%{
        for i=1:length(hhs)
            axes(hhs(i));
            pos=get(gca,'pos');
            set(gca,'pos',[pos(1:2) pos(3)/1.32 pos(4)/1.14])
        end
%}
    
    % Save figure
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%sGrouped_%g.png',save_folder,ff));
    saveas(gcf,sprintf('./result/%sGrouped_%g.fig',save_folder,ff),'fig');
    
end

save(sprintf('./result/%spsyho_scan_gamma_no_noise.mat',save_folder),'group_result','group_thres_boot');
disp('Grouped done');

toc(scan);

%% Plot result for this scan (threshold)
% load psyho_scan_gamma_no_noise;
% xs = [0 0.1 0.3 0.5 1 2 4 8];
figure(1654); clf; hold on;
for gg = 1:length(xs)
    this_boot = group_thres_boot{gg};
    pred_ratio = this_boot(3,:)./this_boot(4,:);
    mean_pred_ratios(gg) = mean(pred_ratio);
    err_pred_ratios(gg) = std(pred_ratio);
end
errorbar(xs,mean_pred_ratios,err_pred_ratios);

% Save figure
export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%sPred_ratio.png',save_folder));
saveas(gcf,sprintf('./result/%sPred_ratio.fig',save_folder),'fig');
%}

%% Save para_scan file
copyfile([mfilename('fullpath') '.m'],sprintf('./result/%s/',save_folder));


