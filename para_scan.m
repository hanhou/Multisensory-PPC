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
clear
scan = tic;
% save_folder = '20170602_drop0.6_lognorm1.2/';
save_folder = '20170609_rescan_with_SVM_0.6_1.8/';

% Coarse
% drop_outs = 0:0.1:0.8;
% lognormal_ratio = 0:0.3:4;

% Fine
% drop_outs = 0.4:0.1:0.9;
% lognormal_ratio = 0:0.2:1.5;

% Fixed
drop_outs = 0.6*ones(1,4);
lognormal_ratio = 1.8*ones(1,5);

% Test
% drop_outs = 0.6;
% lognormal_ratio = 1.2;


initiated = 0;

total_n = length(drop_outs)*length(lognormal_ratio);

for dd = 1:length(drop_outs)
    for ll = 1:length(lognormal_ratio)
        heter_dropout = ones(1,4)*drop_outs(dd);
        heter_lognormal = ones(1,4)*lognormal_ratio(ll);
        result = lip_HH({'heter_dropout',heter_dropout;'heter_lognormal',heter_lognormal;...
            'save_folder',save_folder},{'h_grouped'});
        
        finished = (dd-1)*length(lognormal_ratio)+ll;
        fprintf('=== Grouped progress: %g / %g ===\n',finished,total_n);
        fprintf('=== Estimated remaining time: %s ===\n', datestr(toc(scan)/finished*(total_n-finished)/24/3600,'HH:MM'));
        
        h_grouped = result.h_grouped;
        
        % Initiated once
        if ~initiated
            for ff = 1:length(h_grouped)
                figure(1700+ff); clf; 
                set(gcf,'uni','norm','pos',[0.051       0.068       0.914        0.79]);
                hs{ff} = tight_subplot(length(drop_outs),length(lognormal_ratio),0.02); % Column-wise
            end
            initiated = 1;
        end
        
        % Copy grouped figure
        for ff = 1:length(h_grouped)
            this_f = figure(1700+ff);
            h_copyed = copyobj(h_grouped(ff),this_f);
            linkprop([hs{ff}(dd+(ll-1)*length(drop_outs)),h_copyed],'pos'); % Put the figure into subfigure
            
            if ~(ll==1 && dd==length(drop_outs))
                set(h_copyed,'xtick',[],'ytick',[]);
            end
            
            xlabel(h_copyed,''); ylabel(h_copyed,''); title(h_copyed,'');
            
            if ll == 1;  ylabel(h_copyed,['dropout = ' num2str(drop_outs(dd))]); end
            if dd == 1;  title(h_copyed,['lognorm = ' num2str(lognormal_ratio(ll))]); end
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
    linkaxes(hhs,'x');
%     
%     for i=1:length(hhs)
%         axes(hhs(i));
%         pos=get(gca,'pos');
%         set(gca,'pos',[pos(1:2) pos(3)*1.32 pos(4)*1.14])
%     end
    
    % Save figure
    export_fig('-painters','-nocrop','-m2' ,sprintf('./result/%sGrouped_%g.png',save_folder,ff));
    saveas(gcf,sprintf('./result/%sGrouped_%g.fig',save_folder,ff),'fig');
    
end
disp('Grouped done');

toc(scan);