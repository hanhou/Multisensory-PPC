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
n=0;
for g_w_lip_lip = 3:2:10
    for dc_w_lip_lip = -linspace(g_w_lip_lip/3,g_w_lip_lip-1,10)
        n=n+1;
        fprintf('%g, %g, %g\n',n,g_w_lip_lip,dc_w_lip_lip);
        lip_HH({'g_w_lip_lip',g_w_lip_lip;'dc_w_lip_lip',dc_w_lip_lip;'save_folder','20170510_lipTolip_fixedRT/'});
    end
end
    
%% 20170510_intTolip_fixedRT
% n=0;
% for g_w_lip_int = 5:2:20
%     for dc_w_lip_int= -linspace(g_w_lip_int/3, g_w_lip_int-1,10)
%         n=n+1;
%         fprintf('%g, %g, %g\n',n,g_w_lip_int,dc_w_lip_int);
%         lip_HH({'g_w_lip_int',g_w_lip_int;'dc_w_lip_int',dc_w_lip_int;'save_folder','20170510_intTolip_fixedRT/'});
%     end
% end


