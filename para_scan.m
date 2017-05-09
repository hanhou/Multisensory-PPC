%% Parameter scan

% %%
% n=0
% for g_w_lip_lip = 0:0.5:10
%     for dc_w_lip_lip = 0:-0.2:-g_w_lip_lip/1.5
%         n=n+1;
%     end
% end
% n
%%
% %%
% for g_w_lip_lip = 0:0.5:10
%     for dc_w_lip_lip = 0:-0.2:-g_w_lip_lip/1.5
%         lip_HH({'g_w_lip_lip',g_w_lip_lip;'dc_w_lip_lip',dc_w_lip_lip});
%     end
% end

%%
% for g_w_int_vis = 1:3:30
%     for k_int_vis = [4 6 8 10 15 30]
%         lip_HH({'g_w_int_vis',g_w_int_vis;'k_int_vis',k_int_vis});
%     end
% end
n=0;
group_result = [];
for thr_single = [18:5:60 inf]
    for thr_comb = [18:5:60 inf]
        n=n+1;
        result = lip_HH({'decis_thres',[thr_single thr_single thr_comb]},{'psycho'});
        
        group_result(n,:) = [thr_single thr_comb result.psycho(:)']
       
    end
end
