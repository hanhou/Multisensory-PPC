function result = DoDropOut_seq(knockout_M,knockout_iter_N)
% Regression 2. Fitting model traces with real data. (We decided to do this when Alex came to Shanghai) @HH20170808
% Add (sequential) subsampling methods to test the robustness of MSE
% "Sequential" means that at each time I knock out the top M cells with the largest weights

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_sample_using_ratio = 1; % If 0, using the same number of neurons
sub_sample_ratio = 1;
sub_sample_N = 100;
% knockout_M = 1; % How many cells will be drop out for each iteration
%knockout_iter_N = 20;
permN = 1;
alpha = 1; % Regularization term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ==== Load optimal traces (Fit target) from the trace without heterogeneity ====
% Generated from calling "result = lip_HH({},{'ts','mean_diff_PSTH_correct_allheading'});"
model_mean_PSTH_trace_without_heter = load('mean_trace_without_heter.mat'); % Gamma = 0
% model_mean_PSTH_trace_without_heter = load('mean__trace_without_heter_gamma=1.mat'); % Gamma = 1
% model_mean_PSTH_trace_without_heter = load('mean__trace_without_heter_gamma=8.mat'); % Gamma = 8
%     model_mean_PSTH_trace_without_heter = load('diff_PSTH_trace_without_heter_vestvelocity.mat'); % Vestibular velocity

model_ts = model_mean_PSTH_trace_without_heter.result.ts*1000;
model_PSTH_optimal_M1 = squeeze(model_mean_PSTH_trace_without_heter.result.mean_diff_PSTH_correct_allheading);

% ==== Load basis traces ====
real_diff_PSTH_trace_LIP = load('diff_PSTH_trace_LIP.mat');
%     real_diff_PSTH_trace_MST = load('diff_PSTH_trace_MST_scaling0.78.mat'); % HH20180916  % Scale of sigma
real_diff_PSTH_trace_MST = load('diff_PSTH_trace_MST_scaling1.mat'); % HH20180916

model_PSTH_trace_with_heter_optimal = load('diff_PSTH_trace_with_heter.mat');
% model_PSTH_trace_with_heter_vest10_vis1 = load('diff_PSTH_trace_with_heter_vest10_vis1.mat');

% model_PSTH_trace_with_heter_shortTau = load('diff_PSTH_trace_with_heter_shortTau.mat')
model_PSTH_trace_with_heter_shortTau = load('diff_PSTH_trace_with_heter_shortTau_fixed100ms.mat'); % 20181010 Fixed real dynamics with 100ms tau
model_PSTH_trace_with_heter_shortTau2000 = load('diff_PSTH_trace_with_heter_shortTau_fixed2000ms.mat');


rate_ts = -495 : 10 : 2195; % LIP time stamp

data_to_fit = {%Times  %Data   %Name
    rate_ts, real_diff_PSTH_trace_LIP.tmp; % 'Real LIP data';
    model_ts, model_PSTH_trace_with_heter_optimal.a; % 'Optimal model with heterogeneity';
    model_ts, model_PSTH_trace_with_heter_shortTau.a; %  'ShortTau with heterogeneity';
%     model_ts, model_PSTH_trace_with_heter_shortTau.a; %  'ShortTau with heterogeneity';
    real_diff_PSTH_trace_MST.PSTH_ts_scaled, permute(real_diff_PSTH_trace_MST.group_PSTH_prefMnull_smoothed,[3 1 2]); %'Real MST data';
    };
% data_names = {'LIP','M2','M3_{100ms}','M3_{2000ms}','MST'};
data_names = {'LIP','M2','M3(100ms)','MST'};
colors = {'b','r','g','k'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSEs = nan(permN, knockout_iter_N, size(data_to_fit,1));
% knockout_ind = nan(permN, knockout_iter_N, knockout_M, size(data_to_fit,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1243); clf;
set(gcf,'uni','norm','pos',[0.023       0.268       0.957       0.523]);
h = tight_subplot(1,4,0.02,[0.1 0.15],[0.07 0.25]);
h5 = axes('pos',[0.805 0.149 0.171 0.634]);

for ii = 1 : size(data_to_fit,1)
    
    use_data_to_fit = ii % 1: LIP, 2: M2, 3: M3-100ms, 4: MST
    
    data_to_fit_time = data_to_fit{use_data_to_fit,1};
    data_to_fit_PSTH = data_to_fit{use_data_to_fit,2};
    
    % Times
    find_common_ts_in_rate_ts = find(min(model_ts)<=data_to_fit_time & data_to_fit_time<=max(model_ts));
    common_ts = data_to_fit_time(find_common_ts_in_rate_ts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alpha = 10;
    % fit_time_range = 500 < common_ts & common_ts <= 1200;
    % fit_time_range = rand(1,length(common_ts)) < 0.2;
    
    % Interleaved training and testing
    fit_time_range = false(1,length(common_ts));
    test_time_range = fit_time_range;
    valid_time_range = fit_time_range;
    
    fit_time_range(1:3:end) = true;
    test_time_range(2:3:end) = true;
    valid_time_range(3:3:end) = true; % Use other (time points) to test
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fit_ts = common_ts(fit_time_range);
    test_ts = common_ts(test_time_range);
    valid_ts = common_ts(valid_time_range);
    
    % Optimal model data
    fit_optimal_trace_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),fit_ts);
    test_optimal_trace_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),test_ts);
    valid_optimal_trace_interp = interp1(model_ts,model_PSTH_optimal_M1(:,:),valid_ts);
    
    n_origin = size(data_to_fit{use_data_to_fit,2},1);
    
    parfor_progress(permN);

    % Heter data (basis functions)
%     if sub_sample_using_ratio
%         this_weights = nan(permN, knockout_iter_N, round(n_origin * sub_sample_ratio));
%     else
%         this_weights = nan(permN, knockout_iter_N, sub_sample_N);
%     end
    
    parfor pp = 1:permN
        %%%%%%%%%%%%%%% Subsampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sub_sample_using_ratio
            this_weights{pp} = nan( knockout_iter_N, round(n_origin * sub_sample_ratio));
            sub_sample_index = randperm(n_origin, round(n_origin * sub_sample_ratio));
        else
            this_weights{pp} = nan( knockout_iter_N, sub_sample_N);
            
%             if sub_sample_N <= n_origin % Do actual sub_sample
%                 sub_sample_index = randperm(n_origin, sub_sample_N);
%             else % Do bootstrap for this instead
                sub_sample_index = datasample(1:n_origin, sub_sample_N);
%             end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fit_basis_trace_per_cell{pp} =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(fit_time_range),:); % Should be indices in original rate_ts{1}!
        %     test_basis_trace_per_cell =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(test_time_range),:);
        %     valid_basis_trace_per_cell =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(valid_time_range),:);
        % Early stopping not used here because the test error neves goes up, as confirmed in LinearFit plot function
        
        basis_this = fit_basis_trace_per_cell{pp};
        ind_this = 1:size(basis_this,1);
        
        for kk = 1 : knockout_iter_N % Sequential knockout
            % Fit
            [fitted_w, MSE] = LinearFit(fit_optimal_trace_interp,basis_this,alpha);
            
            % Save
            MSEs(pp, kk, ii) = MSE - alpha * norm(fitted_w);
            this_weights{pp}(kk, ind_this) = fitted_w;
            
            % Sort
            [~, ind] = sort(fitted_w,'descend');
            
            % Knockout
            to_knockout_this = ind(1:knockout_M);
            knockout_ind{pp}(kk, :, ii) = ind_this(to_knockout_this);
            
            basis_this(to_knockout_this,:,:) = [];
            ind_this(to_knockout_this) = [];
            
        end
        
        parfor_progress;
    end
    
    Weights{ii} = this_weights;
    parfor_progress(0);

    % -- Some plotting --
    figure(1243); % Weights of the first permutation and MSEs of all permutation
    knockout_Ns = ((1:knockout_iter_N)-1) * knockout_M;
    
    axes(h5);
    errorbar(knockout_Ns, mean(MSEs(:,:,ii),1),std(MSEs(:,:,ii),[],1),'linew',2,'color',colors{ii}); hold on
    set(gca,'yscale','log'); axis tight; xlabel('Number of dropped out'); title('MSE')
    
    % pos = [2 3 5 6];
    norm_weights = this_weights{1}; % Plot the last permutation
    norm_weights = norm_weights ./(max(norm_weights,[],2)*ones(1,size(norm_weights,2)));
    axes(h(ii));
    %subplot(1,5,ii); 
    imagesc(knockout_Ns,[1 size(norm_weights,2)], norm_weights') ; colormap(gray)
    % xlabel('Number of dropped out'); ylabel('Cell number'); 
    axis xy
    
    if ii == 1
        ylabel('Cell number');
    else
        set(gca,'ytick',[])
    end
    
    if ~isempty(knockout_ind{1})
        hold on; plot(knockout_Ns, squeeze(knockout_ind{1}(:,:,ii)), 'xr','linew',2,'markersize',7)
    end
    
    title(data_names{ii})    % axis equal
    SetFigure(20);
    
    
    set(figure(1323 + ii),'name',data_names{ii}); clf % Order of best weights cells
    n_Topcells = 7; n_Dropout = min(5,knockout_iter_N);
    set(gcf,'uni','norm','pos',[0.073       0.077        0.88       0.829/5*n_Dropout]);
    
    for dd = 1:n_Dropout  % The (dd-1)th dropout
        
        % Sort according to weights
        % [sort_norm_weights,ind] = sort(norm_weights(dd,:),2,'descend');
        
        % Sort according to weights*dynamic range
        dynamic_range = range(range(fit_basis_trace_per_cell{1},2),3);
        contribution = dynamic_range'.* norm_weights(dd,:);
        contribution = contribution/max(contribution);
        [sort_norm_weights,ind] = sort(contribution,2,'descend');
        
        
        to_plot = ind(~isnan(sort_norm_weights));
        to_plot_w = sort_norm_weights(~isnan(sort_norm_weights));
        
        for tt = 1:n_Topcells % The (tt)th cell
            subplot(n_Dropout, n_Topcells, tt + (dd -1)*n_Topcells)
            try
            plot(fit_ts, squeeze(fit_basis_trace_per_cell{1}(to_plot(tt),:,:)))  
            catch
                keyboard
            end
            title(sprintf('#%g, %g',to_plot(tt), to_plot_w(tt)));
            axis tight
            hold on; plot(xlim,[0 0],'k--'); plot([0 0], ylim,'k--')
        end
    end
    
end

result.to_knockout_ori = knockout_ind;
result.MSEs = MSEs;
result.Weights = Weights;

suffix = sprintf('_%g_%g_%g_drop%g_%g',sub_sample_using_ratio, sub_sample_ratio, sub_sample_N,knockout_M,knockout_iter_N);
save(['./result/result_Seq' suffix '.mat'], 'result');
saveas(1243,['./result/DropOutSeq_Weights_MSE' suffix '.fig']);

for ii = 1: size(data_to_fit,1)
    saveas(1323 + ii,['./result/DropOutSeq_Topcells' data_names{ii} suffix '.fig']);
end
