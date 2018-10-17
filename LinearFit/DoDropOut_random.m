function result = DoDropOut_random(sub_sample_using_ratio, sub_sample_ratio, sub_sample_N)
% Regression 2. Fitting model traces with real data. (We decided to do this when Alex came to Shanghai) @HH20170808
% Add (random) subsampling methods to test the robustness of MSE

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sub_sample_using_ratio = 0; % If 0, using the same number of neurons
% sub_sample_raio = 0.5;
% sub_sample_N = 60;
permN = 1000;
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
% model_PSTH_trace_with_heter_shortTau2000 = load('diff_PSTH_trace_with_heter_shortTau_fixed2000ms.mat');


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MSEs = nan(permN, size(data_to_fit,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if sub_sample_using_ratio
        this_weights = nan(permN,  round(n_origin * sub_sample_ratio));
    else
        this_weights = nan(permN, sub_sample_N);
    end
    
    parfor pp = 1:permN
        %%%%%%%%%%%%%%% Subsampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sub_sample_using_ratio
            sub_sample_index = randperm(n_origin, round(n_origin * sub_sample_ratio));
        else
%             if sub_sample_N <= n_origin % Do actual sub_sample
%                 sub_sample_index = randperm(n_origin, sub_sample_N);
%             else % Do bootstrap for this instead
              sub_sample_index = datasample(1:n_origin, sub_sample_N);
%             end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fit_basis_trace_per_cell =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(fit_time_range),:); % Should be indices in original rate_ts{1}!
        %     test_basis_trace_per_cell =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(test_time_range),:);
        %     valid_basis_trace_per_cell =  data_to_fit_PSTH(sub_sample_index,find_common_ts_in_rate_ts(valid_time_range),:);
        % Early stopping not used here because the test error neves goes up, as confirmed in LinearFit plot function
        
        [fitted_w, MSE] = LinearFit(fit_optimal_trace_interp,fit_basis_trace_per_cell,alpha);
        
        % Save data
        MSEs(pp, ii) = MSE - alpha * norm(fitted_w);
        this_weights(pp,:) = fitted_w;
        
        parfor_progress;
    end
    
    Weights{ii} = this_weights;
    parfor_progress(0);

end

result.MSEs = MSEs;
result.Weights = Weights;

%% Plotting
% -- MSEs --
figure(1807); clf;
set(gcf,'uni','norm','pos',[0.036       0.355       0.919       0.417]);
colors = {'b','r','g','k'};

h1 = subplot(1,3,1);
bar_comp_MSE = BarComparison(result.MSEs,'SEM',0,'Pair',0,'Equal',0,'LineStyle','none','Colors',colors,'axes',h1);
set(gca,'yscale','log');

hs = get(gca,'children');
legend(hs(end:-2:2),data_names,'Location','best');
title('MSE');
ylim([0.1 20])

% Calculate two-way p values using definition
xs = result.MSEs;
ps_mse = Perm_pvalue(xs);

% -- Weights --
h2 = subplot(1,3,2);
Kurtosis_all = nan(permN,size(data_to_fit,1)+1);

% -- Gaussian random --
if ~sub_sample_using_ratio
    Weights_rand = sort(abs(randn(permN, sub_sample_N)),2,'descend');
%   Weights_rand = sort(gamrnd(1,1,[permN, sub_sample_N]),2,'descend');
    Weights_rand = Weights_rand ./ (Weights_rand(:,1) * ones(1,size(Weights_rand,2)));
    SeriesComparison(Weights_rand, 1:size(Weights_rand,2),...
                    'SEM',0,'Colors','m','axes',h2);
    Kurtosis_all(:,1) = kurtosis(Weights_rand,[],2);
end

for ii = 1 : size(data_to_fit,1)
    Weights = sort(result.Weights{ii},2,'descend');
    Weights = Weights ./ (Weights(:,1) * ones(1,size(Weights,2)));
    
    SeriesComparison(Weights, 1:size(Weights,2),...
                 'SEM',0,'Colors',colors{ii},'axes',h2);

    Kurtosis_all(:,ii+1) = kurtosis(result.Weights{ii},[],2);
end

legend off
xlim([0 sub_sample_N]); ylim([0 1])
title('Sorted weights');

% -- Kurtosis --
h3 = subplot(1,3,3);
bar_comp_kurt = BarComparison(Kurtosis_all,'SEM',0,'Pair',0,'Equal',0,'LineStyle','none','Colors',['m',colors],'axes',h3);
title('Kurtosis');
% hold on; plot(xlim,[3 3],'m--','linew',2)
SetFigure(20);

ps_Kurt_within = Perm_pvalue(Kurtosis_all); % Within the groups

% ps_Kurt_rand = nan(1,size(Kurtosis_all,2)); % With randon weights (Gaussian, Kur = 3) % Wrong, should be abs(Gaussian)
% for g1 = 1:size(Kurtosis_all,2)
%     p = 2 * sum(Kurtosis_all(:,g1) > 3)/size(Kurtosis_all,1);
%     if p > 1, p = 2-p; end
%     ps_Kurt_rand(g1) = p;
% end

% -- Saving data --
result.bar_comp_MSE = bar_comp_MSE; 
result.bar_comp_kurt = bar_comp_kurt; 
result.Kurtosis_all = Kurtosis_all; 
result.ps_mse = ps_mse;
result.ps_Kurt = ps_Kurt_within;
% result.ps_Kurt_rand = ps_Kurt_rand;

suffix = sprintf('%g_%g_%g',sub_sample_using_ratio, sub_sample_ratio, sub_sample_N);
save(['./result/result_random_' suffix '.mat'], 'result');
saveas(gcf,['./result/DropOutRandom_' suffix '.fig']);

function ps = Perm_pvalue(xs) 
ps = nan(size(xs,2));

% NOTE: ttest or ttest2 cannot be used because they're not "sampling data"!!
% Instead, use the definition of p-values (and different subsamplings can be pooled)
for g1 = 1:size(xs,2)
    for g2 = g1 + 1:size(xs,2)
        p = 2 * sum(xs(:,g1) > xs(:,g2))/size(xs,1);
        if p > 1, p = 2-p; end
        ps(g1,g2) = p;
        ps(g2,g1) = p;
    end
end

