% Using LOLE and BCFisher to calculate linear Fisher information in heterogenerous
% MST population with differential correlation.
%
% 1. Simulate heterogeneous MST populations with f'f'T noise
% 2. Compute Fisher info in three ways
%    2.1  Ground truth
%    2.2  Use LOLE early stopping method (-ds and ds)
%    2.3  Use BCFisher (empirical f' and SIGMA plus bias correction)
%    2.4  My previous SVM decoder (use whole psychometric curve)
%
% First by Han Hou @ 20180428,  houhan@gmail.com
clear;
rng('default')
rng('shuffle');
addpath(genpath('/home/hh/Multisensory-PPC/util'));

% ====== For cluster running ======
hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name

if ~isempty(strfind(hostname,'node')) || ~isempty(strfind(hostname,'clc')) % contains(hostname,{'node','clc'})  % ION cluster;
    if isempty(gcp('nocreate'))
        parpool;
        %         parpool(hostname(1:strfind(hostname,'.')-1),20);
    end
    ION_cluster = 1;
else % My own computer
    if strcmp(version('-release'),'2013a')
        if matlabpool('size') == 0
            matlabpool;
        end
    else
        if isempty(gcp('nocreate'))
            parpool;
        end
    end
    ION_cluster = 0;
end


% ==== Scan nu_neurons ===
% %{
    epsis = 0.0027;
    mean_baselines = 20;    
    std_baselines = 20;
    nu_neurons = round(10.^linspace(1,log10(1000),10)) + 1;
%}


% ==== Scan epsis and baseline ===
%{
    epsis = [eps 10.^(linspace(-4,0,10))];
    std_baselines = 20;
    mean_baselines = linspace(eps,30,10);
    nu_neurons = 500;
%}


% ==== Scan baseline and var_baseline ===
%{
    epsis = 0.0027;
    std_baselines = linspace(eps,30,10);
    mean_baselines = linspace(0,30,10);
    nu_neurons = 500;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% According to my test (mainScanCorrelatedPoisson.m), with the average firing
% rate f = 25 Hz, dt = 5e=3 and rho = 0.4 should lead to rho_max = 0.2, cor90 =
% 0.1, and firing_increase = 30%, fano factor = 0.85. 
% 
% Note that, however, when the baseline is changed, this approximation may not
% always stand...

dt = 7e-3; % s
% rho = 0.4;  % =>> rho_real = 0.2, cor90 = 0.1
rho = 0.1;  % =>> rho_real = 0.1, cor90 = 0.06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_max = 2; % s
time = 0:dt:time_max;

nu_reps_LOLE = 4000;  % 1000 is enough for BCFisher, but 3000+ is needed for LOLE early stopping.
heter = 1;

% ======== Motion profile =======
% Use the same setting as experiments
NumOfSigma = 3.5;
sigma_t = time_max/2/NumOfSigma;
speed_t = exp(-(time-time_max/2).^2/(2*sigma_t^2));
speed_t = speed_t-min(speed_t);
speed_t = speed_t/max(speed_t);

position_t = cumsum(speed_t);
position_t = position_t-min(position_t);
position_t = position_t/max(position_t);

% ====================== Parameter END ======================== %
% ============================================================= %

% if ~ION_cluster; progressbar(0); end
count = 0;
count_all = length(epsis) * length(mean_baselines) * length(std_baselines) * length(nu_neurons);
optimality_matrix = [];

ifplot_corr = length(nu_neurons) == 1;

for ee = 1:length(epsis) % Scan epsi
    
    epsilon_0 = epsis(ee); % Epsilon(t) = epsilon_0 * 1/v(t).

    for mbmb = 1:length(mean_baselines) % Scan baseline
        for sbsb = 1:length(std_baselines) % Scan std_baseline
            
            I_optimal = []; I_optimal_var = [];
            I_ilPPC = []; I_ilPPC_var = [];
            
            for nn = 1:length(nu_neurons)  % Scan nu_neurons
                
                count = count + 1;
                
                % ======== Neurons =========
                nu_neuron = nu_neurons(nn);
                
                % ---- Uniformly distrubuted preferred heading ----
                theta_pref = linspace(-pi,pi,nu_neuron)';
                
                % ========= Spatial-temporal tunings =======
                % Using gamma distribution to introduce heterogeneity
                gamma_para = [ % Mean   Std  (Desired)
                    45,    20;  % amp
                    1.4,   0.3;  % k (controls sigma of spatial tuning)
                    0.7,   0.5; % beta
                    %     eps, eps;
                    mean_baselines(mbmb),  std_baselines(sbsb);  % original 20, 5
                    ];
                
                if ~ heter
                    gamma_para(:,2) = eps; % Turn off heterogeneity for test
                end
                
                gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
                gamma_B = (gamma_para(:,1))./gamma_A;    gamma_B(isnan(gamma_B)) = 0;
                
                tuning_paras = nan(nu_neuron,size(gamma_para,1));
                for pp = 1:size(gamma_para,1)
                    tuning_paras(:,pp) = gamrnd(gamma_A(pp),gamma_B(pp),nu_neuron,1);
                end
                
                Ai = tuning_paras(:,1);
                ki = tuning_paras(:,2);
                betai = tuning_paras(:,3);
                Bi = tuning_paras(:,4);
                
                betai(betai > 1) = 1; % beta should be smaller than 1
                
                % To make sure the min peak is larger than (1 + min_signal_baseline_ratio)*Baseline
                min_signal_baseline_ratio = 0.5;
                to_trunc = find(Bi > Ai./(min_signal_baseline_ratio + betai));
                Bi (to_trunc) = Ai(to_trunc)./(min_signal_baseline_ratio + betai(to_trunc));
                
                % [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
                %  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline
                tuning_theta = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
                    (1 - betai * speed_t ) .* repmat(Bi,1,length(time));
                % [f']: Derivative of tuning
                tuning_theta_der = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)).* ...
                    (ki .* sin(theta_pref - theta)) * speed_t;
                
                % -- Correlations
                k = 2;  fano = 1;  % Exp. Decay correlation
                C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)
                
                % ========= Population responses =========
                %     theta_stim = [-sqrt(epsilon_0) sqrt(epsilon_0)];
                theta_stim = [-3/180*pi 3/180*pi];
                
                % 1). For ground truth (at 0)
                f_at_0_ts = tuning_theta(mean(theta_stim));
                f_der_at_0_ts = tuning_theta_der(mean(theta_stim));
                
                % 2). For spike count (two stimuli aroud 0)
                f_S1_ts = tuning_theta(theta_stim(1));
                f_der_S1_ts = tuning_theta_der(theta_stim(1));
                
                f_S2_ts = tuning_theta(theta_stim(2));
                f_der_S2_ts = tuning_theta_der(theta_stim(2));
                
                % ========= Initialization ========
                %     SIGMA_0_ilPPC = zeros(nu_neuron,nu_neuron);
                SIGMA_epsi_ilPPC = zeros(nu_neuron,nu_neuron);
                
                %      I_0_ts = nan(1,length(time));
                I_epsi_ts = nan(1,length(time));
                
                % Spiking matrix
                
                spike_S1 = nan(nu_neuron,length(time),nu_reps_LOLE);
                spike_S2 = nan(nu_neuron,length(time),nu_reps_LOLE);
                
                % ================== Generate Poisson spike train ====================
                
                % epsilon_0 = 0;
                fprintf('[%g/%g: epsi = %.2g, Base = (%.2g, %.2g), NN = %g]\n   Generating Poisson spike trains...',...
                            count, count_all, epsis(ee), mean_baselines(mbmb), std_baselines(sbsb), nu_neurons(nn));
                        
                parfor tt = 1:length(time)
                    
                    f_at_0 = f_at_0_ts(:,tt);
                    f_der_at_0 = f_der_at_0_ts(:,tt);
                    
                    f_S1 = f_S1_ts(:,tt);   f_der_S1 = f_der_S1_ts(:,tt);
                    f_S2 = f_S2_ts(:,tt);   f_der_S2 = f_der_S2_ts(:,tt);
                    
                    % ======= Time-dependent epsilon (after 20180426 Talk, see OneNote) ====
                    epsi_this = epsilon_0 / (speed_t(tt) + 1e-6);
                    
                    % ---- Correlated firing rate and Poisson spikes
                    % C_epsi = C + epsi_this * (f_der_at_0 * f_der_at_0') ./ sqrt(f_at_0 * f_at_0');
                    % w_epsi = sqrtm(C_epsi);
                    
                    SIGMA_S1 = C .* sqrt(f_S1 * f_S1') + epsi_this * (f_der_S1 * f_der_S1');
                    SIGMA_S2 = C .* sqrt(f_S2 * f_S2') + epsi_this * (f_der_S2 * f_der_S2');
                    
                    % This is my equation: p_i = dt*f_i + sqrt(dt) * sum_over_j (sqrtm(SIGMA - fI)_ij * eta_j)
                    prob_S1 = dt * f_S1 * ones(1,nu_reps_LOLE) + ...
                        sqrt(dt) * real(sqrtm(SIGMA_S1 - diag(f_S1))) * normrnd(0,1,nu_neuron,nu_reps_LOLE);
                    prob_S2 = dt * f_S2 * ones(1,nu_reps_LOLE) + ...
                        sqrt(dt) * real(sqrtm(SIGMA_S2 - diag(f_S2))) * normrnd(0,1,nu_neuron,nu_reps_LOLE);
                    
                    spike_S1(:,tt,:) = rand(size(prob_S1)) < prob_S1 ; % Poisson
                    spike_S2(:,tt,:) = rand(size(prob_S2)) < prob_S2 ; % Poisson
                end
                fprintf('Done\n');
                
                % ========== Verify corrcoef ==========
                t_bin = 0.5;
                t_centers = t_bin/2:t_bin:time_max-t_bin/2;
                
                % Mean corrcoef. Follow Beck, 2008
                aux_mask = ones(nu_neuron,1) * (1:nu_neuron);
                mask = ((1-cos(abs(aux_mask-aux_mask')/nu_neuron*2*pi))/2)<.5;
                mask = mask.*(1-eye(nu_neuron));
                
                if ifplot_corr
                    figure(1323); clf; hold on
                end
                
                for tt = 1:length(t_centers)
                    corr_matrix = corrcoef(squeeze(sum(spike_S1(:, ...
                        round((t_centers(tt) - t_bin/2)/dt)+1 : round((t_centers(tt) + t_bin/2)/dt),...
                        :),2))');
                    corr_this = circAverage(corr_matrix);
                    corr_max(tt) = corr_this(find(corr_this==1)+1);
                    corr_this(corr_this==1) = nan;
                    
                    if ifplot_corr, plot(theta_pref/pi*180, corr_this), end
                    
                    mean_corr_less_than_90_actual(tt) = nanmean(corr_matrix(logical(mask)));
                end
                
                % Mean firing increase
                real_firing_rate = mean(spike_S1,3) / dt;
                aver_firing_increase = mean(mean(real_firing_rate - f_S1_ts))/mean(mean(f_S1_ts));
                desired_aver_firing = mean(mean(f_S1_ts));
                
                fprintf('   ==> cor90 ~ %.3g, max_cor ~ %.3g, FR increase ~ %.2g%% (desired %.2g Hz)\n', ...
                    mean(mean_corr_less_than_90_actual), mean(corr_max), aver_firing_increase*100, desired_aver_firing);
                
                % Plotting
                if ifplot_corr
                    legend(num2str(t_centers'));
                    ylabel('corrcoef')
                    xlabel('\Deltapref')
                    title(sprintf('\\rho = %g, \\deltat = %g\ncor90 \\approx %.3g, max\\_cor ~ %.3g, FR increase \\approx %.2g%% ', ...
                        rho, dt, mean(mean_corr_less_than_90_actual), mean(corr_max), aver_firing_increase*100));
                end
                
                %% ========== Calculate information using LOLE or BCFisher ===========
                %     keyboard;
                
                % ========= 1. Use counting windows =======
                t_info_bin = 0.100; % s
                t_info_step = t_info_bin; % Must be non-overlapping!!
                t_info_centers = t_info_bin/2 : t_info_step : time_max-t_info_bin/2;
                
                BCFisher_this = nan(1,length(t_info_centers));
                varFI = nan(1,length(t_info_centers));
                EarlyStopping_FIVAL = nan(1,length(t_info_centers));
                actual_popul_act = nan(nu_neurons(nn),length(t_info_centers));
                desired_popul_act = actual_popul_act;
                
                fprintf('   Caculate BCFisher:  ')
                
                for TT = 1:length(t_info_centers)
                    
                    %         fprintf('   Info of time bin %g / %g:',TT,length(t_info_centers));
                    
                    t_info_bin_this = floor((t_info_centers(TT) - t_info_bin/2)/dt) + 1 : floor((t_info_centers(TT) + t_info_bin/2)/dt);
                    
                    % ====== Get spike counts =======
                    count_S1 = squeeze(sum(spike_S1(:,t_info_bin_this,:),2)); % nu_neurons * nu_reps
                    count_S2 = squeeze(sum(spike_S2(:,t_info_bin_this,:),2)); % nu_neurons * nu_reps
                    
                    % >>>>>>  1. Compute Fisher info by BCFisher
                    %         fprintf('  BCFisher...')
                    [BCFisher_this(TT), varFI(TT), ~ ] = BCFisher(count_S1', count_S2', diff(theta_stim));
                    %         fprintf('Done... EarlyStopping...')
                    
                    %         fprintf('Skipped...')
                    % >>>>>>  2. Compute Fisher Info by EarlyStopping
                    %         [EarlyStopping_FIVAL(TT), EarlyStopping_FITR(TT)] = EarlyStopping(count_S1', count_S2',diff(theta_stim),1/3,1/3,20,1);
                    
                    fprintf('>')
                    
                    % --- Save for verifying population activity ---
                    actual_popul_act(:,TT) = mean(count_S1,2)/t_info_bin;
                    desired_popul_act(:,TT) = mean(f_S1_ts(:, t_info_bin_this),2);
                    
                end
                
                % ======= 1. Sum of all momentary info (standard of optimality) ========
                I_optimal(nn) = sum(BCFisher_this) ;  % Straight sum of info
                I_optimal_var(nn) = sum(varFI) ; % Var of sum = sum of Var (independent info bins)
                
                % ======= 2. Info of straight sum of spikes (ilPPC) ========
                count_S1_all_trial = squeeze(sum(spike_S1(:,:,:),2)); % nu_neurons * nu_reps
                count_S2_all_trial = squeeze(sum(spike_S2(:,:,:),2)); % nu_neurons * nu_reps
                
                [I_ilPPC(nn), I_ilPPC_var(nn), ~] = BCFisher(count_S1_all_trial', count_S2_all_trial', diff(theta_stim));
                
                fprintf('\n   I_ilPPC = %g, I_optimal = %g, opt. = %.2g%%\nFinished %3.1f%%\n\n',...
                    I_ilPPC(nn),I_optimal(nn),I_ilPPC(nn)/I_optimal(nn)*100,(count/count_all*100));
                
            end
            
            optimality_matrix(ee,mbmb,sbsb) = mean(I_ilPPC./I_optimal)*100; % optimality matrix does not take care of nu_neurons
            
        end
    end
end

optimality = I_ilPPC./I_optimal ;
optimality_std = sqrt(I_ilPPC_var./(I_optimal.^2) + I_optimal_var.* (I_ilPPC.^2) ./ (I_optimal).^4);


%% ====== Plot momentary information (only for the last set of parameters) ========
% %{
figure(1614);  clf;
errorbar(t_info_centers,BCFisher_this,sqrt(varFI)); hold on;
plot(t_info_centers, EarlyStopping_FIVAL);

if ION_cluster
    file_name = sprintf('./1_Info_time');
    %         export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end

%}

% ====== Compare ilPPC with optimal =======
color = 'b';
if_clf = 1;

figure(11);
if if_clf, clf, end
% Optimal
errorbar(nu_neurons, I_optimal, sqrt(I_optimal_var),[color 'o-'],'linew',3,'markerfacecol',color);
hold on;

% ilPPC
errorbar(nu_neurons + 10, I_ilPPC, sqrt(I_ilPPC_var),[color 's--'],'linew',3,'markerfacecol',color);
thres_discrim_epsi = sqrt(2/I_ilPPC(end))/pi*180;
set(text(max(ylim)*0.5,I_ilPPC(end)*0.8,sprintf('\\sigma = %g',thres_discrim_epsi)),'color',color);

xlabel('Number of neurons');
title('BCFisher');
axis tight

if ION_cluster
    file_name = sprintf('./2_TotalInfo');
    %         export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end


% ====== Optimality of ilPPC =========
figure(12); % Recovery of info
if if_clf, clf, end
errorbar(nu_neurons, optimality * 100, optimality_std * 100, [color 'o-'],'linew',3,'markerfacecol',color); hold on;
axis tight; plot(xlim,[100 100],'k--','linew',2);
ylim([50 105]);
xlabel('Number of neurons');
ylabel('% of info recovery by ilPPC');
title('BCFisher');


if ION_cluster
    file_name = sprintf('./3_Optimality_ilPPC_BCFisher');
    %         export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end


%% ========= Plot optimality matrix =========
optimality_matrix = squeeze(optimality_matrix);

if length(mean_baselines)>1 && length(epsis)>1
    figure(13); clf;
    epsi_log_x = [2*log10(epsis(2))-log10(epsis(3)) log10(epsis(2:end))]; % Add one entry for epsi = 0
    imagesc(epsi_log_x,mean_baselines,optimality_matrix'); axis xy;
    hold on;
    contour(epsi_log_x,mean_baselines, optimality_matrix','color','k','linew',1.5,'ShowText','on');
    
    epsi_log_x_tick = [min(epsi_log_x) linspace(round(epsi_log_x(2)),round(max(epsi_log_x)),5)];
    epsi_log_x_label = [-inf epsi_log_x_tick(2:end)];
    set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)
    set(gca,'ytick',0:5:max(mean_baselines))
    colormap(hot); caxis([50 100])
    xlabel('log_{10}(\epsilon_0)');
    ylabel('Mean(Baseline)');
    
    if ION_cluster
        saveas(13,'./4_Optimality_matrix_BCFisher');
    end
    
    
elseif length(mean_baselines)>1 && length(std_baselines)>1
    figure(13); clf;
    imagesc(std_baselines,mean_baselines,optimality_matrix); axis xy;
    hold on;
    contour(std_baselines,mean_baselines, optimality_matrix,'color','k','linew',1.5,'ShowText','on');
    
    %     set(gca,'ytick',0:5:max(mean_baselines))
    colormap(hot); caxis([50 100])
    xlabel('Std(Baseline)');
    ylabel('Mean(Baseline)');
    
    if ION_cluster
        saveas(13,'./4_Optimality_matrix_BCFisher');
    end
    
end




% ======= Other tests =======
%{
% --- Population activity ---

figure(1513);
h1 = subplot(1,2,1);
plot(actual_popul_act);   title('Actual population activity');
h2 = subplot(1,2,2);
plot(desired_popul_act);   title('Desired population activity');
linkaxes([h1 h2],'xy')

if ION_cluster
    saveas(gcf,'./4_Population Activity');
end

%}

%{
% Speed profile
figure(3); clf; hold on;
plot(time,speed_t,'k','linew',2);
set(gca,'colororder',color_order);
for tt = 1:length(example_time)
    plot(time(example_time(tt)),speed_t(example_time(tt))+0.05,'v','linew',2,'markersize',10);
end
set(gca,'color',[0.8 0.8 0.8])

%}

%% ==== Draw Distribution
% mean_b = 20; std_b = 15;
% b = gamrnd((mean_b./std_b)^2,mean_b/(mean_b./std_b)^2,1,500);
% figure(427); clf; hist(b,[0:1:60]); xlim([0 60]); SetFigure()
% xlabel('Baseline'); ylabel('Num of cells')



