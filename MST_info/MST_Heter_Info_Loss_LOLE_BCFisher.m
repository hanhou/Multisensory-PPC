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
rng('shuffle');

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

% ======= Only scan nu_neurons =====
epsilon_0 = 0.0027; % Epsilon(t) = epsilon_0 * 1/v(t).

% nu_neurons = round(10.^linspace(1,log10(1000),10)) + 1;
% nu_neurons = round(10.^linspace(1,log10(500),10)) + 1;
nu_neurons = 100;

time_max = 2000; % ms
dt = 2; % ms
time = 0:dt:time_max;
nu_reps_LOLE = 1000; 

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

% ======== Tuning Parameters =======
% Using gamma distribution to introduce heterogeneity
gamma_para = [ % Mean   Std  (Desired)
    45,    20;  % amp
    1.4,   0.3;  % k (controls sigma of spatial tuning)
    0.7,   0.5; % beta
    20,  15;  % original 20, 5
    ];

% gamma_para(:,2) = eps; % Turn off heterogeneity for test

gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
gamma_B = (gamma_para(:,1))./gamma_A;    gamma_B(isnan(gamma_B)) = 0;

% ====================== Parameter END ======================== %
% ============================================================= %

if ~ION_cluster ;progressbar(0); end
count = 0;

for nn = 1:length(nu_neurons)  % Scan nu_neurons
    count = count + 1;
    
    % ======== Neurons =========
    nu_neuron = nu_neurons(nn);
    
    % ---- Uniformly distrubuted preferred heading ----
    theta_pref = linspace(-pi,pi,nu_neuron)';
    
    % ========= Spatial-temporal tunings =======
    tuning_paras = nan(nu_neuron,size(gamma_para,1));
    for pp = 1:size(gamma_para,1)
        tuning_paras(:,pp) = gamrnd(gamma_A(pp),gamma_B(pp),nu_neuron,1);
    end
    
    tuning_paras(tuning_paras(:,3)>1,3) = 1; % beta should be smaller than 1
    
    % [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
    %  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline
    tuning_theta = @(theta) tuning_paras(:,1).* exp(tuning_paras(:,2).*(cos(theta_pref - theta)-1)) * (speed_t) + ...
        (1 - tuning_paras(:,3) * speed_t ) .* repmat(tuning_paras(:,4),1,length(time));
    % [f']: Derivative of tuning
    tuning_theta_der = @(theta) tuning_paras(:,1).* exp(tuning_paras(:,2).*(cos(theta_pref - theta)-1)).* ...
        (tuning_paras(:,2) .* sin(theta_pref - theta)) * speed_t;
    
    % -- Correlations 
    rho = 0.5;  k = 2;  fano = 1;  % Decay correlation
    C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)

    % ========= Population responses ========= 
    theta_stim = [-sqrt(epsilon_0) sqrt(epsilon_0)];
    
    % 1). For ground truth (at 0)
    resp_at_0 = tuning_theta(mean(theta_stim));
    resp_der_at_0 = tuning_theta_der(mean(theta_stim));
    
    % 2). For spike count (two stimuli aroud 0)
    resp_S1 = tuning_theta(theta_stim(1));
    resp_S2 = tuning_theta(theta_stim(2));
    
    % ========= Initialization ========
    %     SIGMA_0_ilPPC = zeros(nu_neuron,nu_neuron);
    SIGMA_epsi_ilPPC = zeros(nu_neuron,nu_neuron);
    
    %      I_0_ts = nan(1,length(time));
    I_epsi_ts = nan(1,length(time));
    
    % Spiking matrix
    
    spike_S1 = nan(nu_neuron,length(time),nu_reps_LOLE);
    spike_S2 = nan(nu_neuron,length(time),nu_reps_LOLE);
      
    % ================== Generate Poisson spike train ==================== 
    
    epsilon_0 = 0;
    for tt = 1:length(time)   
        %     for tt = round(time_max/2/dt)  % Only one time point
        f = resp_at_0(:,tt);
        f_der = resp_der_at_0(:,tt);
        
        % ======= Time-dependent epsilon (after 20180426 Talk, see OneNote) ====
        epsi_this = epsilon_0 / (speed_t(tt) + 1e-6);
        
        % ---- Correlated firing rate and Poisson spikes        
        % See my notes in OneNote 20180501
        
        % SIGMA_ij = C_ij * sqrt(f_i*f_j) + epsilon * f'_i * f'_j
        SIGMA = diag(sqrt(f)) * C * diag(sqrt(f)) + epsi_this * (f_der * f_der');  
        
        % W_ij = sqrtm (SIGMA - fI)
        W = real(sqrtm(SIGMA - diag(f)));
        
        % p = dt * f_i + sum_j (W_ij * eta_j) * sqrt (dt);
        prob_S1 =  (dt/1000) * resp_S1(:,tt) * ones(1,nu_reps_LOLE) + W * normrnd(0,1,nu_neuron,nu_reps_LOLE) * sqrt((dt/1000));
        spike_S1(:,tt,:) = rand(size(prob_S1)) < prob_S1 ; % Poisson
        
        prob_S2 =  (dt/1000) * resp_S2(:,tt) * ones(1,nu_reps_LOLE) + W * normrnd(0,1,nu_neuron,nu_reps_LOLE) * sqrt((dt/1000));
        spike_S2(:,tt,:) = rand(size(prob_S2)) < prob_S2 ; % Poisson
    end
    
    %%{
    % Validation of "Fano matrix": Cov(r_i,r_j)/sqrt(f_i,f_j)
    time_range = round(time_max * 3/7 / dt) : round(time_max * 4/7 / dt);
    temp_activity = squeeze(sum(spike_S1(:,time_range,:),2));
    cov_matrix = cov(temp_activity');
    fano_matrix = cov_matrix ./ sqrt(mean(temp_activity,2).* mean(temp_activity,2)');
    
    figure(1742); clf; subplot(1,2,1); imagesc(fano_matrix); colorbar(); title('Fano matrix');
    subplot(1,2,2); plot(fano_matrix(round(nu_neuron/2),:)); title(sprintf('rho = %g, epsi_0 = %g, dt = %g ms', rho, epsilon_0, dt));
    set(gcf,'uni','norm','pos',[0.126       0.464       0.561       0.269]);
    
    figure(1743); clf; subplot(1,2,1); imagesc(mean(spike_S1(:,:,:),3)/(dt/1000)); colorbar; title('Rate matrix');
    subplot(1,2,2); 
    cor_matrix = corrcoef(temp_activity');
    imagesc(cor_matrix); colorbar();
    
    % Mean corrcoef, Follow Beck, 2008
    aux_mask = ones(nu_neuron,1) * [1:nu_neuron];
    mask = ((1-cos(abs(aux_mask-aux_mask')/nu_neuron*2*pi))/2)<.5;
    mask = mask.*(1-eye(nu_neuron));
    
    mean_corr_less_than_90 = nanmean(cor_matrix(logical(mask)));
    title(['mean corr = ' num2str(mean_corr_less_than_90)]);
    
    %}
    
    
    % ========== Calculate information for larger windows ===========
    keyboard;
    
    for TT = 1:length(TT)
            % ======= Add correlation ======
        % Exponentially decay correlation (non-differential in the heterogeneous case)
        % SIGMA_0_ts = fano * C .* sqrt(f * f');
        
        % ====== Calculate Momentary Fisher Information ==============        
        % Total covariance matrix = SIGMA_0 + Differential Correlation
        epsi_this = epsilon_0 / (speed_t(tt) + 1e-6);
        
        SIGMA_epsi_this_t = fano * C .* sqrt(f * f') + epsi_this * (f_der * f_der');
        
        % >>>>>>  1. Compute Linear Fisher info directly by groud truth I = f'T SIGMA^-1 f' 
        if sum(f) == 0
            % I_0_ts(tt) = 0;
            I_epsi_ts(tt) = 0;
        else
            % I_0_ts(tt) = f_der' * inv(SIGMA_0_ts) * f_der;
            I_epsi_ts(tt) = f_der' * (SIGMA_epsi_this_t \ f_der) ;   % A\B is faster than inv(A)*B
            % I_epsi_ts(tt) = I_0_ts(tt) / (1 + epsi_this * I_0_ts(tt));
        end
        
        % --- Adding momentary covariance matrix for ilPPC (assuming independent noise over time) for ground truth calculation of ilPPC
        % SIGMA_0_ilPPC = SIGMA_0_ilPPC + SIGMA_0_ts;
        SIGMA_epsi_ilPPC = SIGMA_epsi_ilPPC + SIGMA_epsi_this_t;
    
    end
    
    % >>>>>>>>>>>>>>> Now compare ilPPC with optimal information <<<<<<<<<<<<<<<<<
    
    % ======= 1. Sum of all momentary info (standard of optimality) ========
    % I_0_optimal(nn) = sum(I_0_ts) * dt/1000;  % Area under momentary info (so that total info does not dependent on dt)
    I_epsi_optimal(nn) = sum(I_epsi_ts) * dt/1000;
    
    %     % For checking momentary info temporally
    %         I_0_optimal(nn) = I_0_ts(round(time_max/2/dt));
    %         I_epsi_optimal(nn) = I_epsi_ts(round(time_max/2/dt));
    
    % ======= 2. Info of straight sum of spikes (ilPPC) ========
    f_der_ilPPC = sum(resp_der_at_0,2);
    
    % I_0_ilPPC(nn) = f_der_ilPPC' * inv(SIGMA_0_ilPPC) * f_der_ilPPC * dt/1000;
    I_epsi_ilPPC(nn) = f_der_ilPPC' * inv(SIGMA_epsi_ilPPC) * f_der_ilPPC * dt/1000;
    
    
    
    if ~ION_cluster
        progressbar(count/(length(nu_neurons)));
    else
        fprintf('Finished %3.1f%%\n',(count/(length(nu_neurons))*100));
    end   
end


%% ====== Plot momentary information ========
%%{
%figure(10); plot(time,I_0_ts,'r','linew',2)
hold on; plot(time,I_epsi_ts,'b','linew',2);
hold on; plot(time,speed_t*max(I_epsi_ts),'k--');
%}

% ====== Compare ilPPC with optimal =======
figure(11);
% Optimal
% plot(nu_neurons, I_0_optimal,'ro-','linew',3,'markerfacecol','r');
hold on; plot(nu_neurons, I_epsi_optimal,'bo-','linew',3,'markerfacecol','b');

% ilPPC
% plot(nu_neurons, I_0_ilPPC,'ro--','linew',3,'markerfacecol','r');
plot(nu_neurons, I_epsi_ilPPC,'bo--','linew',3,'markerfacecol','b');
thres_discrim_epsi = sqrt(2/I_epsi_ilPPC(end))/pi*180;
set(text(max(ylim)*0.9,I_epsi_ilPPC(end)*0.8,sprintf('\\sigma = %g',thres_discrim_epsi)),'color','b');

xlabel('Number of neurons');
% plot(xlim,1/epsi*[1 1],'k','linew',2);
axis tight

% ====== Optimality of ilPPC =========
figure(12); % Recovery of info
% plot(nu_neurons,I_0_ilPPC./I_0_optimal*100,'ro-','linew',3,'markerfacecol','r'); hold on
plot(nu_neurons,I_epsi_ilPPC./I_epsi_optimal*100,'bo-','linew',3,'markerfacecol','b'); hold on
axis tight; plot(xlim,[100 100],'k--','linew',2);
ylim([50 105]);
xlabel('Number of neurons');
ylabel('% of info recovery by ilPPC');

if ION_cluster
    file_name = sprintf('./0_Optimality 1D');
    export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end



%%
% ========== Tuning properties =========
% --- Population activity
% figure(1); clf; imagesc(tuning_theta(0)); colorbar

example_time = round((150:200:time_max-200)/dt)+1;
thetas_stim = linspace(-pi,0,10); thetas_stim = [thetas_stim fliplr(-thetas_stim(2:end-1))];

% --- Get all tuning curves. Because we have much more theta than t, I let thetas_stim be the vectorized term
all_tunings = nan(nu_neuron,length(thetas_stim),length(example_time));
for tt = 1:length(example_time)
    all_tunings(:,:,tt) = (tuning_paras(:,1)*ones(1,length(thetas_stim))).* exp( (tuning_paras(:,2)*ones(1,length(thetas_stim))) ...
        .*(cos(theta_pref(:) * ones(1,length(thetas_stim)) - ones(length(theta_pref),1) * thetas_stim)-1)) ...
        * speed_t(example_time(tt)) + repmat((1 - tuning_paras(:,3) * speed_t(example_time(tt))) .* tuning_paras(:,4),1,length(thetas_stim));
end

% --- Plot example spatial-temporal tuning curve ---
n_example = 21;
example_neuron = randperm(nu_neuron,n_example);

figure(2); clf;
color_order = colormap(jet);
color_order = color_order(round(linspace(1,size(color_order,1),size(example_time,2))),:);

set(gcf,'uni','norm','pos',[0.038       0.171       0.934       0.571]);
for i = 1:n_example
    subplot(3,7,i); hold on;
    set(gca,'colororder',color_order);
    plot(thetas_stim/pi*180, squeeze(all_tunings(example_neuron(i),:,:)),'linew',2);
    title(sprintf('%0.1f/',tuning_paras(example_neuron(i),:)));
    axis tight;
    xlim([-180 180])
    set(gca,'xtick',-180:90:180); set(gca,'color',[0.8 0.8 0.8])
end

ylims = cell2mat(get(findall(gcf,'type','axes'),'Ylim'));
ylnew = max(ylims(:,2))*1.1;
set(findall(gcf,'type','axes'),'Ylim',[0 ylnew]);
legend(num2str(time(example_time)'));

if ION_cluster
    file_name = sprintf('./2_spatial-temporal tuning');
    export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end


% --- Averaged tuning curve with pref aligned to 0 ---
n_average = 195;
average_neuron = randperm(nu_neuron,n_average);

all_tunings_align = nan(length(average_neuron),length(thetas_stim),length(example_time));
for i = 1:length(average_neuron)
    to_shift = round((0-theta_pref(average_neuron(i)))/(thetas_stim(2)-thetas_stim(1)));
    all_tunings_align(i,:,:) = circshift(squeeze(all_tunings(average_neuron(i),:,:)),to_shift);
end

mean_tuning_align = squeeze(mean(all_tunings_align,1));
se_tuning_align = squeeze(std(all_tunings_align,[],1)/sqrt(length(average_neuron)));

figure(3); clf;   hold on;  set(gca,'colororder',color_order);
errorbar(repmat(thetas_stim'/pi*180,1,size(mean_tuning_align,2)),mean_tuning_align,se_tuning_align,'linew',2);
set(gca,'xtick',-180:90:180); ylim([0 60])
title(['n=' num2str(length(average_neuron))]);
set(gca,'color',[0.8 0.8 0.8])

if ION_cluster
    file_name = sprintf('./3_average tuning');
    export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end


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



