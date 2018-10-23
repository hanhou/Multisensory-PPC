function result = InfoLossGroundTruth(paraOverride)

hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name
if [strfind(hostname,'node') strfind(hostname,'clc')] % contains(hostname,{'node','clc'})  % ION cluster;
    ION_cluster = 1;
else % My own computer
    ION_cluster = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%% Defult Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Tunings
nu_neuron = 10000;

% mean_baselineDropRatioOfAi = 0.7;  % If not scan, this is real mean_beta. If scan, this will be the "baseline dropping ratio of Ai". 20180628

ReLU = 0;  % Whether using ReLUs for each neuron to remove their minimal tuning curve. HH20180627
           % Add ReLU = 2: remove to B; ReLU = 1: remove to (1-beta)*B
plotTuning = 0;

% Using gamma distribution to introduce heterogeneity
coherence = 100; % Simulate coherence change by linearly scaling A and beta (manually!!!) HH20180811

A_mean = 50 * coherence / 100;   A_std = 30 * coherence / 100;
fwhm_mean = 125;  fwhm_std = 50; % Full width at half maximum
beta_mean = 0.6 * coherence / 100;   beta_std = 0.4 * coherence / 100;
B_mean = 20;    B_std = 20;

localSharpen = 2; % According to Gu 2010 Supplementary Fig. 3

% For keepStdMeanRatio
A_mean_default = A_mean;   A_std_default = A_std;
fwhm_mean_default = fwhm_mean;  fwhm_std_default = fwhm_std;
beta_mean_default = beta_mean;   beta_std_default = beta_std;
B_mean_default = B_mean;    B_std_default = B_std;


% 2. Correlations
rho = 0.1;  k = 2;  fano = 1;  % Decay correlation
epsilon_0 = 0.0015; % f'f'T  Corresponds to threshold ~= 4 degree. Original 0.0027

% 3. Miscs
NumOfSigma = 3.5;
delay = 0.140;
theta_stim = 0; % Calculate Info around 0 degree
keepStdMeanRatio = 0; % Automatically keep std/mean ratio the same as the defaults (HH20180712)

% --- Motion profile ---
time_max = 2; % s
% if ION_cluster
dt = 100e-3; % s. This could be large because we're not simulating Poisson spikes.
% else
%     dt = 100e-3;
% end
time = 0:dt:time_max;

% Use the same setting as experiments
sigma_t = time_max/2/NumOfSigma;
speed_t = exp(-(time-time_max/2-delay).^2/(2*sigma_t^2));
speed_t = speed_t-min(speed_t);
speed_t = speed_t/max(speed_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Override Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    len = size(paraOverride,2);
    for ll = 1:len/2
        if exist(paraOverride{1,ll*2-1},'var')
            if ~isnan(paraOverride{1,ll*2})
                eval([paraOverride{1,ll*2-1} '= paraOverride{1,ll*2};']);
%             fprintf('Overriding %s = %s...\n',paraOverride{1,ll*2-1},num2str(paraOverride{1,ll*2}));
            else
                result.I_0_optimal = nan;
                result.I_epsi_optimal = nan;
                result.I_0_ilPPC = nan;
                result.I_epsi_ilPPC = nan;
                result.I_epsi_Saturation =  nan;
                result.psycho_Saturation = nan; 
                result.psycho_ilPPC = nan;
                return
            end
        else
%             fprintf('Parameter ''%s'' not found...\n',paraOverride{1,ll*2-1});
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if keepStdMeanRatio
    A_std = A_std_default * A_mean / A_mean_default;
    fwhm_std = fwhm_std_default * fwhm_mean / fwhm_mean_default;
    beta_std = beta_std_default * beta_mean / beta_mean_default;
    B_std = B_std_default * B_mean / B_mean_default;
end

% Finalize tuning parameters
gamma_para = [ % Mean   Std  (Desired)
    A_mean,  A_std;  % amp
    fwhm_mean,   fwhm_std;  % Controls sigma of spatial tuning: k = log(0.5)/(cos(degtorad(desired_half_width))-1)
    beta_mean,   beta_std; % beta
    B_mean,   B_std;  % original 20, 5
    ];

gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
gamma_B = (gamma_para(:,1))./gamma_A;    gamma_B(isnan(gamma_B)) = 0;

% ---- Uniformly distrubuted preferred heading ----
%                 theta_pref = linspace(-pi,pi,nu_neuron)';

% ---- Bias to +/- 90 degrees (Gu 2006) ---- [Little changing of % recovery]
theta_pref = [randn(1,round(nu_neuron/2))*(1.2*pi/4)-(pi/2) randn(1,nu_neuron-round(nu_neuron/2))*(1.2*pi/4)+(pi/2)]';
theta_pref = mod(theta_pref,2*pi)-pi;
theta_pref = sort(theta_pref);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Do Simulation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('>');

tuning_paras = nan(nu_neuron,size(gamma_para,1));
for pp = 1:size(gamma_para,1)
    tuning_paras(:,pp) = gamrnd(gamma_A(pp),gamma_B(pp),nu_neuron,1);
end

Ai = tuning_paras(:,1);
ki = log(0.5)./(cos(deg2rad(tuning_paras(:,2))/2)-1);
betai = tuning_paras(:,3);
Bi = tuning_paras(:,4);

% Introduce correlation between Ai and Bi (stupid method)
%{
                    [~ ,sortAiId] = sort(Ai);
                    sortBi = sort(Bi);
                    newBi = sort(Bi);
                    
                    slidingWinWid = round(nu_neuron/5);
                    slidingWinStep = round(nu_neuron/10);
                    slidingBeg = 1;
                    while slidingBeg + slidingWinWid - 1 <= nu_neuron
                        randpermThis = randperm(slidingWinWid);
                        newBi(slidingBeg : slidingBeg + slidingWinWid - 1) = sortBi(slidingBeg + randpermThis - 1);
                        slidingBeg = slidingBeg + slidingWinStep;
                    end
                    Bi(sortAiId) = newBi;
                    
                    r = corr(Ai,Bi,'type','Spearman');
                    figure(5); clf; plot(Ai,Bi,'o'); title(sprintf('r = %g',r))
                    drawnow
%}

betai(betai > 1) = 1; % beta should be smaller than 1

% To make sure the min peak is larger than (1 + min_signal_baseline_ratio)*Baseline
% min_signal_baseline_ratio = 0;
% to_trunc = find(Bi > Ai./(min_signal_baseline_ratio + betai));
% Bi (to_trunc) = Ai(to_trunc)./(min_signal_baseline_ratio + betai(to_trunc));

% To make sure the min peak is larger than Baseline + min_signal_baseline_ratio * Amp. HH20180811
min_signal_baseline_ratio = 0.1;
to_trunc = find(Bi > (1 - min_signal_baseline_ratio) * Ai./ betai);
Bi (to_trunc) = (1 - min_signal_baseline_ratio) * Ai(to_trunc)./ betai(to_trunc);


% ========= Spatial-temporal tunings =======
% [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
%  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline
tuning_theta = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
    (1 - betai * speed_t ) .* repmat(Bi,1,length(time));

% Override tuning baseline (Assuming we use ReLUs for each neuron to remove their minimal tuning curve) HH20180627
% %{
if ReLU == 1 % Remove to (1-beta)*B
    tuning_theta = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
        (1 - betai * speed_t ) .* repmat(Bi,1,length(time)) ...
        - repmat((Ai .* exp(-2*ki) < betai.* Bi).* (Ai.*exp(-2*ki) * max(speed_t) + (1 - betai * max(speed_t)) .* Bi) ...  % If normal baseline Ai*exp(-2ki) < betai*Bi
        + (Ai .* exp(-2*ki) >= betai.* Bi).* Bi, ...   % If the original baseline moves upwards
        1,length(time));  
elseif ReLU == 2 % Remove to B (cut off some of the null tuning)
    tuning_theta = @(theta) max(eps, Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
        (1 - betai * speed_t ) .* repmat(Bi,1,length(time)) ...
        - repmat(Bi, 1, length(time)));
end
%}

% [f']: Derivative of tuning (independent of ReLU)
tuning_theta_der = @(theta) localSharpen * Ai.* exp(ki.*(cos(theta_pref - theta)-1)).* ...
    (ki .* sin(theta_pref - theta)) * speed_t;
% -- Plot f'--
% fp = tuning_theta_prime(0);
% figure(); plot(theta_pref/pi*180,fp(:,1000))

% ========== Adding Correlations =========
resp = tuning_theta(theta_stim);
resp_der = tuning_theta_der(theta_stim);

C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)

% SIGMA_0_ilPPC = zeros(nu_neuron,nu_neuron);
SIGMA_epsi_ilPPC = zeros(nu_neuron,nu_neuron);

I_0_ts = nan(1,length(time));
I_epsi_ts = nan(1,length(time));

for tt = 1:length(time)    % Across all trial
    %     for tt = round(time_max/2/dt)  % Only one time point
    f = resp(:,tt);
    f_der = resp_der(:,tt);
    
    % ======= Time-dependent epsilon (after 20180426 Talk, see OneNote) ====
    epsi_this = epsilon_0 / (speed_t(tt) + eps);
    
    % ======= Add correlation ======
    % Exponentially decay correlation (non-differential in the heterogeneous case)
    SIGMA_0_ts = fano * C .* sqrt(f * f');
    
    % Total covariance matrix = SIGMA_0 + Differential Correlation
    SIGMA_epsi_ts = SIGMA_0_ts + epsi_this * (f_der * f_der');
    
    % ====== Compute Linear Fisher info directly by groud truth I = f'T SIGMA^-1 f' =====
    if norm(f_der) < eps
        I_0_ts(tt) = 0;
        I_epsi_ts(tt) = 0;
    else
        I_0_ts(tt) = f_der' * (SIGMA_0_ts \ f_der); % (A\b) is faster than inv(A)*b
        
        % I_epsi_ts(tt) = f_der' * inv(SIGMA_epsi_ts) * f_der;   % Could be replaced by the fact that I_epsi = I_0 / (1 + eta * I_0)
        I_epsi_ts(tt) = I_0_ts(tt) / (1 + epsi_this * I_0_ts(tt));
    end
    
    % ====== Adding momentary covariance matrix for ilPPC (assuming independent noise over time)
    % SIGMA_0_ilPPC = SIGMA_0_ilPPC + SIGMA_0_ts;
    SIGMA_epsi_ilPPC = SIGMA_epsi_ilPPC + SIGMA_epsi_ts;
end

% ======= 1. Sum of all momentary info (standard of optimality) ========
% result.I_0_optimal = sum(I_0_ts) * dt;  % Area under momentary info (so that total info does not dependent on dt)
result.I_0_optimal = nan;
result.I_epsi_optimal = sum(I_epsi_ts) * dt;

%     % For checking momentary info temporally
%         I_0_optimal(nn) = I_0_ts(round(time_max/2/dt));
%         I_epsi_optimal(nn) = I_epsi_ts(round(time_max/2/dt));

% ======= 2. Info of straight sum of spikes (ilPPC) ========
f_der_ilPPC = sum(resp_der,2);

% result.I_0_ilPPC = f_der_ilPPC' * (SIGMA_0_ilPPC \ f_der_ilPPC) * dt;
result.I_0_ilPPC = nan;
result.I_epsi_ilPPC = f_der_ilPPC' * (SIGMA_epsi_ilPPC \ f_der_ilPPC) * dt;
result.epsi_optimality = result.I_epsi_ilPPC / result.I_epsi_optimal;

% Threshold
% ========= Prediction of psychothreshold ========
result.I_epsi_Saturation =  sum(speed_t) ./ epsilon_0 * dt; % Total epsi info = sum(momentary epsi info)*dt = sum(1/momentary_epsi)*dt = sum(speed_t/epsis)*dt
result.psycho_Saturation = rad2deg(sqrt(1./result.I_epsi_Saturation) * sqrt(2)); % Fine discrimination threshold (84% correct) = sqrt(2) * sqrt(1/I_fisher)
result.psycho_ilPPC = rad2deg(sqrt(1./result.I_epsi_ilPPC) * sqrt(2));

% Other outputs
result.I_0_ts = I_0_ts;
result.I_epsi_ts = I_epsi_ts;
result.time = time;

%% Plot tuning
if plotTuning
    
    % %{
    % --- Population activity
    % figure(1); clf; imagesc(tuning_theta(0)); colorbar
    
    example_time = round((0.15:0.2:time_max-0.2)/dt)+1;
    thetas_stim = -pi:pi/8:0; thetas_stim = [thetas_stim fliplr(-thetas_stim(2:end-1))];
    
    % --- Get all tuning curves. Because we have much more theta than t, I let thetas_stim be the vectorized term
    all_tunings = nan(nu_neuron,length(thetas_stim),length(example_time));
    for tt = 1:length(example_time)
        all_tunings(:,:,tt) = (Ai*ones(1,length(thetas_stim))).* exp( (ki*ones(1,length(thetas_stim))) ...
            .*(cos(theta_pref(:) * ones(1,length(thetas_stim)) - ones(length(theta_pref),1) * thetas_stim)-1)) ...
            * speed_t(example_time(tt)) + repmat((1 - betai * speed_t(example_time(tt))) .* Bi,1,length(thetas_stim));
        
        if ReLU == 1  % Move to (1-beta)*B
            % Override tuning baseline (Assuming we use ReLUs for each neuron to remove their minimal tuning curve) HH20180627
            %%{
            % Note that for some cells whose Ai * exp(-2*ki) > betai * Bi, the minimal baseline to remove should be exactly Bi
            all_tunings(:,:,tt) = all_tunings(:,:,tt) - repmat(...
                (Ai.*exp(-2*ki) < betai.*Bi) .* (Ai .* exp(-2 * ki) * max(speed_t) + (1 - betai * max(speed_t)) .* Bi) + ...
                (Ai.*exp(-2*ki) >= betai.*Bi) .* Bi, ...
                1, length(thetas_stim));
            %}
        elseif ReLU == 2  % Move to B
            all_tunings(:,:,tt) = max(eps, all_tunings(:,:,tt) - repmat(Bi, 1, length(thetas_stim)));
        end
    end
    
    
    % --- Plot example spatial-temporal tuning curve ---
    n_example = 45;
    example_neuron = randperm(nu_neuron,n_example);
    
    figure(2); clf;
    color_order = colormap(jet);
    color_order = color_order(round(linspace(1,size(color_order,1),size(example_time,2))),:);
    
    set(gcf,'uni','norm','pos',[0.013       0.065       0.962       0.841]);
    
    for i = 1:n_example
        
        to_shift = round((0 - theta_pref(example_neuron(i))) / (thetas_stim(2)-thetas_stim(1)));
        this_tuning_align = circshift(squeeze(all_tunings(example_neuron(i),:,:)),to_shift);
        
        subplot(5,9,i); hold on;
        set(gca,'colororder',color_order);
        plot([thetas_stim/pi*180 180], [this_tuning_align; this_tuning_align(1,:)],'linew',2);
        % title(sprintf('%0.1f/',tuning_paras(example_neuron(i),:)));
        title(sprintf('%.1f/%.1f/%.1f/%.1f',Ai(example_neuron(i)),ki(example_neuron(i)),betai(example_neuron(i)),Bi(example_neuron(i))));
        axis tight; ylim([0 max(ylim)*1.1])
        xlim([-180 180])
        set(gca,'xtick',-180:90:180); set(gca,'color',[0.8 0.8 0.8])
        
    end
    
    %{
ylims = cell2mat(get(findall(gcf,'type','axes'),'Ylim'));
ylnew = max(ylims(:,2))*1.1;
set(findall(gcf,'type','axes'),'Ylim',[0 ylnew]);
    %}
    legend(num2str(time(example_time)'));
    
    if ION_cluster
        file_name = sprintf('./2_spatial-temporal tuning');
        % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end
    
    
    % --- Averaged tuning curve with pref aligned to 0 ---
    n_average = 195;
    % n_average = nu_neuron;
    average_neuron = randperm(nu_neuron,n_average);
    
    all_tunings_align = nan(length(average_neuron),length(thetas_stim),length(example_time));
    for i = 1:length(average_neuron)
        to_shift = round((0-theta_pref(average_neuron(i)))/(thetas_stim(2)-thetas_stim(1)));
        all_tunings_align(i,:,:) = circshift(squeeze(all_tunings(average_neuron(i),:,:)),to_shift);
    end
    
    mean_tuning_align = squeeze(mean(all_tunings_align,1));
    se_tuning_align = squeeze(std(all_tunings_align,[],1)/sqrt(length(average_neuron)));
    
    figure(3); clf;
    set(gcf,'uni','norm','pos',[0.092       0.229       0.765       0.526]);
    
    subplot(2,4,[1 2 5 6]);
    hold on;  set(gca,'colororder',color_order);
    errorbar(repmat([thetas_stim pi]'/pi*180,1,size(mean_tuning_align,2)),...
        [mean_tuning_align; mean_tuning_align(1,:)], [se_tuning_align; se_tuning_align(1,:)],'linew',2);
    set(gca,'xtick',-180:90:180); ylim([0 60])
    title(['n=' num2str(length(average_neuron))]);
    set(gca,'color',[0.8 0.8 0.8])
     
    subplot(2,4,3);
    hist(Ai,30); title('Amp');
    
    subplot(2,4,4);
    half_widths = acos(1+log(0.5)./ki)*180/pi*2;
    hist(half_widths, 30); title('Tuning half width');
    
    subplot(2,4,7);
    hist(Bi,30); title('Baseline');
    
    subplot(2,4,8);
    hist(betai,30); title('beta');
    xlim([0 1])
    
    if ION_cluster
        file_name = sprintf('./3_average tuning');
        % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end

    figure(2045); clf
    plot(theta_pref(round(linspace(1,end,1000)))/pi*180,half_widths(round(linspace(1,end,1000))),'ok')
    set(gca,'xtick',-180:90:180,'ytick',45:45:270)
    ylim([45 280])

    
    if ION_cluster
        file_name = sprintf('./3_Tuningwidth_PreferDirection');
        % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end
    %}
end
