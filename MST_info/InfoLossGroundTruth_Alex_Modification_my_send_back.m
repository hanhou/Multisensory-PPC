clear all

%%%%%%%%%%%%%%%%%%%%%%%% Defult Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 0;
count_max = 12;
base_exp = 2;

coherence = 0.1;

% 1. Using gamma distribution to introduce heterogeneity
A_mean = 47;   A_std = 30;

% k_mean = 1.7;  k_std = 0.9*1.5;
halfwidth_mean = 110;  halfwidth_std = 40;

beta_mean = 0.7;   beta_std = 0.3;
B_mean = 20*1;    B_std = 20;

localSharpening = 1; % HH

% 2. Correlations
rho = 0.1;  k = 2;  fano = 0.8;  % Decay correlation
epsilon_0 = 0.00139; % f'f'T  Corresponds to threshold ~= 4 degree. Original 0.0027


for count=1:count_max
    nu_neuron = base_exp^(count+1)

    % mean_baselineDropRatioOfAi = 0.7;  % If not scan, this is real mean_beta. If scan, this will be the "baseline dropping ratio of Ai". 20180628

    ReLU = 0;  % Whether using ReLUs for each neuron to remove their minimal tuning curve. HH20180627
    % Add ReLU = 2: remove to B; ReLU = 1: remove to (1-beta)*B
    plotTuning = 1;

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


    % Finalize tuning parameters
    gamma_para = [ % Mean   Std  (Desired)
        A_mean,  A_std;  % amp
        halfwidth_mean,   halfwidth_std;  % Controls sigma of spatial tuning: k = log(0.5)/(cos(degtorad(desired_half_width))-1)
        beta_mean,   beta_std; % beta
        B_mean,   B_std;  % original 20, 5
        ];

    gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
    gamma_B = (gamma_para(:,1))./gamma_A; gamma_B(isnan(gamma_B)) = 0;

    % ---- Uniformly distrubuted preferred heading ----
    %                 theta_pref = linspace(-pi,pi,nu_neuron)';

    % ---- Bias to �� 90 degrees (Gu 2006) ---- [Little changing of % recovery]
    theta_pref = [randn(1,round(nu_neuron/2))*(1.2*pi/4)-(pi/2) randn(1,nu_neuron-round(nu_neuron/2))*(1.2*pi/4)+(pi/2)]';
    theta_pref = mod(theta_pref,2*pi)-pi;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Do Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tuning_paras = nan(nu_neuron,size(gamma_para,1));
    for pp = 1:size(gamma_para,1)
        tuning_paras(:,pp) = gamrnd(gamma_A(pp),gamma_B(pp),nu_neuron,1);
    end

    %truncates the B distribution to 50 and rescales the B above 50
    max_B = max(tuning_paras(:,4));
    tuning_paras(:,4) = (tuning_paras(:,4)<=50).*tuning_paras(:,4) + ...
    (tuning_paras(:,4)>50).*(tuning_paras(:,4)/max_B)*50;

    Ai = tuning_paras(:,1);
    
    % ki = tuning_paras(:,2);
    ki = log(0.5)./(cos(deg2rad(tuning_paras(:,2))/2)-1);

    betai = tuning_paras(:,3);
    Bi = tuning_paras(:,4);



    betai(betai > 1) = 1; % beta should be smaller than 1

    % To make sure the min peak is larger than (1 + min_signal_baseline_ratio)*Baseline
    min_signal_baseline_ratio = 0.5;
    to_trunc = find(Bi > Ai./(min_signal_baseline_ratio + betai));
    Bi (to_trunc) = Ai(to_trunc)./(min_signal_baseline_ratio + betai(to_trunc));

    % ========= Spatial-temporal tunings =======
    % [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
    %  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline

    tuning_theta = @(theta) coherence* Ai.* exp(ki.*(cos(theta_pref - theta)-1))* (speed_t)+ ...
        (1 - betai * speed_t * (coherence-0.05)) .* repmat(Bi,1,length(time));

    % [f']: Derivative of tuning (independent of ReLU)
    tuning_theta_der = @(theta) localSharpening * coherence* Ai.* exp(ki.*(cos(theta_pref - theta)-1)).* ...
        (ki .* sin(theta_pref - theta)) * speed_t;

   %%% CAREFUL: this part of the code does not use the coherence dependent tuning curve shown above
    % Override tuning baseline (Assuming we use ReLUs for each neuron to remove their minimal tuning curve) HH20180627
%     % %{
%     if ReLU == 1 % Remove to (1-beta)*B
%         tuning_theta = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
%             (1 - betai * speed_t ) .* repmat(Bi,1,length(time)) ...
%             - repmat((Ai .* exp(-2*ki) < betai.* Bi).* (Ai.*exp(-2*ki) * max(speed_t) + (1 - betai * max(speed_t)) .* Bi) ...  % If normal baseline Ai*exp(-2ki) < betai*Bi
%             + (Ai .* exp(-2*ki) >= betai.* Bi).* Bi, ...   % If the original baseline moves upwards
%             1,length(time));
%     elseif ReLU == 2 % Remove to B (cut off some of the null tuning)
%         tuning_theta = @(theta) max(eps, Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
%             (1 - betai * speed_t ) .* repmat(Bi,1,length(time)) ...
%             - repmat(Bi, 1, length(time)));
%     end
%     %}



    % ========== Adding Correlations =========
    resp = tuning_theta(theta_stim);
    resp_der = tuning_theta_der(theta_stim);

    C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)

    SIGMA_0_ilPPC = zeros(nu_neuron,nu_neuron);
    SIGMA_epsi_ilPPC = zeros(nu_neuron,nu_neuron);

    I_0_ts = nan(1,length(time));
    I_epsi_ts = nan(1,length(time));

    info_neuron_t = zeros(nu_neuron,length(time));


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
            %I_epsi_ts(tt) = f_der' * inv(SIGMA_epsi_ts) * f_der;   % Could be replaced by the fact that I_epsi = I_0 / (1 + eta * I_0)
            I_epsi_ts(tt) = I_0_ts(tt) / (1 + epsi_this * I_0_ts(tt));
        end

        % single neuron fisher information at time tt
        info_neuron_t(:,tt) =  f_der.^2./diag(SIGMA_epsi_ts);

        % ====== Adding momentary covariance matrix for ilPPC (assuming independent noise over time)
        SIGMA_0_ilPPC = SIGMA_0_ilPPC + SIGMA_0_ts;
        SIGMA_epsi_ilPPC = SIGMA_epsi_ilPPC + SIGMA_epsi_ts;
    end

    % ======= 1. Sum of all momentary info (standard of optimality) ========
    result.I_0_optimal = sum(I_0_ts) * dt;  % Area under momentary info (so that total info does not dependent on dt)
    result.I_epsi_optimal = sum(I_epsi_ts) * dt;

    % Total single neuron information
    infoperneuron = sum(info_neuron_t,2) * dt;

    %     % For checking momentary info temporally
    %         I_0_optimal(nn) = I_0_ts(round(time_max/2/dt));
    %         I_epsi_optimal(nn) = I_epsi_ts(round(time_max/2/dt));

    % ======= 2. Info of straight sum of spikes (ilPPC) ========
    f_der_ilPPC = sum(resp_der,2);

    result.I_0_ilPPC = f_der_ilPPC' * (SIGMA_0_ilPPC \ f_der_ilPPC) * dt;
    result.I_epsi_ilPPC = f_der_ilPPC' * (SIGMA_epsi_ilPPC \ f_der_ilPPC) * dt;

    % Threshold
    % ========= Prediction of psychothreshold ========
    result.I_epsi_Saturation =  sum(speed_t) ./ epsilon_0 * dt; % Total epsi info = sum(momentary epsi info)*dt = sum(1/momentary_epsi)*dt = sum(speed_t/epsis)*dt
    result.psycho_Saturation = rad2deg(sqrt(1./result.I_epsi_Saturation) * sqrt(2)); % Fine discrimination threshold (84% correct) = sqrt(2) * sqrt(1/I_fisher)
    result.psycho_ilPPC = rad2deg(sqrt(1./result.I_epsi_ilPPC) * sqrt(2));

    info_PPC(count) = result.I_epsi_ilPPC;
    info_saturation(count) = result.I_epsi_Saturation;
    info_optimal(count) = result.I_epsi_optimal;
    percent_info_loss(count) = (1-result.I_epsi_ilPPC/result.I_epsi_optimal)*100;
    psych_thres_optimal(count) = result.psycho_Saturation;
    psych_thres_PPC(count) = result.psycho_ilPPC;

    % Other outputs
    result.I_0_ts = I_0_ts;
    result.I_epsi_ts = I_epsi_ts;
    result.time = time;


end

%%
neuron_count = base_exp.^[2:count_max+1];

figure(1)
plot(resp');
%plot(mean(resp));

figure(2)
plot(neuron_count,info_PPC)
hold on
plot(neuron_count,info_optimal,'r')
hold off

figure(3)
plot(neuron_count,psych_thres_PPC)
hold on
plot(neuron_count,psych_thres_optimal,'r')
hold off

figure(4)
plot(neuron_count,percent_info_loss);

figure(5)
thresperneuron = rad2deg(sqrt(1./infoperneuron) * sqrt(2));

neuron_animal_ratio = log10(thresperneuron/psych_thres_optimal(count_max));

histogram(neuron_animal_ratio,0:0.1:max(neuron_animal_ratio));
title(sprintf('median = %g',median(neuron_animal_ratio)));
% hold on
% plot([1 1],[0 200]);
% hold off

fprintf('thres saturation: %f \n',psych_thres_optimal(count_max));
fprintf('Info saturation: %f \n',info_saturation(count_max));
fprintf('Info optimal: %f \n',info_optimal(count_max));
fprintf('Rho: %f \n',rho);
fprintf('Fano: %f \n',fano);
fprintf('Coherence: %f \n',coherence);
fprintf('Percent loss: %f \n',percent_info_loss(count_max)); 


if plotTuning
   %% Plot tuning  % HH
 
    % %{
    % --- Population activity
    % figure(1); clf; imagesc(tuning_theta(0)); colorbar
    
    example_time = round((0.15:0.2:time_max-0.2)/dt)+1;
    thetas_stim = -pi:pi/8:0; thetas_stim = [thetas_stim fliplr(-thetas_stim(2:end-1))];
    
    % --- Get all tuning curves. Because we have much more theta than t, I let thetas_stim be the vectorized term
    all_tunings = nan(nu_neuron,length(thetas_stim),length(example_time));
    
    for tt = 1:length(example_time)
        
        %     tuning_theta = @(theta) coherence* Ai.* exp(ki.*(cos(theta_pref - theta)-1))* (speed_t)+ ...
        %         (1 - betai * speed_t * (coherence-0.05)) .* repmat(Bi,1,length(time));
        
    
        all_tunings(:,:,tt) = (  coherence* Ai*ones(1,length(thetas_stim))).* exp( (ki*ones(1,length(thetas_stim))) ...
            .*(cos(theta_pref(:) * ones(1,length(thetas_stim)) - ones(length(theta_pref),1) * thetas_stim)-1)) ...
            * speed_t(example_time(tt)) + repmat((1 - betai * speed_t(example_time(tt))* (coherence-0.05)) .* Bi,1,length(thetas_stim));
        
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
    
    figure(20); clf;
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
    
    figure(21); clf;
    set(gcf,'uni','norm','pos',[0.092       0.229       0.765       0.526]);
    
    subplot(2,4,[1 2 5 6]);
    hold on;  set(gca,'colororder',color_order);
    errorbar(repmat([thetas_stim pi]'/pi*180,1,size(mean_tuning_align,2)),...
        [mean_tuning_align; mean_tuning_align(1,:)], [se_tuning_align; se_tuning_align(1,:)],'linew',2);
    set(gca,'xtick',-180:90:180); ylim([0 60])
    title(['n=' num2str(length(average_neuron))]);
    set(gca,'color',[0.8 0.8 0.8])
     
    subplot(2,4,3);
    hist(coherence * Ai,30); title('Effective Amp');  % HH
    
    subplot(2,4,4);
    half_widths = acos(1+log(0.5)./ki)*180/pi*2; % Change to be the same as Gu
    half_widths(imag(half_widths)~=0) = 360;
    hist(half_widths, 30); title('Tuning half width');
    
    subplot(2,4,7);
    hist(Bi,30); title('Baseline');
    
    subplot(2,4,8);
    hist(betai * (coherence - 0.05),30); title('Effective beta'); % HH
    xlim([0 1])
    
    figure(2045); clf
    plot(theta_pref(1:4:end)/pi*180,half_widths(1:4:end),'ok')
    set(gca,'xtick',-180:90:180,'ytick',45:45:270)
    ylim([45 280])

        %}
end
