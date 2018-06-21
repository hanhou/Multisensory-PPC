% Using GROUND TRUTH to calculate linear Fisher information in heterogenerous
% MST population with differential correlation.
%
% First by Han Hou @ 20180426,  houhan@gmail.com

clear;
rng('shuffle');
addpath(genpath('/home/hh/Multisensory-PPC/util'));

% ====== For cluster running 
hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name


if ~isempty(strfind(hostname,'node')) || ~isempty(strfind(hostname,'clc')) % contains(hostname,{'node','clc'})  % ION cluster;
%     if isempty(gcp('nocreate'))
%         parpool;
% %         parpool(hostname(1:strfind(hostname,'.')-1),20);
%     end
    ION_cluster = 1;
else % My own computer
%     if strcmp(version('-release'),'2013a')
%         if matlabpool('size') == 0
%             matlabpool;
%         end
%     else
%         if isempty(gcp('nocreate'))
%             parpool;
%         end
%     end
    ION_cluster = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%% Defult Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Tunings 
nu_neurons = 1000;

mean_baselines = 60;
std_baselines = 40;

% Using gamma distribution to introduce heterogeneity
gamma_para = [ % Mean   Std  (Desired)
                100,    70;  % amp
                1.7,   0.9;  % Controls sigma of spatial tuning: k = log(0.5)/(cos(degtorad(desired_half_width))-1)
                0.3,   0.2; % beta
                mean_baselines,   std_baselines;  % original 20, 5
             ];
         
% 2. Correlations 
rho = 0.2;  k = 2;  fano = 1;  % Decay correlation
epsis = 0.00139; % f'f'T  Corresponds to threshold ~= 4 degree. Original 0.0027

% 3. Miscs
NumOfSigma = 3.5;
delay = 0.140;
theta_stim = 0; % Calculate Info around 0 degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ==== Scan nu_neurons ===
%{
    nu_neurons = round(10.^linspace(1,log10(1000),10)) + 1; 
%}

% ==== Scan epsis and baseline ===
% %{
    epsis = [0 10.^(linspace(-4,-2,10))];
    mean_baselines = mean_baselines; % Only care about Fig.4
    % mean_baselines = linspace(0,2*mean_baselines,5);  % Should be odd number because I'm gonna plot the mean(mean_baselines)!
%} 

% ==== Scan baseline and var_baseline ===
%{
    std_baselines = linspace(eps,40,15);
    mean_baselines = linspace(0,40,15);
%}


if ION_cluster
    if length(nu_neurons)>1 % Scan 1-D nu_neuron
        nu_runs = 10;
    else % Scan 2-D
        nu_runs = 2;
    end
else
    nu_runs = 1;
end

% --- Motion profile ---
time_max = 2; % s
if ION_cluster
    dt = 50e-3; % s. This could be large because we're not simulating Poisson spikes.
else
    dt = 100e-3;
end
time = 0:dt:time_max;

% Use the same setting as experiments
sigma_t = time_max/2/NumOfSigma;
speed_t = exp(-(time-time_max/2-delay).^2/(2*sigma_t^2));
speed_t = speed_t-min(speed_t);
speed_t = speed_t/max(speed_t);

position_t = cumsum(speed_t);
position_t = position_t-min(position_t);
position_t = position_t/max(position_t);


count = 0; 
optimality_matrix = [];

tic
for ee = 1:length(epsis) % Scan epsi
    epsilon_0 = epsis(ee);
    for mbmb = 1:length(mean_baselines) % Scan baseline
        for sbsb = 1:length(std_baselines) % Scan std_baseline
            
            % ======== Parameters =======
            
            % Override baseline
            gamma_para(end,:) = [ mean_baselines(mbmb),  std_baselines(sbsb) ];  % original 20, 5
            
            % Override heterogeneity
            %     gamma_para(:,2) = eps; % Turn off heterogeneity for test

            gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
            gamma_B = (gamma_para(:,1))./gamma_A;    gamma_B(isnan(gamma_B)) = 0;

            I_0_optimal = nan(nu_runs, length(nu_neurons));
            I_0_ilPPC = I_0_optimal;
            I_epsi_optimal = I_0_optimal;
            I_epsi_ilPPC = I_0_optimal;
            
            for nn = 1:length(nu_neurons)  % Scan nu_neurons
                count = count + 1;
                
                % ======== Neurons and Times ========
                nu_neuron = nu_neurons(nn);
                
                % ---- Uniformly distrubuted preferred heading ----
                %                 theta_pref = linspace(-pi,pi,nu_neuron)';
                
                % ---- Bias to +/- 90 degrees (Gu 2006) ---- [Little changing of % recovery]
                theta_pref = [randn(1,round(nu_neuron/2))*(1.2*pi/4)-(pi/2) randn(1,nu_neuron-round(nu_neuron/2))*(1.2*pi/4)+(pi/2)]';
                theta_pref = mod(theta_pref,2*pi)-pi;
                
                fprintf('\nParameter scan %.2g%%: ',(count/(length(epsis)*length(mean_baselines)*length(nu_neurons)*length(std_baselines))*100));
                
                
                %------------  Runs ------------
                for rruns = 1:nu_runs
                    
                    fprintf('>');
                    
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
%                     min_signal_baseline_ratio = 0.5;
%                     to_trunc = find(Bi > Ai./(min_signal_baseline_ratio + betai));
%                     Bi (to_trunc) = Ai(to_trunc)./(min_signal_baseline_ratio + betai(to_trunc));
                    
                    % ========= Spatial-temporal tunings =======
                    % [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
                    %  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline
                    tuning_theta = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)) * (speed_t) + ...
                        (1 - betai * speed_t ) .* repmat(Bi,1,length(time));
                    % [f']: Derivative of tuning
                    tuning_theta_der = @(theta) Ai.* exp(ki.*(cos(theta_pref - theta)-1)).* ...
                        (ki .* sin(theta_pref - theta)) * speed_t;
                    % -- Plot f'--
                    % fp = tuning_theta_prime(0);
                    % figure(); plot(theta_pref/pi*180,fp(:,1000))
                    
                    % ========== Adding Correlations =========                   
                    resp = tuning_theta(theta_stim);
                    resp_der = tuning_theta_der(theta_stim);
                    
                    C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)
                    
                    SIGMA_0_ilPPC = zeros(nu_neuron,nu_neuron);
                    SIGMA_epsi_ilPPC = zeros(nu_neuron,nu_neuron);
                    
                    I_0_ts = nan(1,length(time));
                    I_epsi_ts = nan(1,length(time));
                    
                    parfor tt = 1:length(time)    % Across all trial
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
                        SIGMA_0_ilPPC = SIGMA_0_ilPPC + SIGMA_0_ts;
                        SIGMA_epsi_ilPPC = SIGMA_epsi_ilPPC + SIGMA_epsi_ts;
                    end
                    
                    % ======= 1. Sum of all momentary info (standard of optimality) ========
                    I_0_optimal(rruns, nn) = sum(I_0_ts) * dt;  % Area under momentary info (so that total info does not dependent on dt)
                    I_epsi_optimal(rruns, nn) = sum(I_epsi_ts) * dt;
                    
                    %     % For checking momentary info temporally
                    %         I_0_optimal(nn) = I_0_ts(round(time_max/2/dt));
                    %         I_epsi_optimal(nn) = I_epsi_ts(round(time_max/2/dt));
                    
                    % ======= 2. Info of straight sum of spikes (ilPPC) ========
                    f_der_ilPPC = sum(resp_der,2);
                    
                    I_0_ilPPC(rruns, nn) = f_der_ilPPC' * (SIGMA_0_ilPPC \ f_der_ilPPC) * dt;
                    I_epsi_ilPPC(rruns, nn) = f_der_ilPPC' * (SIGMA_epsi_ilPPC \ f_der_ilPPC) * dt;
                    
                end
                %------------ End Runs ------------
                
            end
            
            if length(nu_neurons) == 1  % Optimality_matrix will never nested with nu_neurons
                optimality_matrix(ee,mbmb,sbsb) = mean(I_epsi_ilPPC./I_epsi_optimal,1)*100;
                optimality_std_matrix(ee,mbmb,sbsb) = std(I_epsi_ilPPC./I_epsi_optimal,[],1)*100;
                threshold_matrix(ee,mbmb,sbsb) =  mean(sqrt(2./I_epsi_ilPPC(:,end))/pi*180,1);
                threshold_std_matrix(ee,mbmb,sbsb) = std(sqrt(2./I_epsi_ilPPC(:,end))/pi*180,[],1);
            end
            
            %             if isnan(optimality_matrix(ee,mbmb,sbsb))
            %                 keyboard;
            %             end
        end
    end
end

fprintf('\n')
toc

% ========= Prediction of psychothreshold ========
info_pred =  sum(speed_t) ./ epsis * dt; % Total epsi info = sum(momentary epsi info)*dt = sum(1/momentary_epsi)*dt = sum(speed_t/epsis)*dt
psycho_pred = rad2deg(sqrt(1./info_pred) * sqrt(2)); % Fine discrimination threshold (84% correct) = sqrt(2) * sqrt(1/I_fisher)



%% Plotting

% ====== Plot momentary information ========

if length(nu_neurons) > 1
    % %{
    figure(10); plot(time,I_0_ts,'r','linew',2)
    hold on; plot(time,I_epsi_ts,'b','linew',2);
    hold on; plot(time,speed_t*max(I_0_ts),'k--');
    %}
    
    % ====== Compare ilPPC with optimal =======
    figure(11);  clf;
    % Optimal
    % plot(nu_neurons, I_0_optimal,'ro-','linew',3,'markerfacecol','r');
    % hold on; plot(nu_neurons, I_epsi_optimal,'bo-','linew',3,'markerfacecol','b');
    errorbar(nu_neurons, mean(I_0_optimal,1),std(I_0_optimal,[],1),'ro-','linew',3,'markerfacecol','r');
    hold on; errorbar(nu_neurons, mean(I_epsi_optimal,1),std(I_epsi_optimal,[],1),'bo-','linew',3,'markerfacecol','b');
    plot(xlim,[info_pred info_pred],'k--');
    
    % ilPPC
    % plot(nu_neurons, I_0_ilPPC,'ro--','linew',3,'markerfacecol','r');
    % plot(nu_neurons, I_epsi_ilPPC,'bo--','linew',3,'markerfacecol','b');
    
    errorbar(nu_neurons, mean(I_0_ilPPC,1), std(I_0_ilPPC,[],1),'ro--','linew',3,'markerfacecol','r');
    errorbar(nu_neurons, mean(I_epsi_ilPPC,1), std(I_epsi_ilPPC,[],1), 'bo--','linew',3,'markerfacecol','b');
    
    thres_discrim_epsi = sqrt(2/mean(I_epsi_ilPPC(:,end),1))/pi*180;
    set(text(nu_neuron,I_epsi_ilPPC(end),sprintf('\\sigma = %g',thres_discrim_epsi)),'color','b');
    fprintf('Saturated sigma ~= %g\n',thres_discrim_epsi);
    
    if ION_cluster
        file_name = sprintf('./0_Raw Info_GroundTruth');
        % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end
    
    
    xlabel('Number of neurons');
    % plot(xlim,1/epsi*[1 1],'k','linew',2);
    axis tight
    
    % ====== Optimality of ilPPC =========
    figure(12);  clf;
    
    % Recovery of info
    % plot(nu_neurons,I_0_ilPPC./I_0_optimal*100,'ro-','linew',3,'markerfacecol','r'); hold on
    % plot(nu_neurons,I_epsi_ilPPC./I_epsi_optimal*100,'bo-','linew',3,'markerfacecol','b');
    
    errorbar(nu_neurons,mean(I_0_ilPPC./I_0_optimal,1)*100, std(I_0_ilPPC./I_0_optimal,[],1)*100, 'ro-','linew',3,'markerfacecol','r'); hold on
    errorbar(nu_neurons,mean(I_epsi_ilPPC./I_epsi_optimal,1)*100, std(I_epsi_ilPPC./I_epsi_optimal,[],1)*100,'bo-','linew',3,'markerfacecol','b');
    
    axis tight; plot(xlim,[100 100],'k--','linew',2);
    ylim([50 105]);
    xlabel('Number of neurons');
    ylabel('% of info recovery by ilPPC');
    
    fprintf('Saturated optimality = %g%%\n', mean(I_epsi_ilPPC(:,end)./I_epsi_optimal(:,end),1))
    
    if ION_cluster
        file_name = sprintf('./1_Optimality 1D_GroundTruth');
        saveas(gcf,[file_name '.fig'],'fig');
    end
    
end


if length(epsis)>1
    
    % ========= Plot optimality matrix =========
    optimality_matrix = squeeze(optimality_matrix);
    threshold_matrix = squeeze(threshold_matrix);
    
    figure(13); clf; set(gcf,'uni','norm','pos',[0.004        0.29       0.987       0.438]);
    epsi_log_x = [3*log10(epsis(2))-2*log10(epsis(3)) log10(epsis(2:end))]; % Add one entry for epsi = 0
    epsi_log_x_tick = [min(epsi_log_x) linspace(round(epsi_log_x(2)),round(max(epsi_log_x)),5)];
    epsi_log_x_label = [-inf epsi_log_x_tick(2:end)];
    
    if length(mean_baselines) > 1
        
        subplot(1,3,1) % Optimality
        imagesc(epsi_log_x,mean_baselines,optimality_matrix'); axis xy; hold on;
        contour(epsi_log_x,mean_baselines, optimality_matrix','color','k','linew',1.5,'ShowText','on');
        
        set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)
        set(gca,'ytick',0:5:max(mean_baselines))
        colormap(parula); caxis([50 100])
        xlabel('log_{10}(\epsilon_0)');
        ylabel('Mean(Baseline)');
        
        
        subplot(1,3,2)  % Threshold
        imagesc(epsi_log_x,mean_baselines, threshold_matrix'); axis xy; hold on;
        contour(epsi_log_x,mean_baselines, threshold_matrix','color','k','linew',1.5,'ShowText','on');
        
        set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)
        set(gca,'ytick',0:5:max(mean_baselines))
        colormap(hot);
        xlabel('log_{10}(\epsilon_0)');
        ylabel('Mean(Baseline)');
    end
    
    subplot(1,3,3)  % Threshold: prediction and simulation 
    to_plot_baseline = mean(mean_baselines); % 20
    [~, ind_baseline] = min(abs(to_plot_baseline - mean_baselines));
    
    optimality_to_plot_baseline = optimality_matrix(:,ind_baseline);
    
    [ax, h1, h2] = plotyy( epsi_log_x, optimality_to_plot_baseline, epsi_log_x, psycho_pred ); 
    set(h2,'color','b', 'linestyle','--', 'marker', 's', 'linew',2)
    
    hold(ax(1),'on') ; delete(h1); 
    errorbar(ax(1), epsi_log_x, optimality_matrix(:,ind_baseline)', optimality_std_matrix(:,ind_baseline)','ro-','linew',2);
    set(ax(1),'ytick',50:5:100); ylim(ax(1), [50 102]);
    
    hold(ax(2),'on');
    errorbar(ax(2), epsi_log_x, threshold_matrix(:,ind_baseline), threshold_std_matrix(:,ind_baseline),'ko-','linew',2)
    errorbar(ax(2), epsi_log_x, threshold_matrix(:,1), threshold_std_matrix(:,1),'bo-','linew',2)
    legend({'optimality', 'thres\_inf&optimal', sprintf('thres (B=%g)',to_plot_baseline),'thres (B=0)', })
    plot(ax(1),xlim,[100 100],'r--','linew',2);
    plot(xlim,[4 4],'--k'); 
   
    epsi_log_x_label = [-inf epsi_log_x_tick(2:end)];
    set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)

    axis(ax(2),'tight'); set(ax(2),'ytick',0:10)
    linkaxes(ax,'x');
    set(ax(1),'ycolor','r')
    set(ax(2),'ycolor','b')
    
    if ION_cluster
        saveas(13,'./4_Optimality_matrix_GroundTruth');
    end
    
    
elseif length(mean_baselines)>1 && length(std_baselines)>1
    
    % ========= Plot optimality matrix =========
    optimality_matrix = squeeze(optimality_matrix);
    threshold_matrix = squeeze(threshold_matrix);
    
    figure(13); clf;  set(gcf,'uni','norm','pos',[0.118       0.333       0.726       0.453]);
    
    subplot(1,2,1)
    imagesc(std_baselines,mean_baselines,optimality_matrix); axis xy; hold on;
    contour(std_baselines,mean_baselines, optimality_matrix,'color','k','linew',1.5,'ShowText','on');
    
    %     set(gca,'ytick',0:5:max(mean_baselines))
    colormap(parula); % caxis([50 100])
    xlabel('Std(Baseline)');
    ylabel('Mean(Baseline)');
    
    subplot(1,2,2)
    imagesc(std_baselines,mean_baselines, threshold_matrix); axis xy; hold on;
    contour(std_baselines,mean_baselines, threshold_matrix,'color','k','linew',1.5,'ShowText','on');
    
    %     set(gca,'ytick',0:5:max(mean_baselines))
    % caxis([50 100])
    xlabel('Std(Baseline)');
    ylabel('Mean(Baseline)');
    
    
    if ION_cluster
        saveas(13,'./4_Optimality_matrix_GroundTruth');
    end
    
end



%%
% ========== Tuning properties =========
if length(nu_neurons) > 1
    
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
        title(sprintf('%0.1f/',tuning_paras(example_neuron(i),:)));
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
    half_widths = acos(1+log(0.5)./ki)*180/pi;
    half_widths(imag(half_widths)~=0) = 180;
    hist(half_widths, 30); title('Tuning half width');
    
    subplot(2,4,7);
    hist(Bi,30); title('Baseline');
    
    subplot(2,4,8);
    hist(betai,30); title('beta');
    
    if ION_cluster
        file_name = sprintf('./3_average tuning');
        % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
        saveas(gcf,[file_name '.fig'],'fig');
    end
    %}
    
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

 

