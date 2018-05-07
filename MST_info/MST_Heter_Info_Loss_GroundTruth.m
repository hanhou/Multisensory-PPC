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
%{
    epsis = 0.0027;
    mean_baselines = 20;    
    std_baselines = 5;
    nu_neurons = round(10.^linspace(1,log10(1000),10)) + 1;
%}


% ==== Scan epsis and baseline ===
% %{
    epsis = [0 10.^(linspace(-4,0,10))];
    std_baselines = 5;
    mean_baselines = linspace(0,30,10);
    nu_neurons = 500;
%}


% ==== Scan baseline and var_baseline ===
%{
    epsis = 0.0027;
    std_baselines = linspace(eps,30,10);
    mean_baselines = linspace(0,30,10);
    nu_neurons = 500;
%}


if ~ION_cluster ;progressbar(0); end
count = 0; 

tic
for ee = 1:length(epsis) % Scan epsi
    for mbmb = 1:length(mean_baselines) % Scan baseline
        for sbsb = 1:length(std_baselines) % Scan std_baseline
            for nn = 1:length(nu_neurons)  % Scan nu_neurons
                count = count + 1;
                
                % ======== Neurons and Times ========
                nu_neuron = nu_neurons(nn);
                
                % ---- Uniformly distrubuted preferred heading ----
                theta_pref = linspace(-pi,pi,nu_neuron)';
                
                % ---- Bias to +/- 90 degrees (Gu 2006) ---- [No changing of % recovery]
                % theta_pref = [randn(1,round(nu_neuron/2))*(1.2*pi/4)-(pi/2) randn(1,nu_neuron-round(nu_neuron/2))*(1.2*pi/4)+(pi/2)]';
                % theta_pref = mod(theta_pref,2*pi)-pi;
                
                time_max = 2000; % ms
                dt = 10; % ms
                time = 0:dt:time_max;
                
                % Use the same setting as experiments
                NumOfSigma = 3.5;
                sigma_t = time_max/2/NumOfSigma;
                speed_t = exp(-(time-time_max/2).^2/(2*sigma_t^2));
                speed_t = speed_t-min(speed_t);
                speed_t = speed_t/max(speed_t);
                
                position_t = cumsum(speed_t);
                position_t = position_t-min(position_t);
                position_t = position_t/max(position_t);
                
                % ======== Parameters =======
                % Using gamma distribution to introduce heterogeneity
                gamma_para = [ % Mean   Std  (Desired)
                    45,    20;  % amp
                    1.4,   0.3;  % k (controls sigma of spatial tuning)
                    0.7,   0.5; % beta
                    mean_baselines(mbmb),  std_baselines(sbsb);  % original 20, 5
                    ];
                
                %     gamma_para(:,2) = eps; % Turn off heterogeneity for test
                
                gamma_A = (gamma_para(:,1)./gamma_para(:,2)).^2;
                gamma_B = (gamma_para(:,1))./gamma_A;    gamma_B(isnan(gamma_B)) = 0;
                tuning_paras = nan(nu_neuron,size(gamma_para,1));
                for pp = 1:size(gamma_para,1)
                    tuning_paras(:,pp) = gamrnd(gamma_A(pp),gamma_B(pp),nu_neuron,1);
                end
                
                tuning_paras(tuning_paras(:,3)>1,3) = 1; % beta should be smaller than 1
                
                % ========= Spatial-temporal tunings =======
                % [f]: Population activity as a function of stimulus theta. Baseline term accounts for what we found in Gu MST data
                %  f(theta,t) = amp * v(t) * exp {k*(cos(theta-theta_pref)-1)} + (1 - beta * v(t)/max(v)) * baseline
                tuning_theta = @(theta) tuning_paras(:,1).* exp(tuning_paras(:,2).*(cos(theta_pref - theta)-1)) * (speed_t) + ...
                    (1 - tuning_paras(:,3) * speed_t ) .* repmat(tuning_paras(:,4),1,length(time));
                % [f']: Derivative of tuning
                tuning_theta_der = @(theta) tuning_paras(:,1).* exp(tuning_paras(:,2).*(cos(theta_pref - theta)-1)).* ...
                    (tuning_paras(:,2) .* sin(theta_pref - theta)) * speed_t;
                % -- Plot f'--
                % fp = tuning_theta_prime(0);
                % figure(); plot(theta_pref/pi*180,fp(:,1000))
                
                % ========== Adding Correlations =========
                theta_stim = 0; % Calculate Info around 0 degree
                rho = 0.1;  k = 2;  fano = 1;  % Decay correlation
                % epsi = 0.0027; % How much f'f'T noise
                
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
                    epsi_this = epsis(ee) / (speed_t(tt) + 1e-6);
                    
                    % ======= Add correlation ======
                    % Exponentially decay correlation (non-differential in the heterogeneous case)
                    SIGMA_0_ts = fano * C .* sqrt(f * f');
                    
                    % Total covariance matrix = SIGMA_0 + Differential Correlation
                    SIGMA_epsi_ts = SIGMA_0_ts + epsi_this * (f_der * f_der');
                    
                    % ====== Compute Linear Fisher info directly by groud truth I = f'T SIGMA^-1 f' =====
                    if sum(f) == 0
                        I_0_ts(tt) = 0;
                        I_epsi_ts(tt) = 0;
                    else
                        I_0_ts(tt) = f_der' * inv(SIGMA_0_ts) * f_der;
                        % I_epsi_ts(tt) = f_der' * inv(SIGMA_epsi_ts) * f_der;   % Could be replaced by the fact that I_epsi = I_0 / (1 + eta * I_0)
                        I_epsi_ts(tt) = I_0_ts(tt) / (1 + epsi_this * I_0_ts(tt));
                    end
                    
                    % ====== Adding momentary covariance matrix for ilPPC (assuming independent noise over time)
                    SIGMA_0_ilPPC = SIGMA_0_ilPPC + SIGMA_0_ts;
                    SIGMA_epsi_ilPPC = SIGMA_epsi_ilPPC + SIGMA_epsi_ts;
                end
                
                % ======= 1. Sum of all momentary info (standard of optimality) ========
                I_0_optimal(nn) = sum(I_0_ts) * dt/1000;  % Area under momentary info (so that total info does not dependent on dt)
                I_epsi_optimal(nn) = sum(I_epsi_ts) * dt/1000;
                
                %     % For checking momentary info temporally
                %         I_0_optimal(nn) = I_0_ts(round(time_max/2/dt));
                %         I_epsi_optimal(nn) = I_epsi_ts(round(time_max/2/dt));
                
                % ======= 2. Info of straight sum of spikes (ilPPC) ========
                f_der_ilPPC = sum(resp_der,2);
                
                I_0_ilPPC(nn) = f_der_ilPPC' * inv(SIGMA_0_ilPPC) * f_der_ilPPC * dt/1000;
                I_epsi_ilPPC(nn) = f_der_ilPPC' * inv(SIGMA_epsi_ilPPC) * f_der_ilPPC * dt/1000;
                
                if ~ION_cluster 
                    progressbar(count/(length(epsis)*length(mean_baselines)*length(nu_neurons)*length(std_baselines)));
                else
                    fprintf('%3.1f\n',(count/(length(epsis)*length(mean_baselines)*length(nu_neurons)*length(std_baselines))*100));
                end
            end
            
            optimality_matrix(ee,mbmb,sbsb) = mean(I_epsi_ilPPC./I_epsi_optimal)*100;
            
            %             if isnan(optimality_matrix(ee,mbmb,sbsb))
            %                 keyboard;
            %             end
        end
    end
end
toc

% ====== Plot momentary information ========
%%{
figure(10); plot(time,I_0_ts,'r','linew',2)
hold on; plot(time,I_epsi_ts,'b','linew',2);
hold on; plot(time,speed_t*max(I_0_ts),'k--');
%}

% ====== Compare ilPPC with optimal =======
figure(11); 
% Optimal
plot(nu_neurons, I_0_optimal,'ro-','linew',3,'markerfacecol','r');
hold on; plot(nu_neurons, I_epsi_optimal,'bo-','linew',3,'markerfacecol','b'); 

% ilPPC
plot(nu_neurons, I_0_ilPPC,'ro--','linew',3,'markerfacecol','r');
plot(nu_neurons, I_epsi_ilPPC,'bo--','linew',3,'markerfacecol','b');
thres_discrim_epsi = sqrt(2/I_epsi_ilPPC(end))/pi*180;
set(text(max(ylim)*0.9,I_epsi_ilPPC(end)*0.8,sprintf('\\sigma = %g',thres_discrim_epsi)),'color','b');

xlabel('Number of neurons');
% plot(xlim,1/epsi*[1 1],'k','linew',2);  
axis tight

% ====== Optimality of ilPPC =========
figure(12); % Recovery of info
plot(nu_neurons,I_0_ilPPC./I_0_optimal*100,'ro-','linew',3,'markerfacecol','r'); hold on
plot(nu_neurons,I_epsi_ilPPC./I_epsi_optimal*100,'bo-','linew',3,'markerfacecol','b');
axis tight; plot(xlim,[100 100],'k--','linew',2);
ylim([50 105]);
xlabel('Number of neurons');
ylabel('% of info recovery by ilPPC');

if ION_cluster
    file_name = sprintf('./0_Optimality 1D');
    export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
    saveas(gcf,[file_name '.fig'],'fig');
end



% ========= Plot optimality matrix =========
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
    colormap(hot);
    xlabel('log_{10}(\epsilon_0)');
    ylabel('Mean(Baseline)');
elseif length(mean_baselines)>1 && length(std_baselines)>1
    figure(13); clf;
    imagesc(std_baselines,mean_baselines,optimality_matrix); axis xy;
    hold on;
    contour(std_baselines,mean_baselines, optimality_matrix,'color','k','linew',1.5,'ShowText','on');
    
%     set(gca,'ytick',0:5:max(mean_baselines))
    colormap(hot);
    xlabel('Std(Baseline)');
    ylabel('Mean(Baseline)');

end


if ION_cluster
    file_name = sprintf('./1_optimality_matrix');
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

 

