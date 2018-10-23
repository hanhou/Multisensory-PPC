% Using GROUND TRUTH to calculate linear Fisher information in heterogenerous
% MST population with differential correlation.
%
% First by Han Hou @ 20180426,  houhan@gmail.com

clear;
rng('default');
rng('shuffle');
addpath(genpath('/home/hh/Multisensory-PPC/util'));

% ====== For cluster running
hostname = char( getHostName( java.net.InetAddress.getLocalHost)); % Get host name

if [strfind(hostname,'node') strfind(hostname,'clc')] % contains(hostname,{'node','clc'})  % ION cluster;
    ION_cluster = 1;
else % My own computer
    ION_cluster = 0;
end


if ION_cluster
    nu_runs = 10;
else
    nu_runs = 5;
end

nu_neuron = 10000;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% doScan = 'Scan Number of Neurons';
%  doScan = 'Scan Number of Neurons and epsilons or rho';
% doScan = 'Scan epsis and baseline'; % Threshold as a function of epsis
doScan = 'Scan general linear 2-D';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch doScan
    case 'Scan Number of Neurons'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nu_neurons = round(10.^linspace(1,log10(nu_neuron),10)) + 1;
        ReLU = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Linearize parameters
        paras = [repmat({'nu_neuron'},length(nu_neurons),1),mat2cell(nu_neurons(:),ones(length(nu_neurons),1))];
        paras = repmat(paras, nu_runs,1);
        
        % Do Parallel
        infos = nan(size(paras,1),7);
        
        parfor_progress(size(paras,1));  tic
        parfor pp = 1:size(paras,1)
            
            if pp == size(paras,1)
                ifplot = 1;
            else
                ifplot = 0;
            end
            
            result = InfoLossGroundTruth([paras(pp,:), 'ReLU', ReLU, 'plotTuning', ifplot]);
            
            infos(pp,:) = cellfun(@(x)result.(x),{'I_0_optimal','I_0_ilPPC','I_epsi_optimal','I_epsi_ilPPC',...
                                                  'I_epsi_Saturation', 'psycho_Saturation', 'psycho_ilPPC'});
            parfor_progress();
        end
        parfor_progress(0);   toc
        
        
        I_0_optimal = reshape(infos(:,1),[],nu_runs)';
        I_0_ilPPC = reshape(infos(:,2),[],nu_runs)';
        I_epsi_optimal = reshape(infos(:,3),[],nu_runs)';
        I_epsi_ilPPC = reshape(infos(:,4),[],nu_runs)';
        info_pred = infos(1,5);
        
        % =========== Plotting =============
        
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
        set(text(nu_neurons(end),I_epsi_ilPPC(end),sprintf('\\sigma = %g',thres_discrim_epsi)),'color','b');
        fprintf('Saturated sigma ~= %g\n',thres_discrim_epsi);
        
        xlabel('Number of neurons');
        % plot(xlim,1/epsi*[1 1],'k','linew',2);
        axis tight

        if ION_cluster
            file_name = sprintf('./0_Raw Info_GroundTruth');
            % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
            saveas(gcf,[file_name '.fig'],'fig');
        end
        
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
        
    case 'Scan Number of Neurons and epsilons or rho'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         epsis_or_rho = [0 0.2 0.4 0.6 0.8];
%         epsis_or_rho_text = 'rho';
%         epsis_or_rho = [0.1 0.5 1 2 4];
%         epsis_or_rho_text = 'k';  fano
        epsis_or_rho = [eps 0.01 0.1 0.5 1 2];
        epsis_or_rho_text = 'fano'; 

%         epsis_or_rho = [0 10.^(linspace(-4,-1,6))];
%         epsis_or_rho_text = 'epsilon_0';
        nu_neurons = round(10.^linspace(1,log10(10000),7)) + 1;
        ReLU = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
        % Linearize parameters
        [xx, yy] = meshgrid(epsis_or_rho, nu_neurons);
        paras = [repmat({epsis_or_rho_text},numel(xx),1), mat2cell(xx(:),ones(numel(xx),1)),...
                 repmat({'nu_neuron'},numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];
        paras = repmat(paras, nu_runs,1);

        % Do Parallel
        infos = nan(size(paras,1),7);
        
        parfor_progress(size(paras,1));  tic
        for pp = 1:size(paras,1)
            result = InfoLossGroundTruth([paras(pp,:), 'ReLU', ReLU]);
            
            infos(pp,:) = cellfun(@(x)result.(x),{'I_0_optimal','I_0_ilPPC','I_epsi_optimal','I_epsi_ilPPC',...
                                                  'I_epsi_Saturation', 'psycho_Saturation', 'psycho_ilPPC'});
            parfor_progress();
        end
        parfor_progress(0);   toc
        
        % Plot example tuning using the last condition
        InfoLossGroundTruth([paras(end,:),'ReLU',ReLU, 'plotTuning',1]);

        
        I_0_optimal = reshape(infos(:,1),size(xx,1),size(xx,2),nu_runs);
        I_0_ilPPC = reshape(infos(:,2),size(xx,1),size(xx,2),nu_runs);
        I_epsi_optimal = reshape(infos(:,3),size(xx,1),size(xx,2),nu_runs);
        I_epsi_ilPPC = reshape(infos(:,4),size(xx,1),size(xx,2),nu_runs);
        I_epsi_optimal_inf = reshape(infos(:,5),size(xx,1),size(xx,2),nu_runs);
        
        %% =========== Plotting =============
        % ====== Compare ilPPC with optimal =======
        figure(11);  clf; hold on;
        set(gca,'yscale','log','xscale','log')

        % Optimal
%         if epsis_or_rho_text == 'rho'
            cols = jet;
            cols = flipud(cols(round(linspace(1,size(cols,1),length(epsis_or_rho))),:));
%         else
%             cols = lines;
%         end
        
        for ee = 1:length(epsis_or_rho)
            errorbar(nu_neurons,mean(I_epsi_optimal(:,ee,:),3)', ...
                                 std(I_epsi_optimal(:,ee,:),[],3)','o-','linew',3,'color',cols(ee,:));
                             
            I_inf_this = I_epsi_optimal_inf(1,ee,1);
            if I_inf_this ~= inf
                plot(xlim, ones(1,2) * I_inf_this, '--','color',cols(ee,:),'linew',2);
                thres_discrim_epsi_inf = sqrt(2/I_inf_this)/pi*180;
                set(text(nu_neurons(end), mean(I_epsi_optimal(end,ee,:)),...
                    sprintf('%s = %.2e\n\\sigma = %.2g',epsis_or_rho_text,epsis_or_rho(ee),thres_discrim_epsi_inf)),'color',cols(ee,:));
            end
        end
        
        % ilPPC
        % plot(nu_neurons, I_0_ilPPC,'ro--','linew',3,'markerfacecol','r');
        % plot(nu_neurons, I_epsi_ilPPC,'bo--','linew',3,'markerfacecol','b');
        
%         thres_discrim_epsi = sqrt(2/mean(I_epsi_ilPPC(:,end),1))/pi*180;
%         set(text(nu_neurons(end),I_epsi_ilPPC(end),sprintf('\\sigma = %g',thres_discrim_epsi)),'color','b');
%         fprintf('Saturated sigma ~= %g\n',thres_discrim_epsi);
        SetFigure(20);
        xlabel('Number of neurons');
        % plot(xlim,1/epsi*[1 1],'k','linew',2);
        axis tight

        if ION_cluster
            file_name = sprintf('./0_Raw Info_GroundTruth');
            % export_fig('-painters','-nocrop','-m1.5',[file_name '.png']);
            saveas(gcf,[file_name '.fig'],'fig');
        end
        
        % ====== Optimality of ilPPC =========
        figure(12);  clf;
        hold on;
        % Recovery of info
        % plot(nu_neurons,I_0_ilPPC./I_0_optimal*100,'ro-','linew',3,'markerfacecol','r'); hold on
        % plot(nu_neurons,I_epsi_ilPPC./I_epsi_optimal*100,'bo-','linew',3,'markerfacecol','b');
        
        % errorbar(nu_neurons,mean(I_0_ilPPC./I_0_optimal,1)*100, std(I_0_ilPPC./I_0_optimal,[],1)*100, 'ro-','linew',3,'markerfacecol','r'); hold on
        for ee = 1:length(epsis_or_rho)
            errorbar(nu_neurons,mean(I_epsi_ilPPC(:,ee,:)./I_epsi_optimal(:,ee,:),3)'*100, ...
                                 std(I_epsi_ilPPC(:,ee,:)./I_epsi_optimal(:,ee,:),[],3)'*100,'o-','linew',3,'color',cols(ee,:));
        end
        
        legend(num2str(epsis_or_rho'))
        set(gca,'xscale','log')
        axis tight; plot(xlim,[100 100],'k--','linew',2);
        ylim([50 105]);
        xlabel('Number of neurons');
        ylabel('% of info recovery by ilPPC');
        SetFigure(20);
        
        fprintf('Saturated optimality = %g%%\n', mean(I_epsi_ilPPC(:,end)./I_epsi_optimal(:,end),1))
        
        if ION_cluster
            file_name = sprintf('./1_Optimality 1D_GroundTruth');
            saveas(gcf,[file_name '.fig'],'fig');
        end

    case 'Scan epsis and baseline'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        epsis = [0 10.^(linspace(-4,-2,10))];
        mean_baselines = 20; % linspace(eps,40,3);  % Should be odd number because I'm gonna plot the mean(mean_baselines)!
        keepStdMeanRatio = 1;
        ReLU = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Linearize parameters
        [xx, yy] = meshgrid(epsis, mean_baselines);
        paras = [repmat({'epsilon_0'},numel(xx),1), mat2cell(xx(:),ones(numel(xx),1)),...
                 repmat({'B_mean'},numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];
        paras = repmat(paras, nu_runs,1);
        
        % Do Parallel
        infos = nan(size(paras,1),7);
        
        parfor_progress(size(paras,1));   tic
        parfor pp = 1:size(paras,1)
            result = InfoLossGroundTruth([paras(pp,:),'nu_neuron', nu_neuron, 'keepStdMeanRatio',keepStdMeanRatio,'ReLU',ReLU]);
            
            infos(pp,:) = cellfun(@(x)result.(x),{'I_0_optimal','I_0_ilPPC','I_epsi_optimal','I_epsi_ilPPC',...
                                                  'I_epsi_Saturation', 'psycho_Saturation', 'psycho_ilPPC'});
            parfor_progress();
        end
        parfor_progress(0);   toc
        
        % Plot example tuning using the "middle" condition
        % InfoLossGroundTruth([paras(round(end/2),:),'nu_neuron',nu_neuron,'keepStdMeanRatio',keepStdMeanRatio,'ReLU',ReLU, 'plotTuning',1]);
        
        % --- Fetch data ---
        I_0_optimal = reshape(infos(:,1),size(xx,1),size(xx,2),nu_runs);
        I_0_ilPPC = reshape(infos(:,2),size(xx,1),size(xx,2),nu_runs);
        I_epsi_optimal = reshape(infos(:,3),size(xx,1),size(xx,2),nu_runs);
        I_epsi_ilPPC = reshape(infos(:,4),size(xx,1),size(xx,2),nu_runs);
        psycho_Saturation = reshape(infos(:,6),size(xx,1),size(xx,2),nu_runs);
        threshold_epsi = reshape(infos(:,7),size(xx,1),size(xx,2),nu_runs);

        optimality_matrix = mean(I_epsi_ilPPC./I_epsi_optimal,3) * 100;
        optimality_std_matrix = std(I_epsi_ilPPC./I_epsi_optimal,[],3) * 100;
        threshold_matrix = mean(threshold_epsi,3);
        threshold_std_matrix = std(threshold_epsi,[],3);
        
        % =========== Plotting =============
        figure(13); clf; set(gcf,'uni','norm','pos',[0.004        0.29       0.987       0.438]);
        epsi_log_x = [3*log10(epsis(2))-2*log10(epsis(3)) log10(epsis(2:end))]; % Add one entry for epsi = 0
        epsi_log_x_tick = [min(epsi_log_x) linspace(round(epsi_log_x(2)),round(max(epsi_log_x)),5)];
        epsi_log_x_label = [-inf epsi_log_x_tick(2:end)];
        
        if length(mean_baselines) > 1
            
            subplot(1,3,1) % Optimality
            imagesc(epsi_log_x,mean_baselines, optimality_matrix); axis xy; hold on;
            contour(epsi_log_x,mean_baselines, optimality_matrix,'color','k','linew',1.5,'ShowText','on');
            
            set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)
            set(gca,'ytick',0:5:max(mean_baselines))
            colormap(parula); caxis([50 100])
            xlabel('log_{10}(\epsilon_0)');
            ylabel('Mean(Baseline)');
            
            
            subplot(1,3,2)  % Threshold
            imagesc(epsi_log_x,mean_baselines, threshold_matrix); axis xy; hold on;
            contour(epsi_log_x,mean_baselines, threshold_matrix,'color','k','linew',1.5,'ShowText','on');
            
            set(gca,'xtick',epsi_log_x_tick, 'xticklabel', epsi_log_x_label)
            set(gca,'ytick',0:5:max(mean_baselines))
            colormap(hot);
            xlabel('log_{10}(\epsilon_0)');
            ylabel('Mean(Baseline)');
        end
        
        % Threshold and optimality as a function of epsilon
        subplot(1,3,3) 
        to_plot_baseline = mean(mean_baselines); % 20
        [~, ind_baseline] = min(abs(to_plot_baseline - mean_baselines));
        
        optimality_to_plot_baseline = optimality_matrix(ind_baseline,:);
        
        [ax, h1, h2] = plotyy( epsi_log_x, optimality_to_plot_baseline, epsi_log_x, psycho_Saturation(1,:,1) );
        set(h2,'color','b', 'linestyle','--', 'marker', 's', 'linew',2)
        
        hold(ax(1),'on') ; delete(h1);
        errorbar(ax(1), epsi_log_x, optimality_matrix(ind_baseline,:), optimality_std_matrix(ind_baseline,:),'ro-','linew',2);
        set(ax(1),'ytick',50:5:100); ylim(ax(1), [50 102]);
        
        hold(ax(2),'on');
        errorbar(ax(2), epsi_log_x, threshold_matrix(ind_baseline,:), threshold_std_matrix(ind_baseline,:),'ko-','linew',2)
        errorbar(ax(2), epsi_log_x, threshold_matrix(1,:), threshold_std_matrix(1,:),'bo-','linew',2)
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

    case 'Scan general linear 2-D'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        selectScan2D = 1;
        ReLU = 0;
        
        % Automatically keep std ratio
        if selectScan2D > 1  % When not scanning std per se
            keepStdMeanRatio = 1;
        else
            keepStdMeanRatio = 0;
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        scan2D = {    % Variable         Scan             Note
            {    % 1. Baseline Mean vs Baseline Std
            'B_std',  linspace(eps,40,10), 'stdBaseline';
            'B_mean', linspace(eps,40,11), 'meanBaseline';
            }
            {    % 2. Baseline/A vs Drop/A
            '', linspace(eps,1,10), 'meanBaseline/A';
            '',  linspace(eps,1,10), 'meanDrop/A';
            }
            {    % 3. Baseline/A vs A
            '',  linspace(eps,1,10), 'meanBaseline/A';
            'A_mean', linspace(10,100,10), 'Amplitude';
            }
            {    % 4. Drop/A vs A
            '',  linspace(eps,1,10), 'meanDrop/A';
            'A_mean', linspace(10,100,10), 'Amplitude';
            }
            {    % 5. beta vs A 
            'beta_mean', linspace(eps,1,10), 'beta';
            'A_mean', linspace(10,100,10), 'Amplitude';
            }
            {    % 6. Baseline vs Amp
            'B_mean', linspace(5,100,10), 'Baseline';
            'A_mean', linspace(5,100,10), 'Amplitude';
            }
            {    % 7. Baseline/A^2 vs A
            '',  linspace(eps,0.3,10), 'meanBaseline/A^2';
            'A_mean', linspace(10,100,10), 'Amplitude';
            }
            };

        % Default values
        A_mean_default = 47;      A_std_default = 30;
        beta_mean_default = 0.7;   beta_std_default = 0.3;
        B_mean_default = 20;   B_std_default = 20;
        
        
        % Linearize parameters
        xs = scan2D{selectScan2D}{1,2};
        ys = scan2D{selectScan2D}{2,2};
        [xx, yy] = meshgrid(xs, ys);
        
        
        if selectScan2D == 2  % 2. Baseline/A vs Drop /A
            B_means = xx * A_mean_default;
            beta_means = yy ./ xx; % We care about baseline drop ratio in respect to Ai, not Bi.
            beta_means(beta_means > 1) = nan;
            
            paras = [repmat({'B_mean'},numel(B_means),1), mat2cell(B_means(:),ones(numel(B_means),1)),...                    
                     repmat({'beta_mean'},numel(beta_means),1), mat2cell(beta_means(:),ones(numel(beta_means),1)),...
                     ];
                                  
        elseif selectScan2D == 3   % 3. Baseline/A vs A
            B_means = xx .* yy;
            paras = [repmat({'B_mean'},numel(B_means),1), mat2cell(B_means(:),ones(numel(B_means),1)),...
                     repmat(scan2D{selectScan2D}(2,1),numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];
                 
        elseif selectScan2D == 4     % 4. Baseline Drop/A vs A
            beta_means = xx .* yy / B_mean_default;
            beta_means(beta_means > 1) = nan;

            paras = [repmat({'beta_mean'},numel(beta_means),1), mat2cell(beta_means(:),ones(numel(beta_means),1)),...
                     repmat(scan2D{selectScan2D}(2,1),numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];
                 
        elseif selectScan2D == 7   % 7. Baseline/A^2 vs A
            B_means = xx .* yy.^2;
            paras = [repmat({'B_mean'},numel(B_means),1), mat2cell(B_means(:),ones(numel(B_means),1)),...
                     repmat(scan2D{selectScan2D}(2,1),numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];

        else  
            paras = [repmat(scan2D{selectScan2D}(1,1),numel(xx),1), mat2cell(xx(:),ones(numel(xx),1)),...
                     repmat(scan2D{selectScan2D}(2,1),numel(yy),1), mat2cell(yy(:),ones(numel(yy),1))];
        end
        
        paras = repmat(paras, nu_runs,1);
        
        % Do Parallel
        infos = nan(size(paras,1),7);
        
        parfor_progress(size(paras,1));  tic
        parfor pp = 1:size(paras,1)
            
            result = InfoLossGroundTruth([paras(pp,:),'nu_neuron',nu_neuron,'keepStdMeanRatio',keepStdMeanRatio,'ReLU',ReLU]);
            
            infos(pp,:) = cellfun(@(x)result.(x),{'I_0_optimal','I_0_ilPPC','I_epsi_optimal','I_epsi_ilPPC',...
                                                  'I_epsi_Saturation', 'psycho_Saturation', 'psycho_ilPPC'});
            parfor_progress();
        end
        parfor_progress(0);   toc
        
        % Plot example tuning using the "middle" condition
        InfoLossGroundTruth([paras(fix(end/2),:),'nu_neuron',nu_neuron,'keepStdMeanRatio',keepStdMeanRatio,'ReLU',ReLU, 'plotTuning',1]);
        
        % --- Fetch data ---
        I_0_optimal = reshape(infos(:,1),size(xx,1),size(xx,2),nu_runs);
        I_0_ilPPC = reshape(infos(:,2),size(xx,1),size(xx,2),nu_runs);
        I_epsi_optimal = reshape(infos(:,3),size(xx,1),size(xx,2),nu_runs);
        I_epsi_ilPPC = reshape(infos(:,4),size(xx,1),size(xx,2),nu_runs);
        psycho_Saturation = reshape(infos(:,6),size(xx,1),size(xx,2),nu_runs);
        threshold_epsi = reshape(infos(:,7),size(xx,1),size(xx,2),nu_runs);

        optimality_matrix = mean(I_epsi_ilPPC./I_epsi_optimal,3) * 100;
        optimality_std_matrix = std(I_epsi_ilPPC./I_epsi_optimal,[],3) * 100;
        threshold_matrix = mean(threshold_epsi,3);
        threshold_std_matrix = std(threshold_epsi,[],3);
        
        % ========= Plotting =========
        figure(13); clf;  set(gcf,'uni','norm','pos',[0.118       0.333       0.726       0.453]);
        
        subplot(1,2,1)
        imagesc(xs, ys, optimality_matrix); axis xy; hold on;
        contour(xs, ys, optimality_matrix,'color','k','linew',1.5,'ShowText','on');
        
        %     set(gca,'ytick',0:5:max(mean_baselines))
        colormap(parula); % caxis([50 100])
        xlabel(scan2D{selectScan2D}{1,3});
        ylabel(scan2D{selectScan2D}{2,3});
        
        subplot(1,2,2)
        imagesc(xs, ys, threshold_matrix); axis xy; hold on;
        contour(xs, ys, threshold_matrix,'color','k','linew',1.5,'ShowText','on');
        
        colormap jet
        %     set(gca,'ytick',0:5:max(mean_baselines))
        % caxis([50 100])
        xlabel(scan2D{selectScan2D}{1,3});
        ylabel(scan2D{selectScan2D}{2,3});
        
        
        if ION_cluster
            saveas(13,['./4_Optimality_matrix_GroundTruth_General2D_selectScan2D=' num2str(selectScan2D)]);
        end
end


%% Plotting

% ====== Plot momentary information ========
% 
% if length(nu_neurons) > 1  % 1-D Scan number of neurons
%     % %{
%     figure(10); plot(time,I_0_ts,'r','linew',2)
%     hold on; plot(time,I_epsi_ts,'b','linew',2);
%     hold on; plot(time,speed_t*max(I_0_ts),'k--');
%     %}
%     
% end
 
%%
% ========== Tuning properties =========
% InfoLossGroundTruth({'plotTuning',1});


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

