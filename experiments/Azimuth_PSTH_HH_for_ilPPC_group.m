%% Importing data

[filename, pathname] = uigetfile( 'Z:\Data\Tempo\Batch\*.mat','MultiSelect','on');

nFile = length(filename);

load([pathname filename{1}]);
ROI = result.ROI;

group_tuning = nan([size(result.mean_firing_matrix_wrap) nFile]);
group_tuning_norm = nan([size(result.mean_firing_matrix_wrap) nFile]);
group_p_value = nan(3,nFile);

progressbar();
for nn = 1:nFile
    load([pathname filename{nn}]);
    group_tuning(:,:,:,nn) = result.mean_firing_matrix_wrap;
    group_tuning_norm(:,:,:,nn) = result.mean_firing_matrix_wrap./...
        repmat(result.mean_firing_matrix_wrap(ceil(size(result.mean_firing_matrix_wrap,1)/2),end,:),...
        size(result.mean_firing_matrix_wrap,1), size(result.mean_firing_matrix_wrap,2));
    group_p_value(:,nn) = result.p_value(:,1);
    progressbar(nn/nFile);
end

for k = 1:3
    this_sign = find(group_p_value(k,:) < 0.05);
    n_sign(k) = length(this_sign);
    mean_group_tuning = mean(group_tuning(:,:,:,this_sign),4);
    se_group_tuning = std(group_tuning(:,:,:,this_sign),[],4)/sqrt(n_sign(k));
    
    mean_group_tuning_norm = mean(group_tuning_norm(:,:,:,this_sign),4);
    se_group_tuning_norm = std(group_tuning_norm(:,:,:,this_sign),[],4)/sqrt(n_sign(k));
    
    % Calculate info loss using Alex's code [See MST_Heter_Info_Loss_GroundTruth for full accounts of this issue]
    %{
    progressbar();
    for cc = 1:length(this_sign)
        ind = this_sign(cc);
        
        % May not be correct
        real_baseline =  min(min(group_tuning(:,:,k,ind),[],1)); % That can be moved by ReLU
        baseline_tc = mean(mean(group_tuning([1:2 8:9],:,k,ind),2)); % All time average of the anti-pref firings
        amp_tc = max(mean(group_tuning(4:6,:,k,ind),1)); % Average of firings around pref direction, then max over time
        group_info_loss{k}(cc) = info_loss(baseline_tc-real_baseline,amp_tc-real_baseline);
        progressbar(cc/length(this_sign));
    end
    %}
end




%% Plotting
%%%%%%%%%%%%%%%%%%%%%
time_range = 3:length(ROI)-3;
%%%%%%%%%%%%%%%%%%%%%

default_color_order = get(0,'defaultaxescolororder');
color_order = colormap(jet);
color_order = color_order(round(linspace(1,size(color_order,1),length(time_range))),:);

set(0,'defaultaxescolororder',color_order);

% Raw
figure(1610); clf;    set(gcf,'uni','norm','pos',[0.023        0.34       0.855       0.365]);
xs = [-180:45:180]';

for k = 1:3
    subplot(1,3,k);
    errorbar(repmat(xs,1,length(time_range)),mean_group_tuning(:,time_range,k),se_group_tuning(:,time_range,k),'linew',2);
    hold on;
    title(['n=' num2str(n_sign(k))]);
    xlim([-200 200]);
    set(gca,'xtick',[-180:90:180]);
    ylim([0 60]);
    %     plot(mean_group_tuning(:,end,k),'k','linew',2);
end
legend(num2str(ROI(time_range,:)))

% Norm
figure(1619); clf;    set(gcf,'uni','norm','pos',[0.023        0.34       0.855       0.365]);
for k = 1:3
    subplot(1,3,k);
    errorbar(repmat(xs,1,length(time_range)),mean_group_tuning_norm(:,time_range,k),se_group_tuning_norm(:,time_range,k),'linew',2);
    hold on;
    title(['n=' num2str(n_sign(k))]);
    xlim([-200 200]);
    set(gca,'xtick',[-180:90:180]);
    %     plot(mean_group_tuning_norm(:,end,k),'k','linew',2);
end
legend(num2str(ROI(time_range,:)))

% Info_loss
%{
figure(1620); clf;    set(gcf,'uni','norm','pos',[0.023        0.34       0.855       0.365]);
for k = 1:3
    subplot(1,3,k);
    hist(group_info_loss{k},20);
    title(['n=' num2str(n_sign(k))]);
    xlim([0 100]);
end
%}


%% Gaussian
paras = {
          2   0.13  6    'Gu 2006 2010 tuning (combined)'
%           2   0.2   4.5  'Gu 2012 heading (combined)'
%           1   0.2   2    'HeTao heading (visual only)'
%           2 0.18 6 'Polo training'
          1.5   0.2  3.5   'Polo heading (original)'
%           1.5   0.1  4  'Polo heading (combined)'
%           2.0   0.2  4.0  'Polo heading (combined)'
%           1.5   0.14  3.8  'Minimoog'
%           1.5   0.1  4.9  'Change sigma 1'
%           1.5   0.25  2.3  'Change sigma 2'
};

i=1;
maxDuration =  max(cell2mat(paras(:,1)));
duration = paras{i,1};
step=0.005;
t = -0.4:step:duration+0.4;

% delay = 0.14;
delay = 0;

ampl = paras{i,2};
num_sigs = paras{i,3};

pos=ampl*0.5*(erf((t-duration/2-delay)/((duration/2/num_sigs)*sqrt(2))) + 1);  % Note that sigma = duration/2/num_of_sigma !!
veloc = diff(pos)/step;

figure(2237); % subplot(1,2,1); 

plot(t(1:length(t)-1),veloc,'k','linew',2);

i = 1;
for rr = time_range
    hold on;
    plot(ROI(rr,:)/1000,[-0.1 -0.1],'color',color_order(i,:),'linew',9)
    i = i + 1;
    % aver_speed(rr) = mean(veloc(ROI(rr,1)/1000<t & t< ROI(rr,2)/1000));
end

%{
%======== Simple Simulation of Gain modulation
% color_order = colormap(jet);
% color_order = color_order(round(linspace(1,size(color_order,1),size(ROI,1))),:);

set(0,'defaultaxescolororder',color_order);

% -- Const. Baseline + Pure gain
subplot(1,2,2); hold on; set(gcf,'uni','norm','pos',[0.131       0.424       0.609       0.363]);
x=-pi:0.1:pi;

F = @(x,n)(exp(n.*x)-1)./n;
spatial_tuning = F(cos(x),1);

manual_gains = [0 0 0 0 0.05 0.4 0.8 0.6 0.4 0.3 0.3 0.3 0.3];

for k=1:length(ROI)-1
    plot(x/pi*180,0.6 + manual_gains(k)*spatial_tuning,'linew',2);
end
set(gca,'xtick',[-180:90:180]);

set(0,'defaultaxescolororder','default');

%}

