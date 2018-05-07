function result = generateCorrelatedPoisson(rho,f0,dt)
% Verify the method of generating correlated Poisson spikes
% HH20180503isdeployed

nu_neuron = 40;
nu_rep = 200;
theta_pref = linspace(-pi,pi,nu_neuron)';

if nargin == 0
    f0 = 50; % in Hz, Constant firing rate to avoid the nonlinearity in p
    rho = 0.2;
    dt = 1e-3; % s
%     dt = 1/(2*f0); % Best dt ???
end

f = f0 * ones(nu_neuron,1); % Uniform f
% f = f0 * exp(2*(cos(theta_pref)-1)); % With tuning


time_max = dt * 100; % Not important how long the time is
times = 0:dt:time_max;

% C = rand(nu_neuron) ; % Arbitrary C
% -- Correlations
k = 2;  

% Desired correlation coefficient (Moreno 2014 NN)
C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); 

desired_SIGMA = C .* sqrt(f * f');
W = real(sqrtm(desired_SIGMA  - diag(f)));

spikes = nan(nu_neuron,length(times),nu_rep);

% parfor_progress(length(times));
parfor tt = 1:length(times)
    
    prob =  dt * f * ones(1,nu_rep) + W * normrnd(0,1,nu_neuron,nu_rep) * sqrt(dt); %%% THIS IS WHAT I FOUND %%%
    %     prob =  dt * f * ones(1,nu_rep) + W * normrnd(0,1,nu_neuron,nu_rep); %%% Beck 2008 %%%
    
    spikes(:,tt,:) = rand(size(prob)) < prob ; % Poisson
%     parfor_progress;
    
end

%% Compare real and desired covariance matrix
bootstramp = 1; 
each_rep = nu_rep/1;

all_activity = squeeze(sum(spikes(:,:,:),2));

fano_matrix_boot = zeros(nu_neuron);
corrrcoef_matrix_boot = zeros(nu_neuron);

% parfor_progress(bootstramp);

for bb = 1:bootstramp
    temp_activity = all_activity(:,randperm(nu_rep,each_rep));
    
%   cov_matrix_this = cov(temp_activity');  % Cov_rate = cov_matrix / dt should be desired_SIGMA when fano = 1.
    
%     fano_matrix = cov_matrix ./ sqrt(mean(temp_activity,2) * mean(temp_activity,2)');  % Var / Mean
%     fano_matrix_boot = fano_matrix_boot + fano_matrix;
    
    cor_matrix_this = corrcoef(temp_activity');  
    corrrcoef_matrix_boot = corrrcoef_matrix_boot + cor_matrix_this;
    
%     parfor_progress;
end

%% == Compare results ==
% Correlation matrices
real_corrcoef_matrix = corrrcoef_matrix_boot / bootstramp;
desired_corrcoef_matrix = C; 
pred_corrcoef_matrix = (W * W' + diag(f)) ./ sqrt(f * f');  % May not be the same as desired because non-zero residuals from sqrtm()
                                                            % Also note that its W*W', not W*W.

% Average correlations
real_aver_corrcoef = circAverage(real_corrcoef_matrix);
desired_aver_corrcoef = circAverage(desired_corrcoef_matrix);
pred_aver_corrcoef = circAverage(pred_corrcoef_matrix);

% Mean corrcoef. Follow Beck, 2008
aux_mask = ones(nu_neuron,1) * (1:nu_neuron);
mask = ((1-cos(abs(aux_mask-aux_mask')/nu_neuron*2*pi))/2)<.5;
mask = mask.*(1-eye(nu_neuron));

mean_corr_less_than_90_actual = nanmean(real_corrcoef_matrix(logical(mask)));
mean_corr_less_than_90_desired = nanmean(desired_corrcoef_matrix(logical(mask)));

% Distortion of rate and fano factor
real_firing_rate = mean(mean(spikes(:,:,:),3),2)/dt;
real_fano = var(all_activity,[],2)./mean(all_activity,2);

% Save errors
result.err_mean_corr = mean(real_aver_corrcoef - desired_aver_corrcoef) ; % Error of mean of all elements of corrcoef matrix (underestimated the error)
result.err_max_corr = mean(real_aver_corrcoef(round(end/2)-1:round(end/2)+1) - desired_aver_corrcoef(round(end/2)-1:round(end/2)+1)) ; % Error of the Corr_(i,i+1)s (directly comparable with rho)
result.err_firing_rate = mean(real_firing_rate - f);
result.err_fano = mean(real_fano - 1);
result.mean_corr_less_than_90 = mean_corr_less_than_90_actual;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% == Plotting ==
if nargin == 0
    
    figure(1742); clf;
    set(gcf,'uni','norm','pos',[0.04       0.125       0.925       0.684]);
        
    %  subplot(1,4,1); plot(fano_matrix);
    subplot(2,3,2);
    plot(theta_pref/pi*180,f,'k-o');  hold on; 
    plot(theta_pref/pi*180,real_firing_rate,'r-o'); 
    xlabel('Pref'); set(gca,'xtick',-180:90:180);
    title('Firing rate'); legend({'Desired','Actual'},'location','best');
    ylim([0 max(ylim)*1.1]);
    
    subplot(2,3,3);
    plot(theta_pref/pi*180,real_fano,'or-'); title('Fano factor');
    hold on; plot(xlim,[1 1],'--k'); ylim([0 max(ylim)*1.1]);
    xlabel('Pref'); set(gca,'xtick',-180:90:180);
    
    % imagesc(fano_matrix_boot / bootstramp - eye(size(C))); colorbar(); title('RealCov/sqrt(RealMean_{i*j}) - I');
    % text(0,0,sprintf('meanDiag = %.2g',mean(diag(fano_matrix_boot / bootstramp - eye(size(C))))));
    % caxis([0 .5])
    
    h1 = subplot(2,3,6);
    imagesc(theta_pref/pi*180,theta_pref/pi*180,real_corrcoef_matrix - eye(size(C))); colorbar();
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
   
    title(sprintf('Actual CorrCoef - I\ncor90 = %.4g',mean_corr_less_than_90_actual));
    axis square
    % text(0,0,sprintf('meanDiag = %.2g',mean(diag(real_corrcoef_matrix - eye(size(C))))));
    % caxis([0 .5])
    
    h2 = subplot(2,3,5);
    imagesc(theta_pref/pi*180,theta_pref/pi*180, desired_corrcoef_matrix - eye(size(C))) ;
    colorbar(); title(sprintf('Desired CorrCoef - I\ncor90 = %.4g', mean_corr_less_than_90_desired));
    set(gca,'xtick',-180:90:180,'ytick',-180:90:180);
    axis square
    
%     subplot(2,3,6)
%     imagesc( pred_corrcoef_matrix ) ;
%     colorbar(); title('Predicted cov matrix / sqrt(f*f'') - I');
%     colormap(jet)
    
    subplot(2,3,4);
%     plot([real_aver_corrcoef' desired_aver_corrcoef' pred_aver_corrcoef'],'linew',2); legend({'real'; 'desired'; 'pred'});
    plot(theta_pref/pi*180, desired_aver_corrcoef,'-k'); hold on;
    plot(theta_pref/pi*180, real_aver_corrcoef,'-r','linew',3); hold on;
    xlabel('\DeltaPref'); title('CorrCoef');  set(gca,'xtick',-180:90:180);
    legend({'Desired'; 'Actual'});
   
    % caxis([0 .5])
    
    set(findall(gcf,'type','line'),'linew',2);
    SetFigure(13);
    
    h = axes('Position',[0.14 0.717 0.168 0.197]);
    set(text(h,0,0,sprintf('\\deltat = %0.0e s\nf_0 = %g Hz\n\\rho = %g',dt, f0, rho)),'FontSize',20);
    axis off;
    
    
    %%{
    ylims = cell2mat(get([h1 h2],'clim'));
    ylnew = [min(ylims(:,1)) max(ylims(:,2))];
    set(findall(gcf,'type','axes'),'clim',ylnew);
    colormap(jet);
    %}    

    file_name = sprintf('./Test_correlated_poisson');
    saveas(gcf,'Test_corr_Poisson.fig');
end
