% Varify the method of generating correlated Poisson spikes

time_max = 1; % s
nu_neuron = 20;
nu_rep = 200;

dt = 1e-3; % s
times = 0:dt:time_max;

% C = rand(nu_neuron) ; % Arbitrary C
% -- Correlations
rho = 0.5;  k = 2;  
theta_pref = linspace(-pi,pi,nu_neuron)';
C = (1-rho)*eye(nu_neuron) + rho*exp(k*(cos(theta_pref * ones(1,nu_neuron) - ones(nu_neuron,1) * theta_pref')-1)); % Correlation coefficient (Moreno 2014 NN)

f = 50 * ones(nu_neuron,1); % in Hz, Constant firing rate to avoid the nonlinearity in p

desired_SIGMA = C .* sqrt(f * f');
W = real(sqrtm(desired_SIGMA  - diag(f)));

spikes = nan(nu_neuron,length(times),nu_rep);

parfor_progress(length(times));

for tt = 1:length(times)
    
    prob =  dt * f * ones(1,nu_rep) + W * normrnd(0,1,nu_neuron,nu_rep) * sqrt(dt); %%% THIS IS WHAT I FOUND %%%
    %     prob =  dt * f * ones(1,nu_rep) + W * normrnd(0,1,nu_neuron,nu_rep); %%% Beck 2008 %%%
    
    spikes(:,tt,:) = rand(size(prob)) < prob ; % Poisson
    parfor_progress;
    
end

%%  Real cov
bootstramp = 1; 
each_rep = nu_rep/1;

all_activity = squeeze(sum(spikes(:,:,:),2));

fano_matrix_boot = zeros(nu_neuron);
cor_matrix_boot = zeros(nu_neuron);

parfor_progress(bootstramp);

for bb = 1:bootstramp
    temp_activity = all_activity(:,randperm(nu_rep,each_rep));
    
    cov_matrix = cov(temp_activity');
    
%     fano_matrix = cov_matrix ./ sqrt(mean(temp_activity,2) * mean(temp_activity,2)');  % Var / Mean
%     fano_matrix_boot = fano_matrix_boot + fano_matrix;
    
    cor_matrix = corrcoef(temp_activity');  
    cor_matrix_boot = cor_matrix_boot + cor_matrix;
    
    parfor_progress;
end


figure(1742); clf;
set(gcf,'uni','norm','pos',[0.04       0.125       0.925       0.684]);

%  subplot(1,4,1); plot(fano_matrix);
subplot(2,3,1); 
aver_rate = mean(mean(spikes(:,:,:),3),2)/dt;
plot(aver_rate,'b-o'); hold on; plot(f,'k-o');  title('RealMean vs f');
ylim([0 max(ylim)*1.1]);

subplot(2,3,2); 
fano = var(all_activity,[],2)./mean(all_activity,2);
plot(fano,'o-'); title('Fano factor'); hold on; plot(xlim,[1 1],'--k'); ylim([0 1.1]);

% imagesc(fano_matrix_boot / bootstramp - eye(size(C))); colorbar(); title('RealCov/sqrt(RealMean_{i*j}) - I');
% text(0,0,sprintf('meanDiag = %.2g',mean(diag(fano_matrix_boot / bootstramp - eye(size(C))))));
% caxis([0 .5])

subplot(2,3,3)
imagesc(cor_matrix_boot / bootstramp - eye(size(C))); colorbar();  title('Real corrcoef - I');
text(0,0,sprintf('meanDiag = %.2g',mean(diag(cor_matrix_boot / bootstramp - eye(size(C))))));
% caxis([0 .5])


subplot(2,3,5)
imagesc( desired_SIGMA ./ sqrt(f * f') - eye(size(C))) ; 
colorbar(); title('Desired cov matrix / sqrt(f*f'') - I');

subplot(2,3,6)
pred_cov = (W * W + diag(f)) ./ sqrt(f * f') - eye(size(C)); 
imagesc( pred_cov ) ; 
colorbar(); title('Predicted cov matrix / sqrt(f*f'') - I');
colormap(jet)

subplot(2,3,4);
real_aver = circ_Aver(cor_matrix_boot / bootstramp);
desired_aver = circ_Aver(desired_SIGMA ./ sqrt(f * f'));
pred_aver = circ_Aver((W * W + diag(f)) ./ sqrt(f * f'));
plot([real_aver' desired_aver' pred_aver'],'linew',2); legend(['real' 'desired' 'pred']);

file_name = sprintf('./Test_correlated_poisson');
saveas(gcf,'Test_corr_Poisson.fig');
% caxis([0 .5])

%{
ylims = cell2mat(get(findall(gcf,'type','axes'),'clim'));
ylnew = [min(ylims(:,1)) max(ylims(:,2))];
set(findall(gcf,'type','axes'),'clim',ylnew);
%}
