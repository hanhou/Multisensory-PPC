%%
population_data = reshape(aver_population_act,N_lip,[])';
correct_ans = repmat([-1 -1 -1 -1 1 1 1 1 1]',size(population_data,1)/9,1);

% population_data = rate_lip_toy_to_train; % Trials * N_lip 
% correct_ans = correct_ans_train;

% population_data = rate_lip_to_train; % Trials * N_lip 
% correct_ans = correct_ans_train;

stim_types = reshape(repmat([1 2 3],length(correct_ans)/3,1),[],1);

[coeff,score,latent] = pca(population_data);  coeff(:,1) = -coeff(:,1); score(:,1) = -score(:,1);
figure(2139);   clf
set(gcf,'uni','norm','pos',[0.018       0.212       0.912       0.403]);
subplot(2,3,1); 
plot(prefs_lip,coeff(:,1:2),'linew',2);
legend('PC1','PC2');

subplot(2,3,[2 5]); 
colors_char = {'b','r','g'};

for ss = 1:3
    plot(score(stim_types==ss & correct_ans < 0,1),score(stim_types==ss & correct_ans < 0,2),'o','color',colors_char{ss}); 
    hold on; 
    plot(score(stim_types==ss & correct_ans > 0,1),score(stim_types==ss & correct_ans > 0,2),...
        'o','color',colors_char{ss},'markerfacecol',colors_char{ss}); 
end

axis tight;

C = 1e-5;
temp_m = svmtrain(population_data,correct_ans,'box',C,'autoscale',1,'tol',1e-7);
temp_weight = -temp_m.SupportVectors'*temp_m.Alpha;
temp_weight = temp_weight./temp_m.ScaleData.scaleFactor';

% Projection on PCA space
p = temp_weight'*coeff(:,1:2);

hold on; 
%  plot(xlim,xlim * p(2)/p(1),'k-','linew',2); % SVM weight
plot(xlim,xlim * -p(1)/p(2),'k-','linew',2); % SVM seperate
xlabel('PC1'); ylabel('PC2');

choices = svmclassify(temp_m,population_data);
title(sprintf('C = %g, correct rate = %g',C, sum(choices==correct_ans)/length(correct_ans)));

subplot(2,3,4); hold on;  plot(prefs_lip,temp_weight,'k','linew',2)
legend('SVM');

% Hist of projection
hhh = subplot(2,3,[3 6]);
proj_minus_group = population_data(correct_ans < 0,:) * temp_weight;
proj_plus_group = population_data(correct_ans > 0,:) * temp_weight;
[n,x] = hist(hhh,proj_minus_group,20); 
bar(x,n,'k'); hold on;
[n,x] = hist(hhh,proj_plus_group,20);
bar(x,n,'r'); 
auc = rocN(proj_plus_group,proj_minus_group);
title(['auc = ' num2str(auc)]);

%% Simplest demo
% load fisheriris
% xdata = meas(51:end,3:4);
% group = species(51:end);

g1 = [(randn(100,1)+1)*0.2 (randn(100,1)+1)];
g2 = [(randn(100,1)+5)*0.2 (randn(100,1)+5)];

xdata = [g1;g2];
group = [zeros(100,1); ones(100,1)];
 
% xdata = [1 2;2 1;4 3; 3 4;5 6];
% group = [-1 -1 1 1 1]';

figure(2226); 
svmStruct = svmtrain(xdata,group,'ShowPlot',true,'box',1e-3,'autoscale',1,'Method','QP');

ww = svmStruct.SupportVectors'*svmStruct.Alpha;
% ww = ww./svmStruct.ScaleData.scaleFactor';
ww_correct_reweight = ww.*svmStruct.ScaleData.scaleFactor';
%sep = [-ww(2) ww(1)];
hold on; plot(xlim,mean(xdata(:,2))+(xlim-mean(xdata(:,1))) * ww_correct_reweight(2)/ww_correct_reweight(1),'c--');

ww_wrong_reweight = ww./svmStruct.ScaleData.scaleFactor'; % This reweight method = 'autoscale = off' with very small box constraint
plot(xlim,mean(xdata(:,2))+(xlim-mean(xdata(:,1))) * ww_wrong_reweight(2)/ww_wrong_reweight(1),'g');
axis equal

svmStruct = svmtrain(xdata,group,'ShowPlot',true,'box',1e-5,'autoscale',0,'tol',1e-7);
ww = svmStruct.SupportVectors'*svmStruct.Alpha;
hold on; plot(xlim,mean(xdata(:,2))+(xlim-mean(xdata(:,1))) * ww(2)/ww(1),'m-');

