%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute average log odds in the model as a function of actual perf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
x = [200:1200];
load monkey_data
coh_monk = coherence;
clear coherence;

nu_trial=50;
constant_RT=290+30;

load log_odds
load pdf_info

load data0-10
load exp_data0-10
RT_mat(:,1)= RT(1:nu_trial);
conf_mat(:,1) = conf_level(1:nu_trial);
conf_proba_mat(:,1) = conf_level_proba(1:nu_trial);
coh(1) = 1;
performance(1)= perf;
RT_cor(1) = mean_RT;
RT_incor(1)=mean_RT_incor;
neuron_traj_mat(:,1) = neuron_traj;
antineuron_traj_mat(:,1) = antineuron_traj;

load data3-10
load exp_data32-10
RT_mat(:,2)= RT(1:nu_trial);
conf_mat(:,2) = conf_level(1:nu_trial);
conf_proba_mat(:,2) = conf_level_proba(1:nu_trial);
coh(2) = coherence;
performance(2)= perf;
RT_cor(2) = mean_RT;
RT_incor(2)=mean_RT_incor;
neuron_traj_mat(:,2) = neuron_traj;
antineuron_traj_mat(:,2) = antineuron_traj;

load data6-10
load exp_data64-10
RT_mat(:,3)= RT(1:nu_trial);
conf_mat(:,3) = conf_level(1:nu_trial);
conf_proba_mat(:,3) = conf_level_proba(1:nu_trial);
coh(3) = coherence;
performance(3)= perf;
RT_cor(3) = mean_RT;
RT_incor(3)=mean_RT_incor;
neuron_traj_mat(:,3) = neuron_traj;
antineuron_traj_mat(:,3) = antineuron_traj;

load data12-10
load exp_data128-10
RT_mat(:,4)= RT(1:nu_trial);
conf_mat(:,4) = conf_level(1:nu_trial);
conf_proba_mat(:,4) = conf_level_proba(1:nu_trial);
coh(4) = coherence;
performance(4)= perf;
RT_cor(4) = mean_RT;
RT_incor(4)=mean_RT_incor;
neuron_traj_mat(:,4) = neuron_traj;
antineuron_traj_mat(:,4) = antineuron_traj;

load data25-10
load exp_data256-10
RT_mat(:,5)= RT(1:nu_trial);
conf_mat(:,5) = conf_level(1:nu_trial);
conf_proba_mat(:,5) = conf_level_proba(1:nu_trial);
coh(5) = coherence;
performance(5)= perf;
RT_cor(5) = mean_RT;
RT_incor(5)=mean_RT_incor;
neuron_traj_mat(:,5) = neuron_traj;
antineuron_traj_mat(:,5) = antineuron_traj;

load data51-10
load exp_data512-10
RT_mat(:,6)= RT(1:nu_trial);
conf_mat(:,6) = conf_level(1:nu_trial);
conf_proba_mat(:,6) = conf_level_proba(1:nu_trial);
coh(6) = coherence;
performance(6)= perf;
RT_cor(6) = mean_RT;
RT_incor(6)=mean_RT_incor;
neuron_traj_mat(:,6) = neuron_traj;
antineuron_traj_mat(:,6) = antineuron_traj;



[aux_max time_thres]=max(neuron_traj_mat>58);
for k=1:6
   diff_ac(k) = neuron_traj_mat(time_thres(k),k)-antineuron_traj_mat(time_thres(k),k);
end

subplot(221)
plot(coh,diff_ac)

subplot(222)
hist(conf_proba_mat)

subplot(223)
plot(sum(conf_proba_mat.*(conf_proba_mat>0))./sum(conf_proba_mat>0),'s-')

subplot(224)
plot(neuron_traj_mat)
hold on
plot(antineuron_traj_mat,'--')
hold off
xlabel('Reaction time (ms)');
ylabel('Firing rate (Hz)');
axis([0 500 0 58]);



