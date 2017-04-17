clear
x = [200:1200];
load monkey_data
coh_monk = coherence;
clear coherence;

constant_RT=50;

load log_odds
load pdf_info

load data0
load exp_data0

nu_trial = size(RT)/2;

RT_mat(:,1)= RT(1:nu_trial);
conf_mat(:,1) = conf_level(1:nu_trial);
coh(1) = 1;
performance(1)= perf;
RT_cor(1) = mean_RT;
RT_incor(1)=mean_RT_incor;
neuron_traj_mat(:,1) = neuron_traj;
antineuron_traj_mat(:,1) = antineuron_traj;

load data3
load exp_data32
RT_mat(:,2)= RT(1:nu_trial);
conf_mat(:,2) = conf_level(1:nu_trial);
coh(2) = coherence;
performance(2)= perf;
RT_cor(2) = mean_RT;
RT_incor(2)=mean_RT_incor;
neuron_traj_mat(:,2) = neuron_traj;
antineuron_traj_mat(:,2) = antineuron_traj;

load data6
load exp_data64
RT_mat(:,3)= RT(1:nu_trial);
conf_mat(:,3) = conf_level(1:nu_trial);
coh(3) = coherence;
performance(3)= perf;
RT_cor(3) = mean_RT;
RT_incor(3)=mean_RT_incor;
neuron_traj_mat(:,3) = neuron_traj;
antineuron_traj_mat(:,3) = antineuron_traj;

load data12
load exp_data128
RT_mat(:,4)= RT(1:nu_trial);
conf_mat(:,4) = conf_level(1:nu_trial);
coh(4) = coherence;
performance(4)= perf;
RT_cor(4) = mean_RT;
RT_incor(4)=mean_RT_incor;
neuron_traj_mat(:,4) = neuron_traj;
antineuron_traj_mat(:,4) = antineuron_traj;

load data25
load exp_data256
RT_mat(:,5)= RT(1:nu_trial);
conf_mat(:,5) = conf_level(1:nu_trial);
coh(5) = coherence;
performance(5)= perf;
RT_cor(5) = mean_RT;
RT_incor(5)=mean_RT_incor;
neuron_traj_mat(:,5) = neuron_traj;
antineuron_traj_mat(:,5) = antineuron_traj;

load data51
load exp_data512
RT_mat(:,6)= RT(1:nu_trial);
conf_mat(:,6) = conf_level(1:nu_trial);
coh(6) = coherence;
performance(6)= perf;
RT_cor(6) = mean_RT;
RT_incor(6)=mean_RT_incor;
neuron_traj_mat(:,6) = neuron_traj;
antineuron_traj_mat(:,6) = antineuron_traj;


conf_mat = abs(conf_mat);
RT_mat= RT_mat+constant_RT;
RT_cor=RT_cor+constant_RT;
RT_incor=RT_incor+constant_RT;

for j=1:6
   aux_x = RT_mat(:,j)-mean(RT_mat(:,j));
   aux_y = conf_mat(:,j)-mean(conf_mat(:,j));
   a(j) = (aux_x'*aux_y)/(aux_x'*aux_x);
   b(j) = mean(conf_mat(:,j))-a(j)*mean(RT_mat(:,j));
   lin_reg(j,:) = x*a(j)+b(j);
end

RT_cor=RT_cor;
RT_incor=RT_incor;
coh_maz = [1 3.2 6.4 12.8 25.6 51.2];
RT_maz =        [820     800    750     670     550     420];
RT_incor_maz =  [860    860     860     750 ];
performance_maz = [50 64 77 94 99 100]/100;


subplot(221)
semilogx(coh,performance,'b*-')
hold on
semilogx(coh_maz,performance_maz,'ro-')
semilogx(coh_monk,choice2_perf,'go-');
semilogx(coh_monk,choice4_perf,'co-');
hold off
xlabel('log Coherence');
ylabel('Probability correct');
axis([1 100 0.0 1]);


subplot(222)
semilogx(coh,RT_cor,'*-b')

hold on
semilogx([3 6 12 24 ],RT_incor(2:5),'b*--')
semilogx(coh,RT_maz,'ro-')
semilogx(coh_monk,choice4_RT,'co-')
semilogx(coh_monk,choice2_RT,'go-')
semilogx([3.2 6.4 12.8 25.6 ], RT_incor_maz,'ro--')
hold off
xlabel('log Coherence');
ylabel('Reaction time (ms)');
%axis([1 53 400 2000]);

subplot(223)
plot(RT_mat,conf_mat,'.')
hold on
% plot(x,lin_reg')
hold off
xlabel('Reaction time (ms)');
ylabel('Confidence (log odds)');
axis tight

subplot(224)
plot(neuron_traj_mat)
hold on
plot(antineuron_traj_mat,'--')
hold off
xlabel('Reaction time (ms)');
ylabel('Firing rate (Hz)');
axis([0 130 35 45]);

