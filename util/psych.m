clear
clear global
tic

%nu_trial = 500
%save nu_trial nu_trial

%nu_coherence=6;
% 
% for j=1:nu_coherence
%     parameters
%     coherence = (2^(j+3)/10);
%     %% set coherence to 0 if equal to 1
%     coherence=coherence*(1-(coherence==1.6));
%     fprintf('coherence %4.1f\n',coherence);  
%     
%     lip2
%     data_name = num2str(coherence*10,3);
%     data_name = strcat('2data',data_name)
%     save(data_name,'coherence');
%     save(data_name,'RT conf_level','conf_level_proba','perf_trial',...
%         'perf','mean_RT','mean_RT_incor','coherence','w_ole_outx',...
%         'info_in','info_train_in','info_outx','info_train_outx');
%     exp_data_name= num2str(coherence*10,3);
%     = strcat('00data',data_name) = strcat('2exp_data',exp_data_name)
%     save(exp_data_name,'spike_time50','spike_time100','RT','coherence',...
%         'perf_trial','proba_out50','proba_out100','spikes_out50',...
%         'spikes_out100','neuron_traj','antineuron_traj');
% end



clear
disp('coherence 0');
coherence = 0;
lip2
save data0-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data0-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 3.2');
coherence = 3.2;
lip2
save data3-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data32-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 6.4');
coherence = 6.4;
lip2
save data6-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data64-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 12.8');
coherence = 12.8;
lip2
save data12-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data128-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 25.6');
coherence = 25.6;
lip2
save data25-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data256-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 51.2');
coherence = 51.2;
lip2
save data51-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data512-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

clear
disp('coherence 76.8');
coherence = 76.8;
lip2
save data76-4-2 RT conf_level conf_level_proba perf_trial perf mean_RT mean_RT_incor coherence w_ole_outx info_in info_train_in info_outx info_train_outx
save exp_data768-4-2 spike_time50 spike_time100 RT coherence perf_trial proba_out50 proba_out75 proba_out100 spikes_out50 spikes_out75 spikes_out100 neuron_traj antineuron_traj

beep
toc

plot_res