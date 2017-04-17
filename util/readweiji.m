clear

load all_exp_data2

nu_coherence = 6;

[nu_trial junk] = size(A);

nu_trial = nu_trial/12;

sp_count50=0;
sp_count100=0;


for j=1:nu_trial
 sp_count50 = sp_count50 + length(spike_time{j+nu_trial*10}); 
 sp_count100 = sp_count100 + length(spike_time{j+nu_trial*11}); 
end

sp_count50 = sp_count50 /nu_trial;
sp_count100 = sp_count100 /nu_trial;


spikes_out50 = zeros(nu_trial,1500); 
spikes_out100 = zeros(nu_trial,1500); 

for k=1:nu_trial
   nu_spikes = length(spike_time{k+nu_trial*10}); 
   for l=1:nu_spikes
       spikes_out50(k,spike_time{k+nu_trial*10}(l))=1;
   end   
   nu_spikes = length(spike_time{k+nu_trial*11}); 
   for l=1:nu_spikes
       spikes_out100(k,spike_time{k+nu_trial*11}(l))=1;
   end
end



plot(sum(spikes_out50),'r')
hold on
plot(sum(spikes_out100))
hold off